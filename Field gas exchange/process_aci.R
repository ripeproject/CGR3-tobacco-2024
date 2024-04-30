###
### PRELIMINARIES:
### Loading packages, defining constants, creating helping functions, etc.
###

# Load required packages
library(PhotoGEA)
library(lattice)

# Check to make sure we have the correct version of PhotoGEA
if (packageVersion('PhotoGEA') != '0.10.0') {
    stop(
        'The `process_aci.R` script requires PhotoGEA version 0.10.0. ',
        'See README.md for installation instructions.'
    )
}

# Clear the workspace
rm(list=ls())

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- 'event'
REP_COLUMN_NAME <- 'replicate'
PLOT_COLUMN_NAME <- 'plot'
PHIPS2_COLUMN_NAME <- 'PhiPS2'

# Describe a few key features of the data
NUM_OBS_IN_SEQ <- 17
MEASUREMENT_NUMBERS_TO_REMOVE <- c(9, 10)

# Decide whether to make certain plots
MAKE_ANALYSIS_PLOTS <- TRUE

# Choose a maximum value of Ci to use when fitting (ppm). Set to Inf to disable.
MAX_CI <- Inf

# Decide which point to use for box plots of A and other quantities
POINT_FOR_BOX_PLOTS <- 1

# Specify gm value to use here
GM_VALUE <- Inf

# Specify a value of Gamma_star to use (this will override the value calculated
# from Arrhenius equations)
GAMMA_STAR <- 50 # ppm

###
### TRANSLATION:
### Creating convenient R objects from raw data files
###

# Define a vector of paths to the files we wish to load
file_paths <- c(
    '2022-07-09 36625 aci field ripe 10_1.xlsx',
    '2022-07-09 36625 aci field ripe 10_2.xlsx',
    '2022-07-09 36625 aci field ripe 11.xlsx',
    '2022-07-09 36625 aci field ripe15.xlsx',
    '2022-07-09 36625 field aci ripe 2.xlsx',
    '2022-07-10 36625 aci field ripe 10.xlsx',
    '2022-07-10 36625 aci field ripe 10_2.xlsx',
    '2022-07-10 36625 aci field ripe 10_3.xlsx',
    '2022-07-10 36625 aci field ripe 11.xlsx',
    '2022-07-10 36625 aci field ripe 15.xlsx',
    '2022-07-10 36625 aci field ripe 2.xlsx'
)

# Load each file, storing the result in a list
licor_exdf_list <- lapply(file_paths, function(fpath) {
  read_gasex_file(fpath, 'time')
})

# Get the names of all columns that are present in all of the Licor files
columns_to_keep <- do.call(identify_common_columns, licor_exdf_list)

# Extract just these columns
licor_exdf_list <- lapply(licor_exdf_list, function(x) {
  x[ , columns_to_keep, TRUE]
})

# Use `rbind` to combine all the data
licor_data <- do.call(rbind, licor_exdf_list)

###
### VALIDATION:
### Organizing the data, checking its consistency and quality, cleaning it
###

# Add a column that combines `plot` and `replicate`
licor_data[, paste0(PLOT_COLUMN_NAME, REP_COLUMN_NAME)] <-
    paste(licor_data[, PLOT_COLUMN_NAME], licor_data[, REP_COLUMN_NAME])

# Use this instead of the original replicate column
REP_COLUMN_NAME <- paste0(PLOT_COLUMN_NAME, REP_COLUMN_NAME)

# Add a column that combines `event` and `replicate` that we can use to identify
# each curve in the data set
licor_data[, 'curve_identifier'] <-
    paste(licor_data[, EVENT_COLUMN_NAME], licor_data[, REP_COLUMN_NAME])

# Factorize ID columns
licor_data <- factorize_id_column(licor_data, EVENT_COLUMN_NAME)
licor_data <- factorize_id_column(licor_data, 'curve_identifier')

# Remove certain events
#licor_data <- remove_points(licor_data, list(event = c('15', '37')))

# Make sure the data meets basic requirements
check_licor_data(licor_data, 'curve_identifier', NUM_OBS_IN_SEQ)

# Remove points with duplicated `CO2_r_sp` values and order by `Ci`
licor_data <- organize_response_curve_data(
    licor_data,
    'curve_identifier',
    MEASUREMENT_NUMBERS_TO_REMOVE,
    'Ci'
)

###
### PROCESSING:
### Extracting new pieces of information from the data
###

# Include gm values (required for apply_gm)
licor_data <- set_variable(
    licor_data,
    'gmc',
    "mol m^(-2) s^(-1) bar^(-1)",
    'c3_co2_response_v2',
    GM_VALUE
)

# Calculate total pressure (required for apply_gm)
licor_data <- calculate_total_pressure(licor_data)

# Calculate additional gas properties
licor_data <- calculate_gas_properties(licor_data)

# Calculate Cc
licor_data <- apply_gm(licor_data)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_data <- calculate_arrhenius(licor_data, c3_arrhenius_sharkey)

# Manually override Gamma_star
licor_data[, 'Gamma_star'] <- GAMMA_STAR

# Calculate intrinsic water-use efficiency
licor_data <- calculate_wue(licor_data)

# Truncate the Ci range for fitting
licor_data_for_fitting <- licor_data[licor_data[, 'Ci'] <= MAX_CI, , TRUE]

# Fit the C3 A-Ci curves
c3_aci_results <- consolidate(by(
  licor_data_for_fitting,                       # The `exdf` object containing the curves
  licor_data_for_fitting[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
  fit_c3_aci,                                   # The function to apply to each chunk of `licor_data`
  Ca_atmospheric = 420,                         # The atmospheric CO2 concentration
  cj_crossover_min = 100,                       # Wj must be > Wc when Cc < this value (ppm)
  cj_crossover_max = 800,                       # Wj must be < Wc when Cc > this value (ppm)
  fixed = c(NA, NA, NA, NA)
))

# Plot the C3 A-Ci fits (including limiting rates)
pdf_print(xyplot(
  A + Ac + Aj + Ap + A_fit ~ Cc | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto.key = list(space = 'right'),
  grid = TRUE,
  xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
  ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')'),
  par.settings = list(
    superpose.line = list(col = multi_curve_line_colors()),
    superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
  ),
  ylim = c(-5, 55),
  curve_ids = c3_aci_results$fits[, 'curve_identifier'],
  panel = function(...) {
    panel.xyplot(...)
    args <- list(...)
    curve_id <- args$curve_ids[args$subscripts][1]
    fit_param <-
      c3_aci_results$parameters[c3_aci_results$parameters[, 'curve_identifier'] == curve_id, ]
    panel.points(
        fit_param$operating_An_model ~ fit_param$operating_Cc,
        type = 'p',
        col = 'black',
        pch = 1
    )
  }
))

# Plot the C3 A-Ci fits
pdf_print(xyplot(
  A + A_fit ~ Ci | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste0('Intercellular CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
  ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')')
))

###
### SYNTHESIS:
### Using plots and statistics to help draw conclusions from the data
###

# Create a few data frames that will be helpful for plots

all_samples <- c3_aci_results$fits$main_data

col_to_average_as <- c(
  'A', 'iWUE', PHIPS2_COLUMN_NAME, 'ETR', 'Ci', 'Cc', 'gsw'
)

all_samples_one_point <- all_samples[all_samples$seq_num == POINT_FOR_BOX_PLOTS, ]

aci_parameters <- c3_aci_results$parameters$main_data

# Make box-whisker plots and bar charts

boxplot_caption <- paste0(
    "Quartiles for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where CO2 setpoint = ",
    all_samples_one_point[, 'CO2_r_sp'][POINT_FOR_BOX_PLOTS],
    ")"
)

fitting_caption <- "Values obtained by fitting A vs. Ci using the FvCB model"

x_s <- all_samples_one_point[, EVENT_COLUMN_NAME]
x_v <- aci_parameters[, EVENT_COLUMN_NAME]
xl <- "Genotype"

plot_param <- list(
  list(Y = all_samples_one_point[, 'A'],    X = x_s, xlab = xl, ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",          ylim = c(0, 35),  main = boxplot_caption),
  list(Y = all_samples_one_point[, 'iWUE'], X = x_s, xlab = xl, ylab = "Intrinsic water use efficiency (micromol CO2 / mol H2O)", ylim = c(0, 100), main = boxplot_caption),
  list(Y = aci_parameters[, 'Vcmax_at_25'], X = x_v, xlab = xl, ylab = "Apparent Vcmax at 25 degrees C (micromol / m^2 / s)",     ylim = c(0, 150), main = fitting_caption)
)

invisible(lapply(plot_param, function(x) {
  pdf_print(do.call(bwplot_wrapper, x))

  pdf_print(do.call(barchart_with_errorbars, x))
}))

# Make average response curve plots

rc_caption <- "Average response curves for each event"

x_ci <- all_samples[, 'Ci']
x_s <- all_samples[, 'seq_num']
x_e <- all_samples[, EVENT_COLUMN_NAME]

ci_lim <- c(-50, 1500)
a_lim <- c(-10, 55)
gsw_lim <- c(0, 0.7)

ci_lab <- "Intercellular [CO2] (ppm)"
a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
gsw_lab <- "Stomatal conductance to H2O (mol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"

avg_plot_param <- list(
    list(all_samples[, 'A'],   x_ci, x_s, x_e, xlab = ci_lab, ylab = a_lab,   xlim = ci_lim, ylim = a_lim),
    list(all_samples[, 'gsw'], x_ci, x_s, x_e, xlab = ci_lab, ylab = gsw_lab, xlim = ci_lim, ylim = gsw_lim)
)

invisible(lapply(avg_plot_param, function(x) {
    pdf_print(
        do.call(xyplot_avg_rc, c(x, list(
            type = 'b',
            pch = 20,
            auto = TRUE,
            grid = TRUE,
            main = rc_caption
        ))),
        width = 8,
        height = 6
    )
}))

# Save to CSV
tmp <- by(
  all_samples,
  all_samples$curve_identifier,
  function(x) {
    tmp2 <- data.frame(
      event = x[1, EVENT_COLUMN_NAME],
      replicate = x[1, 'replicate'],
      plot <- x[1, 'plot'],
      curve_identifier = x[1, 'curve_identifier']
    )

    for (cn in col_to_average_as) {
      tmp3 <- as.data.frame(t(data.frame(a = x[[cn]])))
      colnames(tmp3) <- paste0(cn, '_', x$seq_num)
      tmp2 <- cbind(tmp2, tmp3)
    }
    tmp2
  }
)

tmp <- do.call(rbind, tmp)

write.csv(tmp,                   'for_jmp.csv',               row.names = FALSE)
write.csv(all_samples,           "all_samples.csv",           row.names = FALSE)
write.csv(all_samples_one_point, "all_samples_one_point.csv", row.names = FALSE)
write.csv(aci_parameters,        "aci_parameters.csv",        row.names = FALSE)

# Print average operating point information
cat('\nAverage operating point Ci for each genotype:\n')
print(tapply(
    c3_aci_results$parameters[, 'operating_Ci'],
    c3_aci_results$parameters[, EVENT_COLUMN_NAME],
    mean
))
cat('\n')

cat('\nAverage operating point An (interpolated) for each genotype:\n')
print(tapply(
    c3_aci_results$parameters[, 'operating_An'],
    c3_aci_results$parameters[, EVENT_COLUMN_NAME],
    mean
))
cat('\n')

cat('\nAverage operating point An (modeled) for each genotype:\n')
print(tapply(
    c3_aci_results$parameters[, 'operating_An_model'],
    c3_aci_results$parameters[, EVENT_COLUMN_NAME],
    mean
))
cat('\n')
