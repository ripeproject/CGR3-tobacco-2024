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
        'The `process_variable_j.R` script requires PhotoGEA version 0.10.0. ',
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

MEASUREMENT_NUMBERS_TO_REMOVE <- c(1, 6, 7, 8, 9, 15, 16, 17)

# Choose a maximum value of Ci to use when fitting (ppm). Set to Inf to disable.
MAX_CI <- Inf

# Decide which point to use for box plots of A and other quantities
POINT_FOR_BOX_PLOTS <- 10

# Specify which parameters should be fixed when fitting
FIXED <- c(NA, NA, NA, NA, NA) # J, Rd, tau, TPU, Vcmax

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

# Check data
check_licor_data(licor_data, 'curve_identifier')

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

# Calculate total pressure (required for fit_c3_variable_j)
licor_data <- calculate_total_pressure(licor_data)

# Calculate additional gas properties (required for calculate_c3_limitations_grassi)
licor_data <- calculate_gas_properties(licor_data)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_data <- calculate_arrhenius(licor_data, c3_arrhenius_sharkey)

# Calculate intrinsic water-use efficiency
licor_data <- calculate_wue(licor_data)

# Truncate the Ci range for fitting
licor_data_for_fitting <- licor_data[licor_data[, 'Ci'] <= MAX_CI, , TRUE]

# Set a seed number before fitting to make sure results are reproducible
set.seed(1234)

# Fit the C3 A-Ci curves using the variable J method
c3_aci_results <- consolidate(by(
  licor_data_for_fitting,                       # The `exdf` object containing the curves
  licor_data_for_fitting[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
  fit_c3_variable_j,                            # The function to apply to each chunk of `licor_data`
  Ca_atmospheric = 420,                         # The atmospheric CO2 concentration
  OPTIM_FUN = optimizer_nmkb(),                 # The optimization algorithm to use
  fixed = FIXED,
  cj_crossover_min = 20,                        # Wj must be > Wc when Cc < this value (ppm)
  cj_crossover_max = 800                        # Wj must be < Wc when Cc > this value (ppm)
))

# Plot the C3 A-Cc fits (including limiting rates)
pdf_print(xyplot(
  A + Ac + Aj + Ap + A_fit ~ Cc | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto.key = list(space = 'right'),
  grid = TRUE,
  xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Cc, ')'),
  ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')'),
  par.settings = list(
    superpose.line = list(col = multi_curve_line_colors()),
    superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
  ),
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

# Plot the C3 A-Ci fits (including limiting rates)
pdf_print(xyplot(
  A + Ac + Aj + Ap + A_fit ~ Ci | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto.key = list(space = 'right'),
  grid = TRUE,
  xlab = paste0('Intercellular CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
  ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')'),
  par.settings = list(
    superpose.line = list(col = multi_curve_line_colors()),
    superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
  ),
  curve_ids = c3_aci_results$fits[, 'curve_identifier'],
  panel = function(...) {
    panel.xyplot(...)
    args <- list(...)
    curve_id <- args$curve_ids[args$subscripts][1]
    fit_param <-
      c3_aci_results$parameters[c3_aci_results$parameters[, 'curve_identifier'] == curve_id, ]
    panel.points(
      fit_param$operating_An_model ~ fit_param$operating_Ci,
      type = 'p',
      col = 'black',
      pch = 1
    )
  }
))

# Plot the C3 A-Cc fits
pdf_print(xyplot(
  A + A_fit ~ Cc | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Cc, ')'),
  ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')')
))

# Plot the C3 J_F-Cc fits
pdf_print(xyplot(
  ETR + J_F ~ Cc | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Cc, ')'),
  ylab = paste0('Electron transport rate (', c3_aci_results$fits$units$J_F, ')')
))

# Plot the C3 gmc-Cc fits
pdf_print(xyplot(
  gmc ~ Cc | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Cc, ')'),
  ylab = paste0('Mesophyll conductance (', c3_aci_results$fits$units$gmc, ')')
))

# Plot the C3 gmc-Ci fits
pdf_print(xyplot(
  gmc ~ Ci | curve_identifier,
  data = c3_aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste0('Intercellular CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
  ylab = paste0('Mesophyll conductance (', c3_aci_results$fits$units$gmc, ')')
))

###
### SYNTHESIS:
### Using plots and statistics to help draw conclusions from the data
###

print(paste('Number of rows before removing outliers:', nrow(c3_aci_results$parameters)))

c3_aci_results$parameters <- exclude_outliers(
    c3_aci_results$parameters,
    'Vcmax_at_25',
    c3_aci_results$parameters[, EVENT_COLUMN_NAME]
)

print(paste('Number of rows after removing outliers:', nrow(c3_aci_results$parameters)))

c3_aci_results$fits <-
    c3_aci_results$fits[c3_aci_results$fits[, 'curve_identifier'] %in% c3_aci_results$parameters[, 'curve_identifier'], , TRUE]

# Create a few data frames that will be helpful for plots
all_samples <- c3_aci_results$fits$main_data

col_to_average_as <- c(
  'A', 'iWUE', PHIPS2_COLUMN_NAME, 'ETR', 'Ci', 'Cc', 'gsw', 'gmc'
)

all_samples_one_point <- all_samples[all_samples$seq_num == POINT_FOR_BOX_PLOTS, ]
aci_parameters <- c3_aci_results$parameters$main_data

# Make box-whisker plots and bar charts

boxplot_caption <- paste0(
    'Quartiles for measurement point ',
    POINT_FOR_BOX_PLOTS,
    '\n(where CO2 setpoint = ',
    all_samples_one_point[, 'CO2_r_sp'][POINT_FOR_BOX_PLOTS],
    ')'
)

fitting_caption <- 'Values obtained by fitting A vs. Ci using the Variable J method'

x_s <- all_samples_one_point[, EVENT_COLUMN_NAME]
x_v <- aci_parameters[, EVENT_COLUMN_NAME]
xl <- 'Genotype'

plot_param <- list(
  list(Y = all_samples_one_point[, 'gmc'],  X = x_s, xlab = xl, ylab = 'Mesophyll conductance (mol / m^2 / s / bar)', ylim = c(0, 0.25), main = boxplot_caption),
  list(Y = aci_parameters[, 'Vcmax_at_25'], X = x_v, xlab = xl, ylab = 'Vcmax at 25 degrees C (micromol / m^2 / s)',  ylim = c(0, 450),  main = fitting_caption),
  list(Y = aci_parameters[, 'tau'],         X = x_v, xlab = xl, ylab = 'tau (dimensionless)',                         ylim = c(0, 1),    main = fitting_caption)
)

invisible(lapply(plot_param, function(x) {
  pdf_print(do.call(bwplot_wrapper, x))

  pdf_print(do.call(barchart_with_errorbars, x))
}))

# Make average response curve plots

rc_caption <- 'Average response curves for each event'

x_ci <- all_samples[, 'Ci']
x_cc <- all_samples[, 'Cc']
x_s  <- all_samples[, 'seq_num']
x_e  <- all_samples[, EVENT_COLUMN_NAME]

ci_lim  <- c(0, 800)
cc_lim  <- c(0, 200)
a_lim   <- c(-10, 45)
gmc_lim <- c(0, 0.25)

ci_lab  <- 'Intercellular [CO2] (ppm)'
cc_lab  <- 'Chloroplast [CO2] (ppm)'
a_lab   <- 'Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)'
gmc_lab <- 'Mesophyll conductance (mol / m^2 / s / bar)\n(error bars: standard error of the mean for same CO2 setpoint)'

avg_plot_param <- list(
    list(all_samples[, 'A'],   x_ci, x_s, x_e, xlab = ci_lab, ylab = a_lab,   xlim = ci_lim, ylim = a_lim),
    list(all_samples[, 'A'],   x_cc, x_s, x_e, xlab = cc_lab, ylab = a_lab,   xlim = cc_lim, ylim = a_lim),
    list(all_samples[, 'gmc'], x_cc, x_s, x_e, xlab = cc_lab, ylab = gmc_lab, xlim = cc_lim, ylim = gmc_lim),
    list(all_samples[, 'gmc'], x_ci, x_s, x_e, xlab = ci_lab, ylab = gmc_lab, xlim = ci_lim, ylim = gmc_lim)
)

invisible(lapply(avg_plot_param, function(x) {
    dev.new(width = 8, height = 6)
    print(do.call(xyplot_avg_rc, c(x, list(
        type = 'b',
        pch = 20,
        auto = TRUE,
        grid = TRUE,
        main = rc_caption
    ))))
}))

# Save to CSV
tmp <- by(
  all_samples,
  all_samples$curve_identifier,
  function(x) {
    tmp2 <- data.frame(
      event = x[1, EVENT_COLUMN_NAME],
      plot = x[1, 'plot'],
      replicate = x[1, 'replicate'],
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

write.csv(tmp,                   'vj_for_jmp.csv',               row.names = FALSE)
write.csv(all_samples,           'vj_all_samples.csv',           row.names = FALSE)
write.csv(all_samples_one_point, 'vj_all_samples_one_point.csv', row.names = FALSE)
write.csv(aci_parameters,        'vj_aci_parameters.csv',        row.names = FALSE)
