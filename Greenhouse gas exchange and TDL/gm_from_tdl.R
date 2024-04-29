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

# Define Rd
RESPIRATION_RATE <- 2.1

# Define Rubisco specificity
RUBISCO_SPECIFICITY_AT_TLEAF <- 97.3 # Bernacchi et al. (2002), tobacco, in vivo,  25 C

# Define gm and Cc cutoffs
MIN_GM <- 0
MAX_GM <- 3
MIN_CC <- 0.0

# Names of important columns in the TDL data
TDL_TIMESTAMP_COLUMN_NAME <- 'TIMESTAMP'
TDL_VALVE_COLUMN_NAME <- 'valve_number'
TDL_12C_COLUMN_NAME <- 'Conc12C_Avg'
TDL_13C_COLUMN_NAME <- 'Conc13C_Avg'

# Specify the variables to extract from the TDL data files. Note that when the
# files are loaded, any Unicode characters such as Greek letters will be
# converted into `ASCII` versions, e.g. the character Δ will be become `Delta`.
# The conversion rules are defined in the `UNICODE_REPLACEMENTS` data frame
# (see `read_licor.R`).
TDL_COLUMNS_TO_EXTRACT <- c(
    TDL_TIMESTAMP_COLUMN_NAME,
    TDL_VALVE_COLUMN_NAME,
    TDL_12C_COLUMN_NAME,
    TDL_13C_COLUMN_NAME
)

# Names of important columns in the Licor data
LICOR_TIMESTAMP_COLUMN_NAME <- 'time'

# Choose gas lines to smooth
TDL_VALVES_TO_SMOOTH <- c(2, 20, 21, 23, 26)

# Use splines to smooth the data
spline_smoothing_function <- function(Y, X) {
    ss <- smooth.spline(X, Y)
    return(ss$y)
}

# Define a helping function for getting genotype information from filenames
get_genotype_info_from_licor_filename <- function(licor_exdf) {
    # Add some new columns to the Licor file in preparation for adding the plant
    # information
    licor_exdf <- document_variables(
        licor_exdf,
        c("plant specification", "genotype",                 "NA"),
        c("plant specification", "event",                    "NA"),
        c("plant specification", "replicate",                "NA"),
        c("plant specification", "genotype_event",           "NA"),
        c("plant specification", "event_replicate",          "NA"),
        c("plant specification", "genotype_event_replicate", "NA"),
        c("plant specification", "original_file",            "NA")
    )

    # Get the filename without the path
    name <- basename(licor_exdf[['file_name']])

    # There are two naming conventions we need to check

    # Search for the following pattern: one space, followed by one or
    # more alphanumeric characters, followed by a dash, followed by one or more
    # alphanumeric characters, followed by a dash, followed by one or more
    # alphanumeric characters, followed by a period. Essentially, we expect the
    # filename to end with ' XXX-YYY-ZZZ.xlsx', where `XXX` is the genotype,
    # `YYY` is the event, and `ZZZ` is the replicate.
    plant_specification1 <-
        regmatches(name, regexpr(" [[:alnum:]]+-[[:alnum:]]+-[[:alnum:]]+\\.xlsx", name))

    # Search for the following pattern: one space, followed by one or
    # more alphanumeric characters, followed by a space, followed by one or more
    # alphanumeric characters, followed by a space, followed by one or more
    # alphanumeric characters, followed by a period. Essentially, we expect the
    # filename to end with ' XXX YYY ZZZ.xlsx', where `XXX` is the genotype,
    # `YYY` is the event, and `ZZZ` is the replicate.
    plant_specification2 <-
        regmatches(name, regexpr(" [[:alnum:]]+ [[:alnum:]]+ [[:alnum:]]+\\.xlsx", name))

    # Make sure we found something
    if (length(plant_specification1) == 0 & length(plant_specification2) == 0) {
        msg <- paste0(
            "Could not extract plant specification information from Licor file:\n'",
            licor_exdf[['file_name']],
            "'\nThe filename must end with ` GGG-EEE-RRR.xlsx` or ` GGG EEE RRR.xlsx`, ",
            "where `GGG`, `EEE`, and `RRR` are alphanumeric specifiers for ",
            "the genotype, event, and replicate represented by the file"

        )
        stop(msg)
    }

    # Extract the info
    g <- character(0)
    e <- character(0)
    r <- character(0)
    if (length(plant_specification1) > 0) {
        # Remove the period, the extension, and the whitespace
        plant_specification1 <- sub(" ", "", plant_specification1)
        plant_specification1 <- sub("\\.xlsx", "", plant_specification1)

        # Split the specification by the dashes
        plant_specification1 <- strsplit(plant_specification1, "-")[[1]]

        g <- plant_specification1[1]
        e <- plant_specification1[2]
        r <- plant_specification1[3]
    } else {
        # Remove the period and the extension
        plant_specification2 <- sub("\\.xlsx", "", plant_specification2)

        # Split the specification by the spaces
        plant_specification2 <- strsplit(plant_specification2, " ")[[1]]

        g <- plant_specification2[2]
        e <- plant_specification2[3]
        r <- plant_specification2[4]
    }

    # Store the info in the file and return it
    licor_exdf[,'genotype'] <- g
    licor_exdf[,'event'] <- e
    licor_exdf[,'replicate'] <- r

    licor_exdf[, 'event_replicate'] <-
        paste(licor_exdf[, 'event'], licor_exdf[, 'replicate'])

    licor_exdf[,'genotype_event'] <-
        paste(licor_exdf[, 'genotype'], licor_exdf[, 'event'])

    licor_exdf[,'event_replicate'] <-
        paste(licor_exdf[, 'event'], licor_exdf[, 'replicate'])

    licor_exdf[,'genotype_event_replicate'] <-
        paste(licor_exdf[, 'genotype'], licor_exdf[, 'event'], licor_exdf[, 'replicate'])

    licor_exdf[,'original_file'] <- licor_exdf[['file_name']]

    return(licor_exdf)
}

###
### TRANSLATION:
### Creating convenient R objects from raw data files
###

# Get all the TDL information and process it

tdl_file_paths <- c(
    "20230314 TGA Sampling System_SiteAvg.dat",
    "20230315 TGA Sampling System_SiteAvg.dat",
    "TGA Sampling System_SiteAvg.dat"
)

tdl_files <- lapply(tdl_file_paths, function(fname) {
    read_gasex_file(
        fname,
        rows_to_skip = 1,
        variable_name_row = 2,
        variable_unit_row = 3,
        data_start_row = 5,
        timestamp_colname = TDL_TIMESTAMP_COLUMN_NAME
    )
})

tdl_columns_to_extract <- c(
    TDL_TIMESTAMP_COLUMN_NAME,
    TDL_VALVE_COLUMN_NAME,
    'Conc12C_Avg',
    'Conc13C_Avg'
)

tdl_files <- lapply(tdl_files, function(exdf_obj) {
    exdf_obj[ , tdl_columns_to_extract, TRUE]
})

tdl_files <- do.call(rbind, tdl_files)

tdl_files <- identify_tdl_cycles(
    tdl_files,
    valve_column_name = TDL_VALVE_COLUMN_NAME,
    cycle_start_valve = 20,
    expected_cycle_length_minutes = 3,
    expected_cycle_num_valves = 9,
    timestamp_colname = TDL_TIMESTAMP_COLUMN_NAME
)

tdl_files_smoothed <- tdl_files

for (valve in TDL_VALVES_TO_SMOOTH) {
    for (column in c(TDL_12C_COLUMN_NAME, TDL_13C_COLUMN_NAME)) {
        tdl_files_smoothed <- smooth_tdl_data(
            tdl_files_smoothed,
            column,
            TDL_VALVE_COLUMN_NAME,
            valve,
            spline_smoothing_function
        )
    }
}

processed_tdl_data <- consolidate(by(
    tdl_files_smoothed,
    tdl_files_smoothed[, 'cycle_num'],
    process_tdl_cycle_erml,
    valve_column_name = TDL_VALVE_COLUMN_NAME,
    noaa_valve = 2,
    calibration_0_valve = 20,
    calibration_1_valve = 21,
    calibration_2_valve = 23,
    calibration_3_valve = 26,
    raw_12c_colname = TDL_12C_COLUMN_NAME,
    raw_13c_colname = TDL_13C_COLUMN_NAME,
    noaa_cylinder_co2_concentration = 294.996, # ppm
    noaa_cylinder_isotope_ratio = -8.40,       # ppt
    calibration_isotope_ratio = -11.505        # ppt
))

# Get all the Licor information and process it

licor_file_paths <- c(
    "2023-03-14 long 4 site 11 36625 WT 4.xlsx",
    "2023-03-14 long 4 site 11 36625-10-6.xlsx",
    "2023-03-14 long 4 site 11 36625-14-3.xlsx",
    "2023-03-14 long 4 site 11 36625-14-4.xlsx",
    "2023-03-14 long 4 site 11 36625-8-2.xlsx",
    "2023-03-14 long 4 site 11 36625-8-5.xlsx",
    "2023-03-14 pluto site 13 36625 WT 2.xlsx",
    "2023-03-14 pluto site 13 36625 WT 7.xlsx",
    "2023-03-14 pluto site 13 36625-10-7.xlsx",
    "2023-03-14 pluto site 13 36625-14-6.xlsx",
    "2023-03-14 pluto site 13 36625-8-3.xlsx",
    "2023-03-15 long 4 site 11 36625 WT 10.xlsx",
    "2023-03-15 long 4 site 11 36625-10-3.xlsx",
    "2023-03-15 long 4 site 11 36625-10-9.xlsx",
    "2023-03-15 long 4 site 11 36625-14-8.xlsx",
    "2023-03-15 long 4 site 11 36625-8-8.xlsx",
    "2023-03-15 pluto site 13 36625 WT 3.xlsx",
    "2023-03-15 pluto site 13 36625-10-4.xlsx",
    "2023-03-15 pluto site 13 36625-14-2.xlsx",
    "2023-03-15 pluto site 13 36625-14-9.xlsx",
    "2023-03-15 pluto site 13 36625-8-4.xlsx",
    "2023-03-15 pluto site 13 36625-8-6.xlsx",
    "2023-03-16 long 4 site 11 36625 WT 1.xlsx",
    "2023-03-16 long 4 site 11 36625 WT 5.xlsx",
    "2023-03-16 long 4 site 11 36625-10-1.xlsx",
    "2023-03-16 long 4 site 11 36625-10-8.xlsx",
    "2023-03-16 long 4 site 11 36625-14-1.xlsx",
    "2023-03-16 long 4 site 11 36625-8-9.xlsx",
    "2023-03-16 pluto site 13 36625 WT 6.xlsx",
    "2023-03-16 pluto site 13 36625 WT 9.xlsx",
    "2023-03-16 pluto site 13 36625-10-10.xlsx",
    "2023-03-16 pluto site 13 36625-14-10.xlsx",
    "2023-03-16 pluto site 13 36625-14-5.xlsx",
    "2023-03-16 pluto site 13 36625-14-7.xlsx",
    "2023-03-16 pluto site 13 36625-8-10.xlsx",
    "2023-03-16 pluto site 13 36625-8-7.xlsx"
)

licor_files <- lapply(licor_file_paths, function(fname) {
    read_gasex_file(fname, LICOR_TIMESTAMP_COLUMN_NAME)
})

common_columns <- do.call(identify_common_columns, licor_files)

licor_files <- lapply(licor_files, function(exdf_obj) {
    exdf_obj[ , common_columns, TRUE]
})

licor_files <- lapply(licor_files, get_genotype_info_from_licor_filename)

licor_files <- lapply(licor_files, get_oxygen_from_preamble)

licor_files <- lapply(licor_files, function(x) {
    get_sample_valve_from_filename(x, list(
        '13' = 12,
        '11' = 10
    ))
})

# Combine the Licor and TDL data
licor_files <- lapply(licor_files, function(licor_exdf) {
    pair_gasex_and_tdl(
        licor_exdf,
        processed_tdl_data$tdl_data
    )
})

licor_files <- do.call(rbind, licor_files)

# Factorize ID columns
licor_files <- factorize_id_column(licor_files, 'event')
licor_files <- factorize_id_column(licor_files, 'event_replicate')

###
### PROCESSING:
### Extracting new pieces of information from the data
###

# Specify respiration values
licor_files <- set_variable(
    licor_files,
    'Rd',
    'micromol m^(-2) s^(-1)',
    'gm_from_tdl',
    abs(RESPIRATION_RATE)
)

# Specify Rubisco specificity values
licor_files <- set_variable(
    licor_files,
    'specificity_at_tleaf',
    'M / M',
    'gm_from_tdl',
    RUBISCO_SPECIFICITY_AT_TLEAF
)

# Calculate total pressure (needed for `calculate_gas_properties`)
licor_files <- calculate_total_pressure(licor_files)

# Calculate gbc, gsc, Csurface (needed for `calculate_gm_ubierna`)
licor_files <- calculate_gas_properties(licor_files)

# Calculate Delta_obs_tdl (needed for `calculate_gm_busch`)
licor_files <- calculate_isotope_discrimination(licor_files)

# Calculates t (needed for `calculate_gm_busch`)
licor_files <- calculate_ternary_correction(licor_files)

# Calculate Gamma_star (needed for `calculate_gm_busch`)
licor_files <- calculate_gamma_star(licor_files)

# Here we use Equation 19 for e_star because we don't have values for
# delta_obs_growth
licor_files <- calculate_gm_busch(
    licor_files,
    e_star_equation = 19
)

licor_files <- apply_gm(licor_files)

licor_files <- calculate_wue(
    licor_files,
    calculate_c3 = TRUE
)

###
### VALIDATION:
### Organizing the data, checking its consistency and quality, cleaning it
###

# Check the data for any issues before proceeding with additional analysis
check_licor_data(licor_files, 'event_replicate', -1)

# Exclude any points where gm or Cc is outside the acceptable range
licor_files_no_outliers <- licor_files
cat(paste("Total number of Licor measurements:", nrow(licor_files_no_outliers), "\n"))

licor_files_no_outliers[['main_data']] <-
    licor_files_no_outliers[['main_data']][
        licor_files_no_outliers[['main_data']][['gmc']] > MIN_GM &
        licor_files_no_outliers[['main_data']][['gmc']] < MAX_GM &
        licor_files_no_outliers[['main_data']][['Cc']] > MIN_CC,]

cat(paste("Number of Licor measurements after removing unacceptable gm and Cc:", nrow(licor_files_no_outliers), "\n"))

# Now we remove outliers from each replicate, based on A and gmc
licor_files_no_outliers <- exclude_outliers(
    licor_files_no_outliers,
    'A',
    licor_files_no_outliers[,'event_replicate']
)

licor_files_no_outliers <- exclude_outliers(
    licor_files_no_outliers,
    'gmc',
    licor_files_no_outliers[,'event_replicate']
)

cat(paste("Number of Licor measurements after removing statistical outliers from each rep:", nrow(licor_files_no_outliers), "\n"))

# Get stats for each rep by averaging over all corresponding observations
rep_stats_no_outliers <- basic_stats(
    licor_files_no_outliers,
    'event_replicate'
)

# Get stats for each event by averaging over all corresponding reps
event_stats <- basic_stats(
  rep_stats_no_outliers,
  'event'
)$main_data

# Save results as CSV files
write.csv(processed_tdl_data[['tdl_data']], "tdl_data_processed.csv",                  row.names=FALSE)
write.csv(licor_files,                      "gm_calculations_outliers_included.csv",   row.names=FALSE)
write.csv(licor_files_no_outliers,          "gm_calculations_outliers_excluded.csv",   row.names=FALSE)
write.csv(rep_stats_no_outliers$main_data,  "gm_stats_by_rep_outliers_excluded.csv",   row.names=FALSE)
write.csv(event_stats,                      "gm_stats_by_event_outliers_excluded.csv", row.names=FALSE)

###
### SYNTHESIS:
### Using plots and statistics to help draw conclusions from the data
###

# Remove one remaining outlier before plotting
rep_stats_no_outliers <- remove_points(
    rep_stats_no_outliers,
    list(event = 14, replicate = 5) # Has unreasonably high gmc
)

# Extract data frame from rep stats
rep_stats_no_outliers <- rep_stats_no_outliers$main_data

# Define plotting parameters
x_e <- rep_stats_no_outliers[['event']]
x_g <- rep_stats_no_outliers[['genotype']]
x_er <- rep_stats_no_outliers[['event_replicate']]

gmc_lab      <- "Mesophyll conductance to CO2 (mol / m^2 / s / bar)"
cc_lab       <- "CO2 concentration in chloroplast (micromol / mol)"
drawdown_lab <- "CO2 drawdown across mesophyll (Ci - Cc) (micromol / mol)"
a_lab        <- "Net CO2 assimilation rate (micromol / m^2 / s)"
iwue_lab     <- "Intrinsic water use efficiency (micromol CO2 / mol H2O)"
g_ratio_lab  <- "Ratio of mesophyll / stomatal conductances to CO2 (gm / gs; dimensionless)"
dtdl_lab     <- "Delta13c (ppt)"

gmc_lim      <- c(0, 1)
cc_lim       <- c(0, 275)
drawdown_lim <- c(0, 100)
a_lim        <- c(0, 50)
iwue_lim     <- c(0, 120)
g_ratio_lim  <- c(0, 1.5)
dtdl_lim     <- c(0, 25)

box_plot_param <- list(
  list(Y = rep_stats_no_outliers[['gmc_avg']],           X = x_er, S = x_g, ylab = gmc_lab,      ylim = gmc_lim),
  list(Y = rep_stats_no_outliers[['Cc_avg']],            X = x_er, S = x_g, ylab = cc_lab,       ylim = cc_lim),
  list(Y = rep_stats_no_outliers[['drawdown_cm_avg']],   X = x_er, S = x_g, ylab = drawdown_lab, ylim = drawdown_lim),
  list(Y = rep_stats_no_outliers[['A_avg']],             X = x_er, S = x_g, ylab = a_lab,        ylim = a_lim),
  list(Y = rep_stats_no_outliers[['iWUE_avg']],          X = x_er, S = x_g, ylab = iwue_lab,     ylim = iwue_lim),
  list(Y = rep_stats_no_outliers[['g_ratio_avg']],       X = x_er, S = x_g, ylab = g_ratio_lab,  ylim = g_ratio_lim),
  list(Y = rep_stats_no_outliers[['Delta_obs_tdl_avg']], X = x_er, S = x_g, ylab = dtdl_lab,     ylim = dtdl_lim)
)

box_bar_plot_param <- list(
  list(Y = rep_stats_no_outliers[['gmc_avg']],           X = x_e,  S = x_g, ylab = gmc_lab,      ylim = gmc_lim),
  list(Y = rep_stats_no_outliers[['Cc_avg']],            X = x_e,  S = x_g, ylab = cc_lab,       ylim = cc_lim),
  list(Y = rep_stats_no_outliers[['drawdown_cm_avg']],   X = x_e,  S = x_g, ylab = drawdown_lab, ylim = drawdown_lim),
  list(Y = rep_stats_no_outliers[['A_avg']],             X = x_e,  S = x_g, ylab = a_lab,        ylim = a_lim),
  list(Y = rep_stats_no_outliers[['iWUE_avg']],          X = x_e,  S = x_g, ylab = iwue_lab,     ylim = iwue_lim),
  list(Y = rep_stats_no_outliers[['g_ratio_avg']],       X = x_e,  S = x_g, ylab = g_ratio_lab,  ylim = g_ratio_lim),
  list(Y = rep_stats_no_outliers[['Delta_obs_tdl_avg']], X = x_e,  S = x_g, ylab = dtdl_lab,     ylim = dtdl_lim)
)

# Make all the box and bar charts
invisible(lapply(box_plot_param, function(x) {
  dev.new()
  print(do.call(bwplot_wrapper, x))
}))

invisible(lapply(box_bar_plot_param, function(x) {
  dev.new()
  print(do.call(bwplot_wrapper, x))

  dev.new()
  print(do.call(barchart_with_errorbars, x))
}))
