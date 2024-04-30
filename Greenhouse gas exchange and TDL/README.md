# Greenhouse gas exchange and TDL

## Overview

This directory contains data and analysis scripts used to calculate isotope
discrimination and mesophyll conductance. Values from the output file
`Summary CGR3 greenhouse mesophyll conductance TDL.xlsx` were used to create
Figure 3 and Table S2 (in the Supplemental Information).

## Files

- `20230314 TGA Sampling System_SiteAvg.dat`,
  `20230315 TGA Sampling System_SiteAvg.dat`, and
  `TGA Sampling System_SiteAvg.dat`: Tunable diode laser (TDL) log files that
  were recorded at the same time as the gas exchange files described below.

- Thirty-six Microsoft Excel files whose names begin with `2023`: Licor LI-6800
  log files that were recorded at the same time as the TDL log files.

- `gm_from_tdl.R`: An R script that calculates mesophyll conductance and other
  key parameters from the measurements in the TDL and Licor log files.

- `Summary CGR3 greenhouse mesophyll conductance TDL.xlsx`: A table of mesophyll
  conductance and other values that was originally produced by running the
  `gm_from_tdl.R` script.

## Running the script

### Required R packages

The scripts in this directory require two R packages: `lattice` and `PhotoGEA`
(version `0.10.0`):

- The `lattice` package can be installed from within R using
  `install.packages('lattice')`.

- The required version of `PhotoGEA` can be installed from within R by calling
  `remotes::install_github('eloch216/PhotoGEA', ref = 'v0.10.0')`. Note that
  this command requires the `remotes` package, which can be installed using
  `install.packages('remotes')`.

### Running from R

1. Set the working directory of an R session to this directory.

2. Type the following to run the script:

   ```r
   source('gm_from_tdl.R')
   ```

3. The script will produce several graphs and output CSV files. In particular,
   the contents of `Summary CGR3 greenhouse mesophyll conductance TDL.xlsx` are
   taken from a subset of the columns included in the
   `gm_stats_by_rep_outliers_excluded.csv` file.

*WARNING:* This script will automatically clear your workspace before it runs.
