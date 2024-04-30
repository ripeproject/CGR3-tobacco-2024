# Field gas exchange

## Overview

This directory contains data and analysis scripts used to fit CO2 response
curves measured along with chlorophyll fluorescence. Values from the output
file `CGR3 field gas exchange (ACi and gm variable J).xlsx` were used to create
Figures 4a, 4b, 5a, 5b, and 5d, and Table Table S3 (in the Supplemental
Information).

## Files

- Eleven Microsoft Excel files whose names begin with `2022`: Licor LI-6800
  log files.

- `process_variable_j.R`: An R script that estimates mesophyll conductance and
  photosynthetic parameters from the combined gas exchange and chlorophyll
  fluorescence data

- `CGR3 field gas exchange (ACi and gm variable J).xlsx`: An Excel file with two
  sheets. The sheet called `gm Variable J` contains a table whose values were
  originally produced by running the `process_variable_j.R` script.

## Running the scripts

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

2. Type the following to run the Variable J script:

   ```r
   source('process_variable_j.R')
   ```

3. The `process_variable_j` script will produce several graphs and output CSV
   files. In particular, the contents of the `gm Variable J` sheet in
   `CGR3 field gas exchange (ACi and gm variable J).xlsx` are taken from
   `vj_aci_parameters.csv` and `vj_for_jmp.csv`.

*WARNING:* This script will automatically clear your workspace before it runs.
