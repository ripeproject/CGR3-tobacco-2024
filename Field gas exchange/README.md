# Field gas exchange

## Overview

This directory contains data and analysis scripts used to fit CO2 response
curves measured along with chlorophyll fluorescence.

Values from the output file `CGR3 field gas exchange ACi fits.xlsx` were used to
create Figures 4c and 5c.

Values from the output file `CGR3 field gas exchange gm variable J.xlsx` were
used to create Figures 4a, 4b, 4d, 5a, 5b, and 5d, and Table S3 (in the
Supplemental Information).

## Files

- Eleven Microsoft Excel files whose names begin with `2022`: Licor LI-6800
  log files.

- `process_aci.R`: An R script that estimates effective Vcmax and other
  photosynthetic parameters from the measured gas exchange data (without using
  the chlorophyll fluorescence data, and by setting `Cc = Ci`)

- `process_variable_j.R`: An R script that estimates mesophyll conductance,
  `Cc`, and other photosynthetic parameters from the combined gas exchange and
  chlorophyll fluorescence data

- `CGR3 field gas exchange ACi fits.xlsx`: An Excel file that contains a table
  whose values were originally produced by running the `process_aci.R` script.

- `CGR3 field gas exchange gm variable J.xlsx`: An Excel file that contains a
  table whose values were originally produced by running the
  `process_variable_j.R` script

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

2. Type the following to run the A-Ci script:

   ```r
   source('process_aci.R')
   ```

3. The `process_aci` script will produce several graphs and output CSV
   files. In particular, the contents of
   `CGR3 field gas exchange ACi fits.xlsx` are taken from `aci_parameters.csv`
   and `for_jmp.csv`.

4. Type the following to run the Variable J script:

   ```r
   source('process_variable_j.R')
   ```

5. The `process_variable_j` script will produce several graphs and output CSV
   files. In particular, the contents of
   `CGR3 field gas exchange gm variable J.xlsx` are taken from
   `vj_aci_parameters.csv` and `vj_for_jmp.csv`.

*WARNING:* This script will automatically clear your workspace before it runs.

## Data dictionary

### Outputs

The output files `CGR3 field gas exchange ACi fits.xlsx` and
`CGR3 field gas exchange gm variable J.xlsx` contain columns with names
formatted like `Parameter_X`, where parameter is `A`, `iWUE`, etc and `X`
corresponds to the number along the measurement sequence. The CO2 response
curves were measured using the following sequence of reference CO2 (with units
of `ppm`):

1. 400
2. 300
3. 200
4. 150
5. 100
6. 75
7. 50
8. 20
9. 400
10. 400
11. 500
12. 600
13. 800
14. 1000
15. 1200
16. 1500
17. 1800

So, for example, `A_3` would refer to the measured net CO2 assimilation rate
when the reference CO2 concentration was set to 200 ppm (the third point along
the sequence).
