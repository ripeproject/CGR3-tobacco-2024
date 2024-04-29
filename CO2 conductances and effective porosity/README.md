# CO2 conductances and effective porosity

## Overview

This directory contains data and analysis scripts used to calculate CO2
conductances and effective porosity. Values from the output file
`CGR3 greenhouse CO2 conductances and effective porosity.csv` were used to
create Figures 1e and 1f. For more details about these calculations, see
Methods S2 in the Supplemenal Information.

## Files

- `leaf_properties.csv`: A table of key input values needed to calculate CO2
  conductances and effective porosity.

- `effective_porosity.R`: An R script that calculates CO2 conductances and
  effective porosity from the contents of `leaf_properties.csv`.

- `CGR3 greenhouse CO2 conductances and effective porosity.csv`: A table of CO2
  conductances and effective porosities that was originally produced by running
  the `effective_porosity.R` script.

## Running the script

1. Set the working directory of an R session to this directory.

2. Type the following to run the script:

   ```r
   source('effective_porosity.R')
   ```

3. The script will produce an output file called `effective_porosity.csv` that
   will be created in this directory. It should be identical to
   `CGR3 greenhouse CO2 conductances and effective porosity.csv`, other than
   small changes in formatting and significant digits introduced by Microsoft
   Excel.

Note: The `effective_porosity.R` script only uses commands from base R, so no
packages are required to use it.

## Data dictionary

### Inputs

The columns in `leaf_properties.csv` are defined as follows:

- `Tmes`: Mesophyll thickness. Units: `μm`

- `Tleaf`: Leaf temperature during gas exchange measurements. Units: `degrees C`

- `Tcw`: Cell wall thickness. Units: `nm`

- `Sc`: Surface area of chloroplasts exposed to the mesophyll, normalized by the
  leaf surface area. Units: `dimensionless` from `m^2 chloroplast / m^2 leaf`

- `gmc`: Mesophyll conductance to CO2 diffusion. Units: `mol / m^2 / s/ bar`

- `Fias`: Fraction of intercellular airspace. Units: `dimensionless`

### Outputs

The columns in `CGR3 greenhouse CO2 conductances and effective porosity.csv` are
defined as follows:

 - `conv`: Multiplicative factor used to convert some conductances from `m / s`
   to `mol / m^2 / s/ bar`. Units: `m^3 bar / mol`

 - `eps_membrane`: Membrane conductance enhancement factor.
   Units: `dimensionless`

 - `Dair`: Diffusivity of CO2 in air at 25 degrees C.
   Units: `m^2 / s`

 - `Dwater`: Diffusivity of CO2 in water at measurement temperature.
   Units: `m^2 / s`

 - `effective_porosity`: Effective porosity of the cell wall.
   Units: `dimensionless`

 - `gcw`: Cell wall conductance to CO2 diffusion expressed on a chloroplast area
   basis. Units: `mol / m^2 chloroplast / s / bar`

 - `gcw_leaf`: Cell wall conductance to CO2 diffusion expressed on a leaf area
   basis. Units: `mol / m^2 leaf / s / bar`

 - `gias`: Intercellular airspace conductance to CO2 diffusion expressed on a
   leaf area basis. Units: `mol / m^2 leaf / s / bar`

 - `gmem`: Membrane conductance to CO2 diffusion expressed on a chloroplast area
   basis. Units: `mol / m^2 chloroplast / s / bar`

 - `gmem_bilayer`: Lipid bilayer conductance to CO2 diffusion expressed on a
   chloroplast area basis. Units: `mol / m^2 chloroplast / s / bar`

 - `gmem_leaf`: Membrane conductance to CO2 diffusion expressed on a leaf area
   basis. Units: `mol / m^2 leaf / s / bar`

 - `KCO2`: Partition coefficient for CO2. Units: `dimensionless`

 - `Pmem25`: Membrane permeability to CO2 at 25 degrees C. Units: `m / s`

 - `R`: Ideal gas constant. Units: `J / K / mol`

 - `rhoH`: Product of molar density of water (`rho`) and Henry's constant for
   CO2 in H2O (`H`). Units: `mol / m^3 / bar`

 - `xi`: Gas path tortuousity in the intercellular airspace. Units: `m / m`
