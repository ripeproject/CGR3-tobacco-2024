# Define a function for calculating effective porosity (and other values). Here
# we use equations from the following sources:
# - Nobel (2009)                 https://doi.org/10.1016/B978-0-12-374143-1.X0001-4
# - Xiong (2023)                 https://doi.org/10.1111/tpj.16098
# - von Caemmerer & Evans (2015) https://doi.org/10.1111/pce.12449
# See methods section for more information.
calculate_effective_porosity <- function(
  eps_membrane, # membrane conductance enhancement factor (dimensionless)
  Fias,         # fraction of intercellular airspace      (dimensionless)
  gmc,          # mesophyll conductance to CO2            (mol / m^2 leaf / s / bar)
  Sc,           # chloroplast surface area per leaf area  (m^2 chloroplast / m^2 leaf)
  Tcw,          # cell wall thickness                     (nm)
  Tleaf,        # leaf temperature                        (degrees C)
  Tmes          # mesophyll thickness                     (um)
)
{
  # Set constant values
  C_to_K   <- 273     # conversion from degrees C to Kelvin          (K or degrees C) [von Caemmerer & Evans (2015)]
  Dair     <- 1.51e-5 # diffusivity of CO2 in air at 25 degrees C    (m^2 / s)        [see page 1045 in Xiong (2023)]
  Etobacco <- 66.0    # activation energy for membrane conductance   (kJ / mol)       [Table 1 in von Caemmerer & Evans (2015)]
  KCO2     <- 1       # partition coefficient for CO2                (dimensionless)  [see page 400 in Nobel (2009)]
  Pmem25   <- 1.4e-3  # membrane permeability to CO2 at 25 degrees C (m / s)          [see page 631 in von Caemmerer & Evans (2015)]
  R        <- 8.314   # ideal gas constant                           (J / K / mol)    [von Caemmerer & Evans (2015)]
  Tref     <- 25      # reference temperature for some calculations  (degrees C)      [von Caemmerer & Evans (2015)]
  xi       <- 1.57    # gas path tortuosity                          (m / m)          [see page 1045 in Xiong (2023)]

  # Unit conversions
  Etobacco_J <- Etobacco * 1e3 # J / mol
  Tcw_m <- Tcw * 1e-9          # m
  Tleaf_K <- Tleaf + C_to_K    # K
  Tmes_m <- Tmes * 1e-6        # m
  Tref_K  <- 25 + C_to_K       # K

  # Product of molar density of water (rho) and Henry's constant for CO2 in H2O
  # (H) calculated using Equation 3 from von Caemmerer & Evans (2015)
  rhoH <- 33.06 * exp(2400 * (1 / Tleaf_K - 1 / Tref_K)) # mol / m^3 / bar

  # Diffusivity of CO2 in water calculated using Equation 4 from von Caemmerer &
  # Evans (2015)
  Dwater <- 1.81e-6 * exp(-16900 / (R * Tleaf_K)) # m^2 / s

  # Membrane conductance to CO2 calculated using Equation 5 from von Caemmerer &
  # Evans (2015). Here we also include a transport enhancement due to
  # facilitation processes as described on page 1046 of Xiong (2023).
  tleaf_factor <- exp((Tleaf - Tref) * Etobacco_J / (R * Tref_K * Tleaf_K)) # dimensionless
  gmem_bilayer <- rhoH * Pmem25 * tleaf_factor                              # mol / m^2 chloroplast / s / bar

  # Include a membrane transport enhancement due to facilitation processes as
  # described on page 1046 of Xiong (2023)
  gmem <- (1 + eps_membrane) * gmem_bilayer # mol / m^2 chloroplast / s / bar

  # Conductance through the intercellular airspace calculated using an
  # un-numbered equation on page 1045 of Xiong (2023)
  gias_ms <- Dair * Fias / (0.5 * Tmes_m * xi) # m / s

  # Convert gias_ms from a vapor density basis (where conductance is in m / s)
  # to a partial pressure basis (where conductance is in mol / m^2 / s / bar).
  # The conversion factor is R * T, which has units J / mol. Note that 1 J =
  # 1 Pa m^3 = 1e-5 bar m^3.
  conv <- R * Tleaf_K * 1e-5 # m^3 bar / mol
  gias <- gias_ms / conv     # mol / m^2 leaf / s / bar

  # Calculate gcw
  gcw <- 1 / (Sc / gmc - Sc / gias - 1 / gmem) # mol / m^2 chloroplast / s / bar

  # Calculate effective porosity
  effective_porosity <- Tcw_m * gcw * conv / (Dwater * KCO2) # dimensionless

  data.frame(
    conv = conv,                             # m^3 bar / mol
    eps_membrane = eps_membrane,             # dimensionless
    Dair = Dair,                             # m^2 / s
    Dwater = Dwater,                         # m^2 / s
    effective_porosity = effective_porosity, # dimensionless
    gcw = gcw,                               # mol / m^2 chloroplast / s / bar
    gcw_leaf = gcw * Sc,                     # mol / m^2 leaf / s / bar
    gias = gias,                             # mol / m^2 leaf / s / bar
    gmem = gmem,                             # mol / m^2 chloroplast / s / bar
    gmem_bilayer = gmem_bilayer,             # mol / m^2 chloroplast / s / bar
    gmem_leaf = gmem * Sc,                   # mol / m^2 leaf / s / bar
    KCO2 = KCO2,                             # dimensionless
    Pmem25 = Pmem25,                         # m / s
    R = R,                                   # J / K / mol
    rhoH = rhoH,                             # mol / m^3 / bar
    xi = xi                                  # m / m
  )
}

# Estimate a membrane conductance enhancement factor from Figure 7d of
# Xiong (2023), taking the average value from the two Solanaceae crops (Stu and
# Les in their notation)
eps <- 2.0

# Load data
leaf_properties <- read.csv('leaf_properties.csv')

# Calculate effective porosity and other quantities
results <- with(leaf_properties, {
  calculate_effective_porosity(
    eps,
    Fias,
    gmc,
    Sc,
    Tcw,
    Tleaf,
    Tmes
  )
})

results <- cbind(leaf_properties, results)

# Print average effective porosity
cat('\naverage effective porosity\n')
print(tapply(results$effective_porosity, results$Genotype, mean))
cat('\n')

# Print average cell wall conductance
cat('\naverage gcw_leaf\n')
print(tapply(results$gcw_leaf, results$Genotype, mean))
cat('\n')

# Print average membrane conductance
cat('\naverage gmem_leaf\n')
print(tapply(results$gmem_leaf, results$Genotype, mean))
cat('\n')

# Print average IAS conductance
cat('\naverage gias\n')
print(tapply(results$gmem_leaf, results$Genotype, mean))
cat('\n')

# Save results
write.csv(results, file = 'effective_porosity.csv', row.names = FALSE)
