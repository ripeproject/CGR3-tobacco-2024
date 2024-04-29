# Field growth and weather data

## Overview

This directory contains plant growth traits measured from field grown tobacco
and associated weather data from the summer 2022 field season. Values from
`CGR3_field_growth_traits.xlsx` were used to create Figure 6 and Table S3.
Values from `2022_Field_weather_data.csv` were used to create Figure S5.

## Files

- `2022_Field_weather_data.csv`: A table of key weather values measured at the
  University of Illinois Energy farm field site summer 2022.

- `CGR3_field_growth_traits.xlsx`: Tobaaco plant growth measurements from summer
  2022 field season.

## Data dictionary

The columns in `2022_Field_weather_data.csv` are defined as follows:

 - `doy`: The day of year, where January 1 is day 1

 - `hour`: The time of day expressed as an hour value >= 0 and < 24

 - `precip`: The precipitation rate. Units: `mm / hr`

 - `RECORD`: An observation ID number produced by the weather sensor

 - `rh`: The relative humidity expressed as fractional value between 0 and 1.
   Units: `dimensionless`

 - `solar`: The incident photosynthetically active photon flux density.
   Units: `micromol / m^2 / s`

 - `temp`: The air temperature. Units: `degrees C`

 - `time_zone_offset`: Time zone offset relative to UTC

 - `time`: The time expressed as a fractional day of year
   (`time = doy + hour / 24`)

 - `windspeed`: The wind speed. Units: `m / s`

 - `year`: The year
