##########################
# Define objects for `env`
##########################
# Description
# Define some objects for better visualizations regarding to environmental data

# Author: Yen-Ting Chen
# Date of creation: Unknown
# Date of last modification: 2023/07/06

#######################
# Set up environment
#######################
# clean environment
rm(list = ls())

###################
# 1. Define objects
###################
# define env metadata and variables
env_metadata <- c("Cruise", "Habitat", "Station", "Date", "Latitude", "Longitude")
env_variables <- 
  c("Depth", "DRM", "Temperature", "Salinity", "SigmaTheta", 
    "Density", "Oxygen", "Fluorescence", "Transmission",
    "Sand", "Silt", "Clay", "D50",
    "TOC", "TN", "CN", "delta13C", "Chla", "WC", "Porosity")

# define env variable order
env_variables_order <- 
  c("Depth", "DRM", "Temperature", "Salinity", "SigmaTheta", 
    "Density", "Oxygen", "Fluorescence", "Transmission",
    "Sand", "Silt", "Clay", "D50",
    "TOC", "TN", "CN", "delta13C", "Chla", "WC", "Porosity")

# define env full names and units
env_variables_names <- 
  c("Depth" = "Depth~(m)", 
    "DRM" = "Distance~to~river~mouth~(km)", 
    "Temperature" = "Temperature~({}^o*C)", 
    "Salinity" = "Salinity~(PSU)",
    # "Sigma-Theta", 
    "Density" = "Density~(kg/m^3)",
    "Oxygen" = "Oxygen~(mg/L)", 
    "Fluorescence" = "Fluorescence~(mg/L)",
    "Transmission" = "Transmission~('%')", 
    "D50" = "Median~grain~size~(mu*m)", 
    "Clay" = "Clay~('%')", 
    "Silt" = "Silt~('%')", 
    "Sand" = "Sand~('%')",
    "CN" = "CN~(ratio)", 
    "TOC" = "TOC~('%')", 
    "TN" = "TN~('%')",
    "WC" = "Water~content~('%')",
    "delta13C" = "delta^13*C~({}^o/{}[oo])", 
    "Chla" = "Chlorophyll~a~(ng/g)", 
    "Porosity" = "Porosity~('%')")

# define env abbr
env_variables_abbr <- 
  c("Depth" = "Depth", 
    "DRM" = "DRM", 
    "Temperature" = "Temp", 
    "Salinity" = "Sal",
    # "Sigma-Theta", 
    "Density" = "Den",
    "Oxygen" = "Oxy", 
    "Fluorescence" = "Fluo",
    "Transmission" = "Trans", 
    "D50" = "D50", 
    "Clay" = "Clay", 
    "Silt" = "Silt", 
    "Sand" = "Sand",
    "CN" = "CN", 
    "TOC" = "TOC", 
    "TN" = "TN",
    "WC" = "WC",
    "delta13C" = "delta^13*C", 
    "Chla" = "Chla", 
    "Porosity" = "Por")

###########
# 2. Output
###########
save(env_metadata,
     env_variables, 
     env_variables_abbr, 
     env_variables_names, 
     env_variables_order, 
     file = "data/env_variables.RData")