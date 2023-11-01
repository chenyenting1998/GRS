###################################################
# Set up a masterfile to run all the analysis codes
###################################################

# Author: Yen-Ting Chen
# Date of creation: 2023/10/24
# Date of last modification: 2023/10/24

#######################
# Set up environment
#######################
# clean environment
rm(list = ls())

# orders matter
Rfiles <- c(
  # setup
  "setup_cruise color.R",
  "setup_taxa rank and color.R",
  "setup_env_variables.R",
  # exploratory analysis
  "explore_macrofauna composition.R",
  "explore_env.R",
  "explore_ou.R",
  # map
  "map.R",
  # analysis
  "analysis_ctd.R",
  "analysis_env.R",
  "analysis_count macrofauna composition.R",
  "analysis_biomass macrofauna composition.R",
  
  
)

for(i in Rfiles){
  source(i)
}