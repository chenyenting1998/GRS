#################################
# Gaoping River-shelf CTD profile
#################################
# Description
# Plot the CTD profile of the GRS

# Author: Yen-Ting Chen
# Date of creation: 2023/07/07
# Date of last modification: 2023/07/07

#######################
# Set up environment
#######################
# clean environment
rm(list = ls())

# Load packages
library(dplyr) # data manipulation
library(tidyr)
library(vegan) # data analysis
library(ggplot2) # visualization
library(Polychrome)
library(GRSmacrofauna)  # data package

# load functions
plot_ctd_variable <- function(input, xx, yy) {
  ggplot(input, aes(.data[[xx]], .data[[yy]], color = - Pressure)) +
    geom_point() +
    facet_grid(Cruise~Station, scales = "free") +
    theme_bw()
}

# load variables
load("data/cruise_color.RData")
load("data/env_variables.RData")

#####################
# 1. Data preparation
#####################
# convert ctd to long format
ctd_variables <- c("Temperature", "Salinity", "Density", "Oxygen", "Fluorescence", "Transmission")
ctd_long <-
  ctd %>% 
  pivot_longer(cols = all_of(ctd_variables),
               names_to = "Variables",
               values_to = "Values")
# Extract bottom water
bottom_water <- 
  ctd %>%
  group_by(Cruise, Station) %>% 
  filter(Pressure == max(Pressure)) %>% 
  pivot_longer(cols = all_of(ctd_variables),
               names_to = "Variables",
               values_to = "Values")

################################
# 2. Plot variable relationships
################################
sal_temp <- plot_ctd_variable(ctd, "Salinity", "Temperature")
sal_dens <- plot_ctd_variable(ctd, "Salinity", "Density")
sal_trans <- plot_ctd_variable(ctd, "Salinity", "Transmission")
sal_fluo <- plot_ctd_variable(ctd, "Salinity", "Fluorescence")
temp_Fluo <- plot_ctd_variable(ctd, "Temperature", "Fluorescence")

# output
ggsave("figure/ctd/ctd_sal_temp.png", scale = 1.2, plot = sal_temp)
ggsave("figure/ctd/ctd_sal_dens.png", scale = 1.2, plot = sal_dens)
ggsave("figure/ctd/ctd_sal_trans.png", scale = 1.2, plot = sal_trans)
ggsave("figure/ctd/ctd_sal_fluo.png", scale = 1.2, plot = sal_fluo)
ggsave("figure/ctd/ctd_temp_Fluo.png", scale = 1.2, plot = temp_Fluo)

#################
# 3. Plot profile
#################
# set station color
station_color <- kelly.colors(8)[-1]
names(station_color) <- paste0('S', 1:7)

# plot profile
ctd_facets <- 
  c(env_variables_names,
    "OR1-1219" = "OR1-1219",
    "OR1-1242" = "OR1-1242") 
ctd_profile <- 
  ggplot(ctd_long) +
  geom_line(aes(x = Pressure, y = Values, color = Station)) +
  geom_point(data = bottom_water, 
             aes(x = Pressure, y = Values, color = Station)) +
  facet_grid(Cruise ~ factor(Variables, ctd_variables), 
             labeller = as_labeller(ctd_facets, label_parsed),
             scales = "free") +
  coord_flip()+
  scale_x_reverse() +
  scale_color_manual(values = station_color) +
  theme_bw()

ggsave("figure/ctd/ctd_profile.png", scale = 1.4, plot = ctd_profile)
