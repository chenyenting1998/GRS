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

#################
# 2. Plot profile
#################
# set station color
station_color <- kelly.colors(8)[-1]
names(station_color) <- paste0('S', 1:7)

# plot profile
ctd_profile <- 
  ggplot(ctd_long) +
  geom_line(aes(x = Pressure, y = Values, color = Station)) +
  geom_point(data = bottom_water, 
             aes(x = Pressure, y = Values, color = Station)) +
  facet_grid(Cruise ~ factor(Variables, ctd_variables), scales = "free") +
  coord_flip()+
  scale_x_reverse() +
  scale_color_manual(values = station_color) +
  theme_bw()

ggsave("figure/ctd_profile.png", plot = ctd_profile)
