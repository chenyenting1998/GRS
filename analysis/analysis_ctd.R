#################################
# Gaoping River-shelf CTD profile
#################################
# Description
# Plot the CTD profile of the GRS

# Author: Yen-Ting Chen
# Date of creation: 2023/07/07
# Date of last modification: 2023/10/24

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
    facet_grid(Month~Station, scales = "free") +
    theme_bw()
}

# load variables
load("data/cruise_color.RData")
load("data/env_variables.RData")

# add function
add_month <- function(x){
  x$Month <- if_else(x$Cruise == "OR1-1219", "March", "October")
  return(x)
}

#####################
# 1. Data preparation
#####################
# convert ctd to long format
ctd_variables <- c("Temperature", "Salinity", "Density", "Oxygen", "Fluorescence", "Transmission")
ctd_month <- add_month(ctd)
ctd_long <-
  ctd_month %>% 
  pivot_longer(cols = all_of(ctd_variables),
               names_to = "Variables",
               values_to = "Values")

# Extract bottom water
bottom_water <- 
  ctd_month %>%
  group_by(Cruise, Station) %>% 
  filter(Pressure == max(Pressure)) %>% 
  pivot_longer(cols = all_of(ctd_variables),
               names_to = "Variables",
               values_to = "Values")

################################
# 2. Plot variable relationships
################################
sal_temp <- plot_ctd_variable(ctd_month, "Salinity", "Temperature")
sal_dens <- plot_ctd_variable(ctd_month, "Salinity", "Density")
sal_trans <- plot_ctd_variable(ctd_month, "Salinity", "Transmission")
sal_fluo <- plot_ctd_variable(ctd_month, "Salinity", "Fluorescence")
temp_Fluo <- plot_ctd_variable(ctd_month, "Temperature", "Fluorescence")

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
    "March" = "March",
    "October" = "October")

ctd_profile <- 
  ggplot(ctd_long) +
  geom_line(aes(x = Pressure, y = Values, color = Station)) +
  geom_point(data = bottom_water, 
             aes(x = Pressure, y = Values, color = Station)) +
  facet_grid(Month ~ factor(Variables, ctd_variables), 
             labeller = as_labeller(ctd_facets, label_parsed),
             scales = "free",
             switch = "y") +
  scale_color_manual(values = station_color) +
  scale_x_continuous(position = "top", trans = "reverse") +
  coord_flip()+
  theme_bw()

save(ctd_profile, ctd_variables, file = "data/ctd.RData")
ggsave("figure/ctd/ctd_profile.png", scale = 1.4, plot = ctd_profile)
