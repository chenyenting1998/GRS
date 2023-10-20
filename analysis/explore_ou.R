#########################################################
# Explore the Gaoping River-shelf oxygen utilization data
#########################################################
# Description
# Explores the GRS oxygen utilization data 

# Author: Yen-Ting Chen
# Date of creation: 2023/07/15
# Date of last modification: 2023/07/15

#######################
# Set up environment
#######################
# clean environment
rm(list = ls())

# Load packages
library(dplyr) # data manipulation
library(tidyr)
library(vegan) # data analysis
library(rstatix) # pipe-friendly stats package
library(ggplot2) # visualization
library(GGally)
library(ggrepel)
library(writexl) # xlsx
library(GRSmacrofauna)  # data package

# load variables
load("data/cruise_color.RData")
load("data/env_variables.RData")

add_month <- function(x){
  x$Month <- if_else(x$Cruise == "OR1-1219", "March", "October")
  return(x)
}

###############################
# 1. Oxygen utilization by core
###############################

range(tou$Incub_temperature) # 20.79450 26.08633
range(tou$In_situ_temperature) # 23.68867 28.14290
mean(tou$Incub_temperature) # 23.33776
mean(tou$In_situ_temperature) # 25.89264

# q10 equation: q10 = (r2/r1)^(10/(t2-t1)); t2 > t1
# rearrange 1: r1 = r2 / (q10^[(t2-t1)/10])
# rearrange 2: r2 = r1 * q10^(t2-t1)/10)
T25 <- 25
q10 <- 2

# The values varied a little when t2 < t1 
tou_t25 <- 
  tou %>% 
  mutate(T25_TOU = if_else(T25 > Incub_temperature,
                           q10 ^ ((T25 - Incub_temperature)/10) * Incub_TOU,
                           Incub_TOU / q10 ^ ((Incub_temperature - T25)/10)))

dou_core <- 
  dou %>% 
  group_by(Cruise, Station, Deployment, Tube) %>% 
  mutate(T25_DOU = if_else(T25 > Incub_temperature,
                           q10 ^ ((T25 - Incub_temperature)/10) * Incub_DOU,
                           Incub_DOU / q10 ^ ((Incub_temperature - T25)/10))) %>% 
  summarise(Electrodes = n(),
            Incub_DOU_mean = mean(Incub_DOU),
            Incub_DOU_sd = sd(Incub_DOU),
            In_situ_DOU_mean = mean(In_situ_DOU),
            In_situ_DOU_sd = sd(In_situ_DOU),
            T25_DOU_mean = mean(T25_DOU),
            T25_DOU_sd = sd(T25_DOU),
            OPD_mean = mean(OPD),
            OPD_sd = sd(OPD)) %>% 
  ungroup()

ou_core <-
  tou_t25 %>% 
  left_join(dou_core) %>% 
  relocate(Date, .before = Station) %>% 
  relocate(In_situ_temperature, .after = Incub_temperature) %>%  
  relocate(Electrodes, .after = Tube)

# TOU
ggplot(ou_core, aes(x = Station)) +
  geom_point(aes(y = Incub_TOU), color = "red") +
  geom_point(aes(y = T25_TOU), color = "green") +
  geom_point(aes(y = In_situ_TOU), color = "blue") +
  facet_grid(~Cruise, scales = "free") +
  ylab("TOU") +
  theme_bw()

# DOU
ggplot(ou_core, aes(x = Station)) +
  geom_point(aes(y = Incub_DOU_mean), color = "red") +
  geom_point(aes(y = T25_DOU_mean), color = "green") +
  geom_point(aes(y = In_situ_DOU_mean), color = "blue") +
  facet_grid(~Cruise, scales = "free") +
  ylab("TOU") +
  theme_bw()

##################################
# 2. Oxygen utilization by station
##################################
dou_station <- 
  dou %>% 
  group_by(Cruise, Station) %>% 
  summarise(Electrodes = n(),
            Incub_DOU_mean = mean(In_situ_DOU),
            Incub_DOU_sd = sd(In_situ_DOU),
            In_situ_DOU_mean = mean(In_situ_DOU),
            In_situ_DOU_sd = sd(In_situ_DOU),
            OPD_mean = mean(OPD),
            OPD_sd = sd(OPD)) %>% 
  add_month()

ou_station <-
  tou_t25 %>% 
  group_by(Cruise, Station) %>% 
  summarise(In_situ_TOU_mean = mean(In_situ_TOU),
            In_situ_TOU_sd = sd(In_situ_TOU),
            T25_TOU_mean = mean(T25_TOU),
            T25_TOU_sd = sd(T25_TOU),
            Incub_TOU_mean = mean(Incub_TOU),
            Incub_TOU_sd = sd(Incub_TOU)) %>% 
  left_join(dou_station, 
            by = c("Cruise", "Station")) %>% 
  add_month

# station vs. TOU====
plot_station_TOU <- function(object, mean, sd){
  object %>% 
    ggplot() +
    geom_pointrange(aes(x = Station, 
                        y = .data[[mean]], 
                        ymax = .data[[mean]] + .data[[sd]], 
                        ymin = .data[[mean]] - .data[[sd]],
                        color = Month),
                    position = position_dodge(width = 0.5)) +
    scale_color_manual(values = cruise_color) +
    theme_bw()
}

plot_station_TOU(ou_station, "In_situ_TOU_mean", "In_situ_TOU_sd")
plot_station_TOU(ou_station, "Incub_TOU_mean", "Incub_TOU_sd")
plot_station_TOU(ou_station, "T25_TOU_mean", "T25_TOU_sd")

# DRM vs. TOU====
plot_env_TOU <- function(object, mean, sd, x){
  object %>% 
    left_join(env) %>% 
    ggplot() +
    geom_pointrange(aes(x = .data[[x]], 
                        y = .data[[mean]], 
                        ymax = .data[[mean]] + .data[[sd]], 
                        ymin = .data[[mean]] - .data[[sd]],
                        color = Month)) +
    geom_text_repel(aes(x = .data[[x]], 
                        y = .data[[mean]], 
                        color = Month,
                        label = Station)) +
    scale_color_manual(values = cruise_color) +
    theme_bw()
}

plot_env_TOU(ou_station, "In_situ_TOU_mean", "In_situ_TOU_sd", "DRM")
plot_env_TOU(ou_station, "Incub_TOU_mean", "Incub_TOU_sd", "DRM")
plot_env_TOU(ou_station, "T25_TOU_mean", "T25_TOU_sd", "DRM")

# Depth vs. TOU ====
plot_env_TOU(ou_station, "In_situ_TOU_mean", "In_situ_TOU_sd", "Depth")
plot_env_TOU(ou_station, "Incub_TOU_mean", "Incub_TOU_sd", "Depth")
plot_env_TOU(ou_station, "T25_TOU_mean", "T25_TOU_sd", "Depth")

# Temperature vs. TOU
plot_env_TOU(ou_station, "In_situ_TOU_mean", "In_situ_TOU_sd", "Temperature")
plot_env_TOU(ou_station, "Incub_TOU_mean", "Incub_TOU_sd", "Temperature")
plot_env_TOU(ou_station, "T25_TOU_mean", "T25_TOU_sd", "Temperature")

########
# 2. DOU
########
ou_station %>% 
  ggplot() +
  geom_pointrange(aes(x = Station, 
                      y = In_situ_DOU_mean, 
                      ymax = In_situ_DOU_mean + In_situ_DOU_sd, 
                      ymin = In_situ_DOU_mean - In_situ_DOU_sd,
                      color = Cruise)) +
  scale_color_manual(values = cruise_color) +
  theme_bw()



save(ou_station, ou_core, file = "data/ou.Rdata")
