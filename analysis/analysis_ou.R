############################################################
# Sediment community oxygen consumption patterns and drivers  
############################################################
# Description
# Draft plots, generate tables, and conduct analysis for SCOC

# Author: Yen-Ting Chen
# Date of creation: 2023/08/30
# Date of last modification: 2023/08/30

#######################
# Set up environment
#######################
# clean environment
rm(list = ls())

# Load packages
library(dplyr)
library(tidyr)
library(vegan)
library(Hmisc)
library(MuMIn)
library(corrplot)
library(ggplot2)
library(GGally)
library(ggrepel)
library(patchwork)
library(GRSmacrofauna)
library(writexl)

# Load data
load("data/cruise_color.Rdata") 
load("data/env_variables.Rdata")
load("data/env_selected.Rdata")
load("data/ou.Rdata")
load("data/standing stock.Rdata")

# source function
source("source/lm_model.avg.R")

#############################
# 1. Spatial patterns of SCOC 
#############################
# set up data
sp <- 
  env %>% 
  select(Cruise, Station, Depth, DRM) 

# scale Depth and DRM
sp$Depth <- (sp$Depth - mean(sp$Depth)) / sd(sp$Depth)
sp$DRM <- (sp$DRM - mean(sp$DRM)) / sd(sp$DRM)

ou_sp <- left_join(ou_core, sp)

###############
# in situ TOU #
###############
# create full model
insituTOU_sp_fullmodel <- 
  lm(In_situ_TOU ~
       DRM * Depth + 
       DRM * Cruise + 
       Depth * Cruise, 
     data = ou_sp,
     na.action = "na.fail")

insituTOU_sp_full_sum <- summary(insituTOU_sp_fullmodel)
get_lm_summary(insituTOU_sp_full_sum, "table/lm/In.situ.TOU_sp_lm.xlsx")
# plot(insituTOU_sp_fullmodel) # good residual patterns

# dredge; add adj.R2; add nested info.
insituTOU_sp_d <- dredge(insituTOU_sp_fullmodel, extra = c("R^2"))
insituTOU_sp_d <- calculate_adjR2(insituTOU_sp_fullmodel, insituTOU_sp_d)
insituTOU_sp_d$nested <- nested(insituTOU_sp_d)
head(insituTOU_sp_d)
# store dredge results
write_xlsx(list(dredge = insituTOU_sp_d), 
           path = "table/dredge/in.situ.TOU_sp_d.xlsx")

# extract model with delta < 6 and un-nested model
insituTOU_sp_subset <- get.models(insituTOU_sp_d, subset = delta < 6 & !nested(insituTOU_sp_d))
insituTOU_sp_avg <- model.avg(insituTOU_sp_subset)
get_model.avg_results(insituTOU_sp_avg, "table/model.avg/in.situ.TOU_sp_model.avg.xlsx")


##################################
# 2. Environmental drivers of SCOC 
##################################
# scaling selected env variables
env_scaled <- env_selected[c("Cruise", "Station", env_variables_selected)]
env_scaled[env_variables_selected] <- scale(env_scaled[env_variables_selected])
# Join OU data with scaled env
ou_env <- 
  ou_core %>% 
  select(-Habitat, -Date) %>% 
  left_join(env_scaled)

#############
# In situ TOU
#############
InsituTOU_env_fullmodel <- 
  lm(In_situ_TOU ~ 
     Temperature + 
     Fluorescence +
     D50 +
     TOC +
     CN +
     Chla +
     Porosity,
   data = ou_env,
   na.action = "na.fail")

# 
# plot(InsituTOU_env_fullmodel)
InsituTOU_env_fullmodel_sum <- summary(InsituTOU_env_fullmodel)
get_lm_summary(InsituTOU_env_fullmodel_sum, "table/lm/In.situ.TOU_env_lm.xlsx")

# dredge
InsituTOU_env_d <- dredge(InsituTOU_env_fullmodel, extra = "R^2")
InsituTOU_env_d <- calculate_adjR2(InsituTOU_env_fullmodel, InsituTOU_env_d)
InsituTOU_env_d$nested <- nested(InsituTOU_env_d)
write_xlsx(list(dredge = InsituTOU_env_d), 
           "table/dredge/In.situ.TOU_env_d.xlsx")

# model avg
InsituTOU_env_subset <- get.models(InsituTOU_env_d, subset = delta < 6 & !nested(InsituTOU_env_d))
InsituTOU_env_avg <- model.avg(InsituTOU_env_subset)
get_model.avg_results(InsituTOU_env_avg, "table/model.avg/In.situ.TOU_env_model.avg.xlsx")

##########
# T25 TOU
##########
T25TOU_env_fullmodel <- 
  lm(T25_TOU ~ 
       Fluorescence +
       D50 +
       TOC +
       CN +
       Chla +
       Porosity,
     data = ou_env,
     na.action = "na.fail")

# 
# plot(T25TOU_env_fullmodel)
T25TOU_env_fullmodel_sum <- summary(T25TOU_env_fullmodel)
get_lm_summary(T25TOU_env_fullmodel_sum, "table/lm/T25.TOU_env_lm.xlsx")

# dredge
T25TOU_env_d <- dredge(T25TOU_env_fullmodel, extra = "R^2")
T25TOU_env_d <- calculate_adjR2(T25TOU_env_fullmodel, T25TOU_env_d)
T25TOU_env_d$nested <- nested(T25TOU_env_d)
write_xlsx(list(dredge = T25TOU_env_d), 
           "table/dredge/T25.TOU_env_d.xlsx")

# model avg
T25TOU_env_subset <- get.models(T25TOU_env_d, subset = delta < 6 & !nested(T25TOU_env_d))
T25TOU_env_avg <- model.avg(T25TOU_env_subset)
get_model.avg_results(T25TOU_env_avg, "table/model.avg/T25.TOU_env_model.avg.xlsx")

#########################
# 3. Casual relationships
#########################
########################
# core_level correlation
########################
ss_ou_core_data <- ss_core %>% left_join(ou_core) %>% select(-Date) %>% left_join(env, c("Cruise", "Station"))
ss_ou_station_data <- ss_station %>% left_join(ou_station) %>% left_join(env)

# casual spatial relationships with abundance, biomass, and SCOC
# vs. DRM
lm(Abundance ~ DRM, data = ss_ou_core_data[ss_ou_core_data$Cruise == "OR1-1219", ]) %>% summary
lm(Biomass ~ DRM, data = ss_ou_core_data[ss_ou_core_data$Cruise == "OR1-1219", ]) %>% summary
lm(In_situ_TOU ~ DRM, data = ss_ou_core_data[ss_ou_core_data$Cruise == "OR1-1219", ]) %>% summary

lm(Abundance ~ DRM, data = ss_ou_core_data[ss_ou_core_data$Cruise == "OR1-1242", ]) %>% summary
lm(Biomass ~ DRM, data = ss_ou_core_data[ss_ou_core_data$Cruise == "OR1-1242", ]) %>% summary
lm(In_situ_TOU ~ DRM, data = ss_ou_core_data[ss_ou_core_data$Cruise == "OR1-1242", ]) %>% summary

##################################
# Core-level correlation matrix
##################################
ss_ou_core_corr_data <-
  ss_core %>% 
  left_join(ou_core) %>% 
  select(all_of(c("Cruise", "Station",
                  "Abundance", "Biomass", "In_situ_TOU"))) %>%
  left_join(env[c("Cruise", "Station", 
                  env_variables_spatial, 
                  env_variables_selected)])

ss_ou_pairplot <- 
  ggpairs(ss_ou_core_corr_data, 
          columns = c("Abundance", "Biomass", "In_situ_TOU"),
          aes(color = Cruise, fill = Cruise)) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  theme_bw()

ggsave("figure/pairplot/ss_ou_pairplot.png",
       plot = ss_ou_pairplot,
       width = 5,
       height = 4,
       scale = 1)

##################################
# Station-level correlation matrix
##################################
ss_ou_station_corr_data <-
  ss_station %>% 
  left_join(ou_station) %>% 
  select(all_of(c("Cruise", "Station",
                  "Abundance_mean", "Biomass_mean", "In_situ_TOU_mean"))) %>%
  left_join(env[c("Cruise", "Station", 
                  env_variables_spatial, 
                  env_variables_selected)])

ss_ou_pairplot <- 
  ggpairs(ss_ou_station_corr_data, 
          columns = c("Abundance_mean", "Biomass_mean", "In_situ_TOU_mean"),
          aes(color = Cruise, fill = Cruise)) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  theme_bw()

ggsave("figure/pairplot/ss_ou_pairplot.png",
       plot = ss_ou_pairplot,
       width = 5,
       height = 4,
       scale = 1)

ss_ou_sp_env_pairplot <- 
  ggpairs(ss_ou_station_corr_data, 
          columns = c("Abundance_mean", "Biomass_mean", "In_situ_TOU_mean",
                      env_variables_spatial, 
                      env_variables_selected),
          aes(color = Cruise, fill = Cruise)) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  theme_bw()

ggsave("figure/pairplot/ss_ou_sp_env_pairplot.png",
       plot = ss_ou_sp_env_pairplot,
       width = 10,
       height = 6,
       scale = 1.6)

# plots
ss_ou_names <- c("Abundance" = "Abundance~(ind.~m^{-2})",
                 "Biomass" = "Biomass~(g~m^{-2})",
                 "In_situ_TOU" = "TOU~(mmol~m^{-2}~d^{-1})",
                 "OR1-1219" = "OR1-1219",
                 "OR1-1242" = "OR1-1242")
ss_ou_plot_data <-
  ss_core %>% 
  left_join(ou_core) %>% 
  select(all_of(c("Cruise", "Station",
                  "Abundance", "Biomass", "In_situ_TOU"))) %>%
  pivot_longer(cols = c("Abundance", "Biomass", "In_situ_TOU"),
               names_to = "Variable",
               values_to = "Value") %>% 
  group_by(Cruise, Station, Variable) %>% 
  summarise(mean = mean(Value),
            sd = sd(Value)) %>% 
  left_join(env)

plot_ss_ou_env <- function(x){
  ggplot(ss_ou_plot_data, 
         aes(x = .data[[x]], 
             y = mean, 
             ymax = mean + sd, 
             ymin = mean - sd,
             color = Cruise,
             label = Station)) +
    geom_pointrange() +
    geom_text_repel() +
    facet_grid(Variable ~ Cruise, 
               scales = "free",
               labeller = as_labeller(ss_ou_names, label_parsed)) +
    scale_color_manual(values = cruise_color) +
    theme_bw()
}
ss_ou_drm <- 
  plot_ss_ou_env("DRM") + 
  xlab(Distance~to~River~mouth~(km)) +
  ylab("Value")
ss_ou_depth<- 
  plot_ss_ou_env("Depth") + xlab(Depth~(m)) + 
  xlab(Depth~(m)) +
  ylab("Value")

ggsave("figure/scatterplot/DRM_vs_vars.png", plot = ss_ou_drm, width = 8, height = 6)
ggsave("figure/scatterplot/Depth_vs_vars.png", plot = ss_ou_depth, width = 8, height = 6)


save(ss_ou_drm, ss_ou_depth, file = "data/ss_ou_drm_depth.RData")
