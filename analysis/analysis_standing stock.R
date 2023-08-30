##############################################
# Macrofauna abundance, biomass, and diveristy 
##############################################
# Description
# Create plots and tables for macrofauna density and biomass 

# Author: Yen-Ting Chen
# Date of creation: Unknown
# Date of last modification: 2023/07/06

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
load("data/macrofauna composition.Rdata")
load("data/wide_data.Rdata")
load("data/cruise_color.Rdata") 
load("data/env_variables.Rdata")
load("data/env_selected.Rdata")
load("data/ou.Rdata")

# define function
plot_ss <- function(object, x, y, sd, text_repel = FALSE){
  plot <- 
    object %>% 
    ggplot(aes(x = .data[[x]], 
               y = .data[[y]], 
               ymin = .data[[y]] - .data[[sd]], 
               ymax = .data[[y]] + .data[[sd]],
               color = Cruise)) +
    geom_pointrange(position = position_dodge(width = .5))
  
  if(text_repel == TRUE){
    plot <-
      plot + 
      geom_text_repel(data = object,
                      aes(x = .data[[x]],
                          y = .data[[y]],
                          label = Station,
                          color = Cruise),
                      seed = 10)
  }
  plot <- 
    plot +
    scale_color_manual(values = cruise_color) +
    facet_wrap(~Variable, scales = "free") +
    ylab("Value") +
    theme_bw()
  
  return(plot)
}
get_lm_summary <- function(x, path){
  # extract coefs ====
  coefficients <- as.data.frame(x$coefficients)
  coefficients <- cbind(" " = row.names(coefficients), coefficients)
  
  # extract full model statistics ====
  r.squared <- x$r.squared
  adj.r.squared <- x$adj.r.squared
  fstatistic <- x$fstatistic
  df <- x$df
  
  ## calculate p-value ====
  p_value <- pf(fstatistic[1], fstatistic[2], fstatistic[3], lower.tail = FALSE)
  
  # set up a dataframe for full model statistics ====
  fullmodel_statistcs <- 
    data.frame("r squared" = r.squared,
               "adj. r squared" = adj.r.squared,
               "f statistic" = fstatistic[1], 
               "df" = df[2], # get n-p only
               "p value" = p_value)
  write_xlsx(list(Coefficients = coefficients,
                  Fullmodel_statistics = fullmodel_statistcs),
             path = path)
}
calculate_adjR2 <- function(fullmodel, modelset){
  # set model 
  n <- nrow(fullmodel$model)
  k <- modelset$df - 2 # minus parameters of intercept and standard error
  # calculate adj.R2
  modelset$`adj.R^2` <- 1 - ((1 - modelset$`R^2`) * (n - 1) / (n - k - 1))
  return(modelset)
}
get_model.avg_results <- function(avg_model, path){
  # summary() model.avg object
  if(!class(avg_model) %in% "summary.averaging"){
    sum_object <- summary(avg_model)
  }

  # extract objects
  coefficients <- as.data.frame(sum_object$coefficients)
  msTable <- as.data.frame(sum_object$msTable)
  sw <- as.data.frame(sum_object$sw)
  coefmat.full <- as.data.frame(sum_object$coefmat.full)
  coefmat.subset <- as.data.frame(sum_object$coefmat.subset)
  
  # add rownames
  coefficients <- cbind(" " = rownames(coefficients), coefficients)
  msTable <- cbind(" " = rownames(msTable), msTable)
  sw <- cbind(" " = rownames(sw), sw)
  coefmat.full <- cbind(" " = rownames(coefmat.full), coefmat.full)
  coefmat.subset <- cbind(" " = rownames(coefmat.subset), coefmat.subset)
  
  # write xlsx
  write_xlsx(list(coefficients = coefficients, 
                  msTable = msTable, 
                  sw = sw, 
                  coefmat.full = coefmat.full, 
                  coefmat.subset = coefmat.subset),
             path = path)
}

# define multiple corer area; diameter = 10 cm
mc_area <- (0.1 / 2) ^ 2 * pi

###############################################
# 1. Calculate unit area abundance and biomass
###############################################
ss_core <-
  composition %>%
  group_by(Cruise, Station, Deployment, Tube) %>%
  summarise(
    Abundance = sum(Count) / mc_area,
    Biomass = sum(Biomass) / 1000 / mc_area) %>%  # mg to g
  ungroup()

# output .xlsx table
ss_station <-
  ss_core %>% 
  group_by(Cruise, Station) %>%  
  # get the mean and sd of abundance and biomass of each station
  summarise(Abundance_mean = mean(Abundance),
            Abundance_sd = sd(Abundance),
            Biomass_mean = mean(Biomass),
            Biomass_sd = sd(Biomass)) 

write_xlsx(ss_station, "table/data/abundance_biomass.xlsx")

############################################################
# 2. Plot standing stock relationships with station and DRM
############################################################
# set up data.frame for ggplot
ss_plot_data <-
  ss_core %>%
  pivot_longer(cols = c("Abundance", "Biomass"),
              names_to = "Variable",
              values_to = "Value") %>%
  group_by(Cruise, Station, Variable) %>%
  summarise(mean = mean(Value),
            sd = sd(Value)) %>%
  left_join(env) %>% 
  ungroup()
  
# vs station (discret var.)
ss_vs_Station <- plot_ss(ss_plot_data, "Station", "mean", "sd")
ss_vs_DRM <- plot_ss(ss_plot_data, "DRM", "mean", "sd", TRUE) + xlab("Distance to river mouth (km)")
ss_vs_depth <- plot_ss(ss_plot_data, "Depth", "mean", "sd", TRUE) + xlab("Depth (m)")

ggsave("figure/standing_stock_station.png",
       plot = ss_vs_Station,
       scale = 1.2,
       width = 6,
       height = 3)

ggsave("figure/standing_stock_drm_depth.png",
       plot = ss_vs_depth/ss_vs_DRM,
       scale = 1.2,
       width = 8,
       height = 6)

###################################################################
# 3. Statistical relationships between ss_core and DRM*Depth*Cruise
###################################################################
# extract spatiotemporal data.frame
sp <- 
  env %>% 
  select(Cruise, Station, Depth, DRM) 

# scale Depth and DRM
sp$Depth <- (sp$Depth - mean(sp$Depth)) / sd(sp$Depth)
sp$DRM <- (sp$DRM - mean(sp$DRM)) / sd(sp$DRM)

ss_sp <- left_join(ss_core, sp)

# generate a fullmodel regarding spatiotemporal relationships
# Develop hypotheses for variable inputs:
# Single variables:
#    DRM    : Areas closer to the river is prone to sedimentation
#    Depth  : deeper waters can act as a refuge for the benthos
#    Cruise : before and after flood season + a typhoon

# Two-way Interaction:
#    DRM:Depth    : further/deeper waters seem to be less prone to sedimentation
#    Cruise:Depth : shallower waters might be more abu. before flood season.
#    Cruise:DRM   : closer waters might be more abu. before flood season.

# Three-way interaction:
#    Cruise:DRM:Depth : The interaction between water depth and distance varies at flood season.

#############
# Abundance #
#############
abu_sp_fullmodel <- 
  lm(log10(Abundance) ~
       DRM * Depth + 
       DRM * Cruise +
       Depth * Cruise,
     data = ss_sp,
     na.action = "na.fail")

# extract full model coef. and statistics
abu_sp_full_sum <- summary(abu_sp_fullmodel)
get_lm_summary(abu_sp_full_sum, "table/lm/abu_sp_lm.xlsx")
# plot residual patterns
# plot(abu_sp_fullmodel)

# dredge; add adj.R2; add nested info.
abu_sp_d <- dredge(abu_sp_fullmodel, extra = c("R^2"))
abu_sp_d <- calculate_adjR2(abu_sp_fullmodel, abu_sp_d)
abu_sp_d$nested <- nested(abu_sp_d)
# store model selection result 
write_xlsx(list(dredge = abu_sp_d),path = "table/dredge/abu_sp_d.xlsx")

# extract models following Richard 2005, 2008
abu_sp_d_subset <- get.models(abu_sp_d, subset = delta < 6 & !nested(abu_sp_d))
# model averaging
abu_sp_avg <- model.avg(abu_sp_d_subset)
# plot coefficients
plot(abu_sp_avg)
# store averaging result
get_model.avg_results(abu_sp_avg, "table/model.avg/abu_sp_model.avg.xlsx")

###########
# Biomass #
###########
# create full model
bio_sp_fullmodel <- 
  lm(log10(Biomass) ~
       DRM * Depth + 
       DRM * Cruise + 
       Depth * Cruise, 
     data = ss_sp,
     na.action = "na.fail")

bio_sp_full_sum <- summary(bio_sp_fullmodel)
get_lm_summary(bio_sp_full_sum, "table/lm/bio_sp_lm.xlsx")
# plot(bio_sp_fullmodel) # good residual patterns

# dredge; add adj.R2; add nested info.
bio_sp_d <- dredge(bio_sp_fullmodel, extra = c("R^2"))
bio_sp_d <- calculate_adjR2(bio_sp_fullmodel, bio_sp_d)
bio_sp_d$nested <- nested(bio_sp_d)

# store dredge results
write_xlsx(list(dredge = bio_sp_d), 
           path = "table/dredge/bio_sp_d.xlsx")

# extract model with delta < 6
bio_sp_subset <- get.models(bio_sp_d, subset = delta < 6 & !nested(bio_sp_d))
bio_sp_avg <- model.avg(bio_sp_subset)
get_model.avg_results(bio_sp_avg, "table/model.avg/bio_sp_model.avg.xlsx")

###################################################
# 4. Environmental drivers of Abundance and biomass
###################################################
# scaling selected env variables
env_scaled <- env_selected
env_scaled[env_variables_selected] <- scale(env_scaled[env_variables_selected])

ss_env <- 
  left_join(ss_core, env_scaled) %>% 
  select(all_of(c("Cruise", "Station", 
                  "Abundance", "Biomass", 
                  env_variables_selected)))

#############
# Abundance #
#############
# abundance full model
abu_env_fullmodel <- 
  lm(log10(Abundance) ~
       Temperature +
       Fluorescence +
       D50 + 
       TOC +
       CN +
       Chla +
       Porosity,
     data = ss_env,
     na.action = "na.fail")

# full model statistics
abu_env_fullmodel_sum <- summary(abu_env_fullmodel)
get_lm_summary(abu_env_fullmodel_sum, "table/lm/abu_env_lm.xlsx")
# full model residual variations
# plot(abu_env_fullmodel)

# dredge all combinations of parameter inputs
abu_env_d <- dredge(abu_env_fullmodel, extra = "R^2")
abu_env_d <- calculate_adjR2(abu_env_fullmodel, abu_env_d)
abu_env_d$nested <- nested(abu_env_d)
# show top results
head(abu_env_d)
# store dredge results
write_xlsx(list(dredge = abu_env_d), 
           path = "table/dredge/abu_env_d.xlsx")

# extract top models 
abu_env_subset <- get.models(abu_env_d, subset = delta < 6 & !nested(abu_env_d))

# model average
abu_env_avg <- model.avg(abu_env_subset)
abu_env_avg_sum <- summary(abu_env_avg)
get_model.avg_results(abu_env_avg, "table/model.avg/abu_env_model.avg.xlsx")
plot(abu_env_avg)

###########
# Biomass #
###########
# Biomass full model
bio_env_fullmodel <- 
  lm(log10(Biomass) ~
       Temperature +
       Fluorescence +
       D50 + 
       TOC +
       CN +
       Chla +
       Porosity,
     data = ss_env,
     na.action = "na.fail")

# full model statistics
bio_env_fullmodel_sum <- summary(bio_env_fullmodel)
get_lm_summary(bio_env_fullmodel_sum, path = "table/lm/bio_env_lm.xlsx")
# full model residual variations
# plot(bio_env_fullmodel)

# dredge; add adj.R2; add nested info.
bio_env_d <- dredge(bio_env_fullmodel, extra = "R^2")
bio_env_d <- calculate_adjR2(bio_env_fullmodel, bio_env_d)
bio_env_d$nested <- nested(bio_env_d)
# show top results
head(bio_env_d)
# store dredge results
write_xlsx(bio_env_d, path = "table/dredge/bio_env_d.xlsx")

# extract top models
bio_env_subset <- get.models(bio_env_d, subset = delta < 6 & !nested(bio_env_d))
# model average
bio_env_avg <- model.avg(bio_env_subset)
bio_env_avg_sum <- summary(bio_env_avg)
plot(bio_env_avg)
get_model.avg_results(bio_env_avg, "table/model.avg/bio_env_model.avg.xlsx")

