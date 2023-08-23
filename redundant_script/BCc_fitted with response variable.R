###################################################################
# Seek Box-Cox-chord exponent fits best with the response variables
###################################################################
# Description
# Find a best Box-Cox-chord exponent (BCc exponent) that could best fit the response variable
# This searching method is recommended by Legendre and Bocard (2019) under the condition
# that shows multinormality (which is checked with Dagnelie's test).

# Reference:
# Legendre, P. , Bocard, D. (2018). Box-Cox-chord transformations for community composition 
# data prior to beta diversity analysis. Ecography. 41(11):1820-1824.

# Author: Yen-Ting Chen
# Date of creation: 2023/07/14
# Date of last modification: 2023/08/10

# set up environment
rm(list = ls())
library(dplyr)
library(vegan)
library(GRSmacrofauna)

# load data
load("data/env_selected.Rdata")
load("data/wide_data.Rdata")

# source function
source("source/box.cox.chord.R")

#######################
# 1.  Exhaustive search
#######################
# Exponents within [0,1] does not yield normal distribution of distances
# Use RDA to fit env variables onto various BCc transformed community matrix
# Pick the exponent that yields the highest adjR2, maximizes linear fit
# match env_selected data.frame with biomass_wide
env_selected_expand <-
  left_join(biomass_wide, env_selected) %>%
  select(all_of(env_variables_selected))

# set an empty df
bc.exp_fit <- NULL

for(bc.exp in c(0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)){
  # BCc transformation
  biomass_chord <-
    biomass_wide %>%
    select(-all_of(c("Cruise", "Station", "Deployment", "Tube"))) %>%
    box.cox.chord(bc.exp)

  # full RDA model
  rda <- rda(biomass_chord ~ ., data = as.data.frame(scale(env_selected_expand)))

  # extract rsquared
  r2 <- RsquareAdj(rda)

  # store info
  temp <-
    data.frame(`Box-Cox exp` = bc.exp,
               r.squared = r2$r.squared,
               adj.r.squared = r2$adj.r.squared)
  # merge
  bc.exp_fit <- rbind(bc.exp_fit, temp)
}
bc.exp_fit

# visualize the best BCc exponent across bc.exp
plot(x = bc.exp_fit$Box.Cox.exp, y = bc.exp_fit$adj.r.squared,
     xlab = "Box-Cox-chord exponent",
     ylab = "Adjusted R squared")

##############
# 2. Visualize
##############
# locate exp with the maximum adj.R2 
max_adjr2_which <- which(bc.exp_fit$adj.r.squared == max(bc.exp_fit$adj.r.squared))
# extract bc_exp
max_bc_exp <- bc.exp_fit$Box.Cox.exp[max_adjr2_which]

# BCc transformation
biomass_chord <-
  biomass_wide %>%
  select(-all_of(c("Cruise", "Station", "Deployment", "Tube"))) %>%
  box.cox.chord(max_bc_exp)

# full RDA model
rda <- rda(biomass_chord ~ ., data = as.data.frame(scale(env_selected_expand)))

# check residual
ordiresids(rda) # roughly normally distributed
