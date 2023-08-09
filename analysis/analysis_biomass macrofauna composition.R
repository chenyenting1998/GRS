#########################
# Macrofauna composition
#########################
# Description
# Macrofauna assemblage data exploration

# Author: Yen-Ting Chen
# Date of creation: 2023/07/14
# Date of last modification: 2023/07/14

#####################
# Set up environment
#####################
rm(list = ls())

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(GGally)
library(ggrepel)
library(GRSmacrofauna)
library(writexl)

# Load data
load("data/taxa_rank.Rdata")
load("data/taxa_color.Rdata")
load("data/cruise_color.Rdata")
load("data/env_variables.Rdata")
load("data/env_selected.Rdata")
load("data/wide_data.Rdata")
load("data/BC_exp.Rdata")

# source functions
source("source/box.cox.chord.R")
source("source/plot_pca.R")
source("source/plot_rda.R")
source("source/plot_scree.R")

##################################################
# 1. Seek best fit Box-Cox transformation exponent
##################################################
# Exponents within [0,1] does not yield normal distribution of distances
# Use RDA to fit env variables onto various Box-Cox-chord transformed biomass matrix
# Pick the exponent that yields the highest adjR2, maximizes linear fit
# match env_selected data.frame with biomass_wide
# env_selected_expand <- 
#   left_join(biomass_wide, env_selected) %>% 
#   select(all_of(env_variables_selected))
# 
# # set an empty df
# bc.exp_fit <- NULL
# 
# for(bc.exp in c(0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)){
#   # Box-Cox-chord transformation
#   biomass_chord <- 
#     biomass_wide %>% 
#     select(-all_of(c("Cruise", "Station", "Deployment", "Tube"))) %>% 
#     box.cox.chord(bc.exp)
#   
#   # full RDA model
#   rda <- rda(biomass_chord ~ ., data = as.data.frame(scale(env_selected_expand)))
#   
#   # extract rsquared 
#   r2 <- RsquareAdj(rda)
#   
#   # store info
#   temp <- 
#     data.frame(`Box-Cox exp` = bc.exp,
#                r.squared = r2$r.squared,
#                adj.r.squared = r2$adj.r.squared)
#   # merge
#   bc.exp_fit <- rbind(bc.exp_fit, temp)
# }
# bc.exp_fit
# plot(x = bc.exp_fit$Box.Cox.exp, y = bc.exp_fit$adj.r.squared)


#################################################
# 2. Principal component analysis -- Biomass data
#################################################
# Box-Cox-chord transformation
biomass_chord <- 
  biomass_wide %>% 
  select(-all_of(c("Cruise", "Station", "Deployment", "Tube"))) %>% 
  box.cox.chord(biomass_exp)

# tb-PCA
biomass_pca <- rda(biomass_chord)

# get output scaling = 1
biomass_sc1 <- 
  get_pca_output(biomass_pca, 
                 metadata = biomass_wide[1:2],
                 scaling = 1)
# get output scaling = 2
biomass_sc2 <- 
  get_pca_output(biomass_pca, 
                 metadata = biomass_wide[1:2], 
                 scaling = 2)

# plot PC eigenvalues
biomass_pca_screeplot <- plot_scree(biomass_pca)

# plot pca plots scaling = 1
biomass_pca_sc1 <- 
  plot_pca(biomass_sc1$pca_sites, 
           biomass_sc1$pca_species, 
           scaling = 1, 
           eig_vector = biomass_sc1$pca_eig,
           stretch = .3)
# plot pca plots scaling = 2
biomass_pca_sc2 <- 
  plot_pca(biomass_sc2$pca_sites, 
           biomass_sc2$pca_species, 
           scaling = 2, 
           eig_vector = biomass_sc2$pca_eig, 
           stretch = 1)
# output
ggsave(filename = "figure/polished/biomass_pca_sc1.png", 
       plot = biomass_pca_sc1, 
       scale = 1,
       width = 8,
       height = 6)

ggsave(filename = "figure/polished/biomass_pca_sc2.png", 
       plot = biomass_pca_sc2, 
       scale = 1,
       width = 8,
       height = 6)


######################################
# 3. PERMANOVA and PERMDISP -- Biomass
######################################
set.seed(100)
biomass_permanova <- 
  adonis2(biomass_chord ~ Cruise / Station,
          data = biomass_wide[,c("Cruise", "Station")],
          method = "euclidean",
          permutations = 99999)
biomass_disp <- 
  betadisper(vegdist(biomass_chord, method = "euclidean"),
             group = biomass_wide$Cruise)
biomass_permdisp <- 
  permutest(biomass_disp, nperm = 99999)

write_xlsx(list(PERMANOVA = as.data.frame(biomass_permanova),
                PERMDISP = as.data.frame(biomass_permdisp$tab)),
           path = "table/biomass_permanova.xlsx")

##################################################
# 4. Canonical redundancy analysis -- Biomass data
##################################################
# match env_selected data.frame with biomass_wide
env_selected_expand <- 
  left_join(biomass_wide, env_selected) %>% 
  select(all_of(env_variables_selected))

# full RDA model
biomass_rda <- rda(biomass_chord ~ ., data = as.data.frame(scale(env_selected_expand)))
# backward selection
biomass_rda_back <- ordistep(biomass_rda, method = "backward")
# compare full and reduced model
anova(biomass_rda, biomass_rda_back) # no sig. diff. btw the full and reduced
# rsquared
RsquareAdj(biomass_rda)
RsquareAdj(biomass_rda_back) # reduced r2
# vif
vif.cca(biomass_rda) # TOC and porosity have high vif
vif.cca(biomass_rda_back)
# residual plots
ordiresids(biomass_rda)
ordiresids(biomass_rda_back)

# extract reduced model statistics
set.seed(10)
biomass_rda_axis <- anova.cca(biomass_rda_back, by = "axis", permutations = 9999)
# first two rda axises are sig.
biomass_rda_margin <- anova.cca(biomass_rda_back, by = "margin", permutations = 9999)
# all variables are sig.
write_xlsx(list(biomass_rda_axis = biomass_rda_axis,
                biomass_rda_margin = biomass_rda_margin),
           path = "table/biomass_rda_anova.xlsx")

# 
biomass_rda_output <- get_rda_output(biomass_rda_back, biomass_wide[1:4], env_variables_abbr)
biomass_rda_plot<- 
  plot_rda_sc1(biomass_rda_output$rda_sites,
               biomass_rda_output$rda_env,
               biomass_rda_back)
ggsave("figure/polished/biomass_rda_plot.png", 
       plot = biomass_rda_plot, 
       scale = 1,
       width = 8,
       height = 6)

