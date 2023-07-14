#########################
# Macrofauna composition
#########################
# Description
# Macrofauna assemblage data exploration

# Author: Yen-Ting Chen
# Date of creation: 2023/07/05
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

# source functions
source("analysis/box.cox.chord.R")
source("analysis/plot_pca.R")
source("analysis/plot_rda.R")
source("analysis/plot_scree.R")

# Set up variables
count_exp <- 0.5

###############################################
# 1. Principal component analysis -- Count data
###############################################
# Box-Cox-chord transformation
count_chord <- 
  count_wide %>% 
  select(-all_of(c("Cruise", "Station", "Deployment", "Tube"))) %>% 
  box.cox.chord(count_exp)

# tb-PCA
count_pca <- rda(count_chord)

# get output scaling = 1
count_sc1 <- 
  get_pca_output(count_pca, 
                 metadata = count_wide[1:2],
                 scaling = 1)
# get output scaling = 2
count_sc2 <- 
  get_pca_output(count_pca, 
                 metadata = count_wide[1:2], 
                 scaling = 2)

# plot PC eigenvalues
count_pca_screeplot <- plot_scree(count_pca)

# plot pca plots scaling = 1
count_pca_sc1 <- 
  plot_pca(count_sc1$pca_sites, 
           count_sc1$pca_species, 
           scaling = 1, 
           eig_vector = count_sc1$pca_eig)
# plot pca plots scaling = 2
count_pca_sc2 <- 
  plot_pca(count_sc2$pca_sites, 
           count_sc2$pca_species, 
           scaling = 2, 
           eig_vector = count_sc2$pca_eig, 
           stretch = 2)

# output
ggsave(filename = "figure/polished/count_pca_screeplot.png", plot = count_pca_screeplot, scale = 1.5)
ggsave(filename = "figure/polished/count_pca_sc1.png", plot = count_pca_sc1, scale = 1.5)
ggsave(filename = "figure/polished/count_pca_sc2.png", plot = count_pca_sc2, scale = 1.5)

####################################
# 3. PERMANOVA and PERMDISP -- Count
####################################
# count
set.seed(100)
count_permanova <- 
  adonis2(count_chord ~ Cruise / Station,
          data = count_wide[,c("Cruise", "Station")],
          method = "euclidean",
          permutations = 9999)
count_disp <- 
  betadisper(vegdist(count_chord, method = "euclidean"),
             group = count_wide$Cruise)
count_permdisp <- 
  permutest(count_disp, nperm = 9999)

write_xlsx(list(PERMANOVA = as.data.frame(count_permanova),
                PERMDISP = as.data.frame(count_permdisp$tab)),
           path = "table/count_permanova.xlsx")

################################################
# 4. Canonical redundancy analysis -- Count data
################################################
# match env_selected data.frame with count_wide
env_selected_expand <- 
  left_join(count_wide, env_selected) %>% 
  select(all_of(env_variables_selected))

# full RDA model
count_rda <- rda(count_chord ~ ., data = as.data.frame(scale(env_selected_expand)))
# backward selection
count_rda_back <- ordistep(count_rda, method = "backward")
# compare full and reduced model
anova(count_rda, count_rda_back) # no sig. diff. btw the full and reduced
# rsquared
RsquareAdj(count_rda)
RsquareAdj(count_rda_back) # reduced r2
# vif
vif.cca(count_rda) # TOC and porosity have high vif
vif.cca(count_rda_back)
# residual plots
ordiresids(count_rda)
ordiresids(count_rda_back)

# extract reduced model statistics
set.seed(10)
count_rda_axis <- anova.cca(count_rda_back, by = "axis", permutations = 9999)
# first two rda axises are sig.
count_rda_margin <- anova.cca(count_rda_back, by = "margin", permutations = 9999)
# all variables are sig.
write_xlsx(list(count_rda_axis = count_rda_axis,
                count_rda_margin = count_rda_margin),
           path = "table/count_rda_anova.xlsx")

# 
count_rda_output <- get_rda_output(count_rda_back, count_wide[1:4], env_variables_abbr)
count_rda_plot<- 
  plot_rda_sc1(count_rda_output$rda_sites,
               count_rda_output$rda_env,
               count_rda_back)
ggsave("figure/polished/count_rda_plot.png", plot = count_rda_plot, scale = 1.5)

# Variance partitioning ------------
# match env_spatial dataframe wtih count_wide
env_spatial_expand <- 
  left_join(count_wide, env) %>% 
  select(all_of(env_variables_spatial))
# RDA
count_rda_spatial <- rda(count_chord ~ scale(env_spatial_expand))

# r.square
RsquareAdj(count_rda_spatial)
# backward selection
count_rda_spatial_back <- ordistep(count_rda_spatial, method = "backward")
# anova.cca
set.seed(10)
count_rda_spatial_axis <- anova(count_rda_spatial, by = "axis", permutations = 9999)
count_rda_spatial_margin <- anova(count_rda_spatial, by = "margin", permutations = 9999)

# plot_RDA
count_varpart <- 
  varpart(count_chord,
          ~ scale(env_spatial_expand),
          ~ scale(env_selected_expand))

showvarparts(2, bg = c("hotpink","skyblue"))
plot(count_varpart, bg = c("hotpink","skyblue"))

