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
load("data/BC_exp.Rdata")

# source functions
source("source/box.cox.chord.R")
source("source/plot_pca.R")
source("source/plot_rda.R")
source("source/plot_scree.R")
source("source/extract_goodness.R")

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
summary(count_pca)

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
# extract goodness species
count_pca_goodness <- extract_goodness(count_pca, "CA")
count_pca_goodness_plot <-
  ggplot(count_pca_goodness) +
  geom_point(aes(x = PC2, 
                 y = Taxon)) +
  xlab("Cummulative variance explained to PC2") +
  theme_bw()
  
ggsave("figure/pca/count_pca_goodness_plot.png", 
       plot = count_pca_goodness_plot,
       width = 6,
       height = 5)

# plot pca plots scaling = 1
count_pca_sc1 <- 
  plot_pca(count_sc1$pca_sites, 
           count_sc1$pca_species, 
           scaling = 1, 
           eig_vector = count_sc1$pca_eig,
           stretch = 1)
# plot pca plots scaling = 2
count_pca_sc2 <- 
  plot_pca(count_sc2$pca_sites, 
           count_sc2$pca_species, 
           scaling = 2, 
           eig_vector = count_sc2$pca_eig, 
           stretch = 1)

# output
ggsave(filename = "figure/pca/count_pca_screeplot.png", 
       plot = count_pca_screeplot, scale = 1.5)
ggsave(filename = "figure/pca/count_pca_sc1.png", 
       plot = count_pca_sc1,
       scale = 1,
       width = 8,
       height = 6)

ggsave(filename = "figure/pca/count_pca_sc2.png", 
       plot = count_pca_sc2, 
       scale = 1,
       width = 8,
       height = 6)


####################################
# 3. PERMANOVA and PERMDISP -- Count
####################################
# sp
env_sp <- env[, c("Cruise", "Station", env_variables_spatial)]
env_sp$Depth <- scale(env_sp$Depth)
env_sp$DRM <- scale(env_sp$DRM)
data_temp <- left_join(count_wide, env_sp)

# permanova
set.seed(14)
count_permanova <- 
  adonis2(count_chord ~ Depth * DRM * Cruise,
          data = data_temp,
          method = "euclidean",
          permutations = 9999)

write_xlsx(list(PERMANOVA = cbind(" " = rownames(count_permanova), count_permanova)),
           path = "table/permanova/count_permanova.xlsx")

################################################
# 4. Canonical redundancy analysis -- Count data
################################################
# match env_selected data.frame with count_wide
env_selected_expand <- 
  left_join(count_wide, env_selected) %>% 
  select(all_of(env_variables_selected))

# null and full RDA model
count_rda_full <- rda(count_chord ~ ., data = as.data.frame(scale(env_selected_expand)))

# backward selection
count_rda_back <- ordistep(count_rda_full, method = "backward")

# summary
summary(count_rda_full)
summary(count_rda_back)

# test difference between full and reduced model
anova(count_rda_full) # the full model is significant
anova(count_rda_back) # the reduced model is significant
anova(count_rda_full, count_rda_back) # no sig. diff.

# rsquared
RsquareAdj(count_rda_full)
RsquareAdj(count_rda_back) # slight reduction in r2

# vif
vif.cca(count_rda_full) # TOC and porosity have high vif
vif.cca(count_rda_back)

# residual plots
ordiresids(count_rda_full)
ordiresids(count_rda_back)

# species goodness
count_rda_for_goodness <- extract_goodness(count_rda_back, "CCA")
count_rda_for_goodness_plot <-
  ggplot(count_rda_for_goodness) +
  geom_point(aes(x = RDA2, 
                 y = Taxon)) +
  xlab("Cummulative variance explained to RDA2") +
  theme_bw()

ggsave("figure/rda/count_rda_for_goodness_plot.png", 
       plot = count_rda_for_goodness_plot,
       width = 6,
       height = 5)

# extract reduced model statistics
set.seed(10)
count_rda_axis <- anova.cca(count_rda_back, by = "axis", permutations = 9999)

# first two rda axises are sig.
count_rda_margin <- anova.cca(count_rda_back, by = "margin", permutations = 9999)
# all variables are sig.
write_xlsx(list(count_rda_axis = cbind(" " = rownames(count_rda_axis), count_rda_axis),
                count_rda_margin = cbind(" " = rownames(count_rda_margin), count_rda_margin)),
           path = "table/rda/count_rda_anova.xlsx")

# 
# scaling = 1 
count_rda_output_sc1 <- 
  get_rda_output(count_rda_back, 
                 count_wide[1:4], 
                 env_variables_abbr,
                 scaling = 1)
count_rda_plot_sc1 <- 
  plot_rda(rda_sites   = count_rda_output_sc1$rda_sites, 
           rda_env     = count_rda_output_sc1$rda_env,
           rda_species = count_rda_output_sc1$rda_species,
           rda_result  = count_rda_back,
           scaling = 1)
ggsave("figure/rda/count_rda_plot_sc1.png", 
       plot = count_rda_plot_sc1, 
       scale = 1,
       width = 8,
       height = 6)

# scaling = 2
count_rda_output_sc2 <- 
  get_rda_output(count_rda_back, 
                 count_wide[1:4], 
                 env_variables_abbr,
                 scaling = 2)
count_rda_plot_sc2 <- 
  plot_rda(rda_sites   = count_rda_output_sc2$rda_sites, 
           rda_env     = count_rda_output_sc2$rda_env,
           rda_species = count_rda_output_sc2$rda_species,
           rda_result  = count_rda_back,
           scaling = 2)
ggsave("figure/rda/count_rda_plot_sc2.png", 
       plot = count_rda_plot_sc2, 
       scale = 1,
       width = 8,
       height = 6)

