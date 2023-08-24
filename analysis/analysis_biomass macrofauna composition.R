#########################
# Macrofauna composition
#########################
# Description
# Macrofauna assemblage data exploration

# Author: Yen-Ting Chen
# Date of creation: 2023/07/14
# Date of last modification: 2023/08/10

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

# species goodness
biomass_pca_goodness <- extract_goodness(biomass_pca, "CA")
biomass_pca_goodness_plot <-
  ggplot(biomass_pca_goodness) +
  geom_point(aes(x = PC2, 
                 y = Taxon)) +
  xlab("Cummulative variance explained to PC2") +
  theme_bw()

ggsave("figure/pca/biomass_pca_goodness_plot.png", 
       plot = biomass_pca_goodness_plot,
       width = 6,
       height = 5)

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
ggsave(filename = "figure/pca/biomass_pca_sc1.png", 
       plot = biomass_pca_sc1, 
       scale = 1,
       width = 8,
       height = 6)

ggsave(filename = "figure/pca/biomass_pca_sc2.png", 
       plot = biomass_pca_sc2, 
       scale = 1,
       width = 8,
       height = 6)


######################################
# 3. PERMANOVA and PERMDISP -- Biomass
######################################
# sp
env_sp <- env[, c("Cruise", "Station", env_variables_spatial)]
env_sp$Depth <- scale(env_sp$Depth)
env_sp$DRM <- scale(env_sp$DRM)
data_temp <- left_join(biomass_wide, env_sp)

# permanova
set.seed(100)
biomass_permanova <- 
  adonis2(biomass_chord ~ Depth*Cruise*DRM,
          data = data_temp,
          method = "euclidean",
          permutations = 9999)

write_xlsx(list(PERMANOVA = cbind(" " = rownames(biomass_permanova), biomass_permanova)),
           path = "table/permanova/biomass_permanova.xlsx")

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
biomass_rda_for <- ordistep(biomass_rda, method = "forward")

# compare full and reduced model
anova(biomass_rda)
anova(biomass_rda_for)
anova(biomass_rda, biomass_rda_for) # no sig. diff. btw the full and reduced
# rsquared
RsquareAdj(biomass_rda)
RsquareAdj(biomass_rda_for) # reduced r2
# vif
vif.cca(biomass_rda) # TOC and porosity have high vif
vif.cca(biomass_rda_for)
# residual plots
ordiresids(biomass_rda)
ordiresids(biomass_rda_for)

# goodness 
biomass_rda_for_goodness <- extract_goodness(biomass_rda_for, "CCA")
biomass_rda_for_goodness_plot <-
  ggplot(biomass_rda_for_goodness) +
  geom_point(aes(x = RDA2, 
                 y = Taxon)) +
  xlab("Cummulative variance explained to RDA2") +
  theme_bw()

ggsave("figure/rda/biomass_rda_goodness_plot.png", 
       plot = biomass_rda_for_goodness_plot,
       width = 6,
       height = 5)

# extract reduced model statistics
set.seed(10)
biomass_rda_axis <- anova.cca(biomass_rda_for, by = "axis", permutations = 9999)
# first two rda axises are sig.
biomass_rda_margin <- anova.cca(biomass_rda_for, by = "margin", permutations = 9999)
# all variables are sig.
write_xlsx(list(biomass_rda_axis = cbind(rownames(biomass_rda_axis), biomass_rda_axis),
                biomass_rda_margin = cbind(rownames(biomass_rda_margin), biomass_rda_margin)),
           path = "table/rda/biomass_rda_anova.xlsx")

# scaling = 1 
biomass_rda_output_sc1 <- 
  get_rda_output(biomass_rda_for, 
                 biomass_wide[1:4], 
                 env_variables_abbr,
                 scaling = 1)
biomass_rda_plot_sc1 <- 
  plot_rda(rda_sites = biomass_rda_output_sc1$rda_sites, 
           rda_env = biomass_rda_output_sc1$rda_env,
           rda_species = biomass_rda_output_sc1$rda_species,
           rda_result = biomass_rda_for,
           scaling = 1)
ggsave("figure/rda/biomass_rda_plot_sc1.png", 
       plot = biomass_rda_plot_sc1, 
       scale = 1,
       width = 8,
       height = 6)

# scaling = 2
biomass_rda_output_sc2 <- 
  get_rda_output(biomass_rda_for, 
                 biomass_wide[1:4], 
                 env_variables_abbr,
                 scaling = 2)
biomass_rda_plot_sc2 <- 
  plot_rda(rda_sites = biomass_rda_output_sc2$rda_sites, 
           rda_env = biomass_rda_output_sc2$rda_env,
           rda_species = biomass_rda_output_sc2$rda_species,
           rda_result = biomass_rda_for,
           scaling = 2)
ggsave("figure/rda/biomass_rda_plot_sc2.png", 
       plot = biomass_rda_plot_sc2, 
       scale = 1,
       width = 8,
       height = 6)

