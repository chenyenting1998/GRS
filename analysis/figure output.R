#################
# figure output #
#################
rm(list = ls())

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(GGally)
library(patchwork)

# load
load("data/env.RData")
load("data/env_variables.RData")
load("data/env_selected.RData")
load("data/cruise_color.RData")

# source function
source("source/plot_pca.R")
######################
# Environment PCA ####
######################
load("data/env_pca.RData")

ggsave("figure/publish/env_pca_plot.png", 
       plot = env_pca_plot,
       width = 5,
       height = 5,
       scale = 1.25)

##################
# CTD profile ####
##################
load("data/ctd.RData")

ggsave("figure/publish/ctd_profile.png",
       plot = ctd_profile,
       width= 8,
       height = 6,
       scale = 1.2)

#############################
# Macrofauna composition ####
#############################
load("data/composition_figure.RData")

composition_figure <- 
  density_percentage_composition -
  biomass_percentage_composition +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggsave("figure/publish/composition_figure.png",
       plot = composition_figure,
       width= 10,
       height = 5,
       scale = 1)

#####################
# Macrofauna PCA ####
#####################
load("data/count_ord_sc1.RData")
load("data/biomass_ord_sc1.RData")

## original ####
count_pca_sc1 <- 
  count_pca_sc1 +
  theme(legend.position = "none")+
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05))

biomass_pca_sc1 <-
  biomass_pca_sc1 +
  scale_x_continuous(limits = c(-0.25, 0.3))+
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05))


ggsave(filename = "figure/publish/macrofauna_pca_sc1.png", 
       plot = count_pca_sc1 - biomass_pca_sc1+ plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")"),
       scale = 1,
       width = 10,
       height = 5)

## depth ####
count_pca_sc1_depth <- 
  count_pca_sc1_depth +
  theme(legend.position = "none")+
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05))

biomass_pca_sc1_depth <-
  biomass_pca_sc1_depth +
  scale_x_continuous(limits = c(-0.25, 0.3))+
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05))


ggsave(filename = "figure/publish/macrofauna_depth_pca_sc1.png", 
       plot = count_pca_sc1_depth - biomass_pca_sc1_depth+ plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")"),
       scale = 1,
       width = 10,
       height = 5)

## drm ####
count_pca_sc1_drm <- 
  count_pca_sc1_drm +
  theme(legend.position = "none")+
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05))

biomass_pca_sc1_drm <-
  biomass_pca_sc1_drm +
  scale_x_continuous(limits = c(-0.25, 0.3))+
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05))


ggsave(filename = "figure/publish/macrofauna_drm_pca_sc1.png", 
       plot = count_pca_sc1_drm - biomass_pca_sc1_drm+ plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")"),
       scale = 1,
       width = 10,
       height = 5)

#####################
# Macrofauna RDA ####
#####################
count_rda_plot_sc1_adj <- 
  count_rda_plot_sc1 + 
  theme(legend.position = "none")+
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05)) +
  coord_fixed()

biomass_rda_plot_sc1_adj <- 
  biomass_rda_plot_sc1 + 
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05))

ggsave(filename = "figure/publish/macrofauna_rda_sc1.png", 
       plot = count_rda_plot_sc1_adj - biomass_rda_plot_sc1_adj + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")"),
       scale = 1,
       width = 10,
       height = 5)

## depth ####
count_rda_plot_sc1_depth <- 
  count_rda_plot_sc1_depth + 
  theme(legend.position = "none")+
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05)) +
  coord_fixed()

biomass_rda_plot_sc1_depth<- 
  biomass_rda_plot_sc1_depth + 
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05))

ggsave(filename = "figure/publish/macrofauna_rda_sc1_depth.png", 
       plot = count_rda_plot_sc1_depth - biomass_rda_plot_sc1_depth + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")"),
       scale = 1,
       width = 10,
       height = 5)

## drm ####
count_rda_plot_sc1_drm <- 
  count_rda_plot_sc1_drm + 
  theme(legend.position = "none")+
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05)) +
  coord_fixed()

biomass_rda_plot_sc1_drm<- 
  biomass_rda_plot_sc1_drm + 
  scale_x_continuous(expand = expansion(0, .05)) +
  scale_y_continuous(expand = expansion(0, .05))

ggsave(filename = "figure/publish/macrofauna_rda_sc1_drm.png", 
       plot = count_rda_plot_sc1_drm - biomass_rda_plot_sc1_drm + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")"),
       scale = 1,
       width = 10,
       height = 5)

##########################
# Abundance, Biomass, SCOC
##########################
load("data/ss_ou_drm_depth.RData")
ss_ou_dd <- ss_ou_drm + theme(legend.position = "none") + ss_ou_depth + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave("figure/publish/ss_ou_drm_depth.png", ss_ou_dd, width = 14, height = 6)

