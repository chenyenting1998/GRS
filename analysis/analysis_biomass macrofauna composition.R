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
library(ggnewscale)
# library(GRSmacrofauna)
library(writexl)

# Load data
load("data/taxa_rank.Rdata")
load("data/taxa_color.Rdata")
load("data/cruise_color.Rdata")
load("data/env.Rdata")
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

add_month <- function(x){
  x$Month <- if_else(x$Cruise == "OR1-1219", "March", "October")
  return(x)
}

ord.rsquare <- 
  function(fullmodel, backward_model){
    data.frame(Model = c("Full model", "Backward selected"),
               R.squared = c(RsquareAdj(biomass_rda_full)$r.squared, 
                             RsquareAdj(biomass_rda_back)$r.squared),
               Adj.R.squared = c(RsquareAdj(biomass_rda_full)$adj.r.squared, 
                                 RsquareAdj(biomass_rda_back)$adj.r.squared))
    
}

# set up left_joining spatial data
env_sp_metadata <- env[c("Month", "Station", "Depth", "DRM")]
point_size_range = c(3, 6.5)
text_size_range = c(1.7, 3.5)
size_depth_breaks <- c(30,50,70,90)
size_drm_breaks = c(15, 20, 25)
stretch = 0.2

######################################################
# 2. Principal component analysis -- Biomass data ####
######################################################
# Box-Cox-chord transformation
biomass_chord <- 
  biomass_wide %>% 
  select(-all_of(c("Station", "Deployment", "Tube", "Month"))) %>% 
  box.cox.chord(biomass_exp)

# tb-PCA
biomass_pca <- rda(biomass_chord)

## get output scaling = 1 ####
biomass_sc1 <- 
  get_pca_output(biomass_pca, 
                 metadata = left_join(biomass_wide[1:4], env_sp_metadata),
                 # metadata = biomass_wide[1:4],
                 scaling = 1)
## get output scaling = 2 ####
biomass_sc2 <- 
  get_pca_output(biomass_pca, 
                 metadata = left_join(biomass_wide[1:4], env_sp_metadata),
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

## plot pca plots scaling = 1 ####
biomass_pca_sc1 <- 
  plot_pca(biomass_sc1$pca_sites, 
           biomass_sc1$pca_species, 
           scaling = 1, 
           eig_vector = biomass_sc1$pca_eig,
           stretch = .2)

## plot pca plots scaling = 2 ####
biomass_pca_sc2 <- 
  plot_pca(biomass_sc2$pca_sites, 
           biomass_sc2$pca_species, 
           scaling = 2, 
           eig_vector = biomass_sc2$pca_eig, 
           stretch = 1)

## plot pca_sc1 depth ####
biomass_pca_sc1_depth <- 
  ggplot() +
  
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # species
  geom_segment(data = biomass_sc1$pca_species[biomass_sc1$pca_species$Show == TRUE,],
               aes(x = 0, 
                   y = 0,
                   xend = PC1 * stretch, 
                   yend = PC2 * stretch),
               size = .4, 
               color = "purple")+
  geom_label(data = biomass_sc1$pca_species[biomass_sc1$pca_species$Show == TRUE,],
             aes(x = PC1 * stretch, 
                 y = PC2 * stretch, 
                 label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +
  # sites
  geom_point(data = biomass_sc1$pca_sites,
             aes(x = PC1,
                 y = PC2,
                 color = Month,
                 size = Depth)) +
  scale_size_binned("Depth (m)", range = point_size_range, breaks = size_depth_breaks) +
  new_scale("size") +
  geom_text(data = biomass_sc1$pca_sites,
            aes(x = PC1,
                y = PC2,
                color = Month,
                size = Depth,
                label = Station),
            color = "white") +
  scale_size_binned("Depth (m)", range = text_size_range, breaks = size_depth_breaks) +
  # change axis label
  xlab(paste0("PC1 (", biomass_sc1$pca_eig[1], "% of the total variance)")) +
  ylab(paste0("PC2 (", biomass_sc1$pca_eig[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

## plot pca sc1 DRM ####
biomass_pca_sc1_drm <- 
  ggplot() +
  
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # species
  geom_segment(data = biomass_sc1$pca_species[biomass_sc1$pca_species$Show == TRUE,],
               aes(x = 0, 
                   y = 0,
                   xend = PC1 * stretch, 
                   yend = PC2 * stretch),
               size = .4, 
               color = "purple")+
  geom_label(data = biomass_sc1$pca_species[biomass_sc1$pca_species$Show == TRUE,],
             aes(x = PC1 * stretch, 
                 y = PC2 * stretch, 
                 label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +
  # sites
  geom_point(data = biomass_sc1$pca_sites,
             aes(x = PC1,
                 y = PC2,
                 color = Month,
                 size = DRM)) +
  scale_size_binned("DRM (km)", range = point_size_range, breaks = size_drm_breaks) +
  new_scale("size") +
  geom_text(data = biomass_sc1$pca_sites,
            aes(x = PC1,
                y = PC2,
                label = Station,
                size = DRM),
            color = "white") +
  scale_size_binned("DRM (km)", range = text_size_range, breaks = size_drm_breaks) +
  # change axis label
  xlab(paste0("PC1 (", biomass_sc1$pca_eig[1], "% of the total variance)")) +
  ylab(paste0("PC2 (", biomass_sc1$pca_eig[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

## output ####
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

ggsave(filename = "figure/pca/biomass_pca_sc1_depth.png", 
       plot = biomass_pca_sc1_depth, 
       scale = 1,
       width = 8,
       height = 6)

ggsave(filename = "figure/pca/biomass_pca_sc1_drm.png", 
       plot = biomass_pca_sc1_drm, 
       scale = 1,
       width = 8,
       height = 6)

###########################################
# 3. PERMANOVA and PERMDISP -- Biomass ####
###########################################
# sp
env_sp <- env[, c("Month", "Station", env_variables_spatial)]
env_sp$Depth <- scale(env_sp$Depth)
env_sp$DRM <- scale(env_sp$DRM)
data_temp <- left_join(biomass_wide, env_sp)

## permanova ####
set.seed(100)
biomass_permanova <- 
  adonis2(biomass_chord ~ Depth * Month + DRM * Month + Depth * DRM,
          data = data_temp,
          method = "euclidean",
          permutations = 9999)

write_xlsx(list(PERMANOVA = cbind(" " = rownames(biomass_permanova), biomass_permanova)),
           path = "table/permanova/biomass_permanova.xlsx")

#######################################################
# 4. Canonical redundancy analysis -- Biomass data ####
#######################################################
# match env_selected data.frame with biomass_wide
env_selected_expand <- 
  left_join(biomass_wide, env_selected) %>% 
  select(all_of(env_variables_selected))

# full RDA model
biomass_rda_full <- rda(biomass_chord ~ ., data = as.data.frame(scale(env_selected_expand)))

# backward selection
biomass_rda_back <- ordistep(biomass_rda_full, method = "backward", permutations = 9999)

# compare full and reduced model
anova(biomass_rda_full) # significant
anova(biomass_rda_back) # significant
anova(biomass_rda_full, biomass_rda_back) # no sig. diff. btw the full and reduced

# rsquared
RsquareAdj(biomass_rda_full)
RsquareAdj(biomass_rda_back) # reduced r2

# vif
vif.cca(biomass_rda_full) # TOC and porosity have high vif
vif.cca(biomass_rda_back)

# residual plots
ordiresids(biomass_rda_full)
ordiresids(biomass_rda_back)

summary(biomass_rda_back)

# goodness 
biomass_rda_for_goodness <- extract_goodness(biomass_rda_back, "CCA")
biomass_rda_for_goodness_plot <-
  ggplot(biomass_rda_for_goodness) +
  geom_point(aes(x = RDA5, 
                 y = Taxon)) +
  xlab("Cummulative variance explained to RDA2") +
  theme_bw()

ggsave("figure/rda/biomass_rda_goodness_plot.png", 
       plot = biomass_rda_for_goodness_plot,
       width = 6,
       height = 5)

# extract reduced model statistics
set.seed(10)
biomass_rda_axis <- anova.cca(biomass_rda_back, by = "axis", permutations = 9999)
# first two rda axises are sig.
biomass_rda_margin <- anova.cca(biomass_rda_back, by = "margin", permutations = 9999)
# all variables are sig.
write_xlsx(list(biomass_rda_statistics = ord.rsquare(biomass_rda_full, biomass_rda_back),
                biomass_rda_axis = cbind(" " = rownames(biomass_rda_axis), biomass_rda_axis),
                biomass_rda_margin = cbind(" " = rownames(biomass_rda_margin), biomass_rda_margin)),
           path = "table/rda/biomass_rda_anova.xlsx")

## scaling = 1 ####
biomass_rda_output_sc1 <- 
  get_rda_output(biomass_rda_back, 
                 # left_join(biomass_wide[1:4], env_sp),
                 left_join(biomass_wide[1:4], env_sp_metadata),
                 env_variables_abbr,
                 scaling = 1)

biomass_rda_plot_sc1 <- 
  plot_rda(rda_sites = biomass_rda_output_sc1$rda_sites, 
           rda_env = biomass_rda_output_sc1$rda_env,
           rda_species = biomass_rda_output_sc1$rda_species,
           rda_result = biomass_rda_back,
           scaling = 1,
           stretch = 0.2)

ggsave("figure/rda/biomass_rda_plot_sc1.png", 
       plot = biomass_rda_plot_sc1, 
       scale = 1,
       width = 8,
       height = 6)

## biomass_rda_plot_sc1_depth ####
biomass_rda_plot_sc1_depth <- 
  ggplot() +
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # # plot env
  geom_segment(data = biomass_rda_output_sc1$rda_env,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 0),
               size = .5, color = "black")+
  geom_text(data = biomass_rda_output_sc1$rda_env,
            aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = abbr),
            size = 3,
            parse = TRUE) +
  # species vectors
  geom_segment(data = biomass_rda_output_sc1$rda_species[biomass_rda_output_sc1$rda_species$Show == TRUE,],
               aes(x = 0, y = 0, xend = RDA1 * stretch, yend = RDA2 * stretch),
               arrow = arrow(angle = 0),
               size = .5, 
               color = "purple")+
  geom_label(data = biomass_rda_output_sc1$rda_species[biomass_rda_output_sc1$rda_species$Show == TRUE,], 
             aes(x = RDA1 * stretch, y = RDA2 * stretch, label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +   
  # plot stations
  geom_point(data = biomass_rda_output_sc1$rda_sites,
            aes(x = RDA1,
                y = RDA2,
                color = Month,
                size = Depth)) +
  scale_size_binned("Depth (m)", range = point_size_range, breaks = size_depth_breaks) +
  new_scale("size") +
  geom_text(data = biomass_rda_output_sc1$rda_sites,
            aes(x = RDA1,
                y = RDA2,
                label = Station,
                size = Depth),
            color = "white") +
  scale_size_binned("Depth (m)", range = text_size_range, breaks = size_depth_breaks) +
  
  # add R2
  xlab(paste0("RDA1 (", rda_eig_percent(biomass_rda_back)[1], "% of the total variance)")) +
  ylab(paste0("RDA2 (", rda_eig_percent(biomass_rda_back)[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

ggsave("figure/rda/biomass_rda_plot_sc1_depth.png", 
       plot = biomass_rda_plot_sc1_depth, 
       scale = 1,
       width = 8,
       height = 6)

## biomass_rda_plot_sc1_drm ####
biomass_rda_plot_sc1_drm <- 
  ggplot() +
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # # plot env
  geom_segment(data = biomass_rda_output_sc1$rda_env,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 0),
               size = .5, color = "black")+
  geom_text(data = biomass_rda_output_sc1$rda_env,
            aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = abbr),
            size = 3,
            parse = TRUE) +
  # species vectors
  geom_segment(data = biomass_rda_output_sc1$rda_species[biomass_rda_output_sc1$rda_species$Show == TRUE,],
               aes(x = 0, y = 0, xend = RDA1 * stretch, yend = RDA2 * stretch),
               arrow = arrow(angle = 0),
               size = .5, 
               color = "purple")+
  geom_label(data = biomass_rda_output_sc1$rda_species[biomass_rda_output_sc1$rda_species$Show == TRUE,], 
             aes(x = RDA1 * stretch, y = RDA2 * stretch, label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +   
  # plot stations
  geom_point(data = biomass_rda_output_sc1$rda_sites,
             aes(x = RDA1,
                 y = RDA2,
                 color = Month,
                 size = DRM)) +
  scale_size_binned("DRM (km)", range = point_size_range, breaks = size_drm_breaks) +
  new_scale("size") +
  geom_text(data = biomass_rda_output_sc1$rda_sites,
            aes(x = RDA1,
                y = RDA2,
                label = Station,
                size = DRM),
            color = "white") +
  scale_size_binned("DRM (km)", range = text_size_range, breaks = size_drm_breaks) +
  
  # add R2
  xlab(paste0("RDA1 (", rda_eig_percent(biomass_rda_back)[1], "% of the total variance)")) +
  ylab(paste0("RDA2 (", rda_eig_percent(biomass_rda_back)[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

ggsave("figure/rda/biomass_rda_plot_sc1_drm.png", 
       plot = biomass_rda_plot_sc1_drm, 
       scale = 1,
       width = 8,
       height = 6)

## scaling = 2 ####
biomass_rda_output_sc2 <- 
  get_rda_output(biomass_rda_back, 
                 biomass_wide[1:4], 
                 env_variables_abbr,
                 scaling = 2)
biomass_rda_plot_sc2 <- 
  plot_rda(rda_sites = biomass_rda_output_sc2$rda_sites, 
           rda_env = biomass_rda_output_sc2$rda_env,
           rda_species = biomass_rda_output_sc2$rda_species,
           rda_result = biomass_rda_back,
           scaling = 2)
ggsave("figure/rda/biomass_rda_plot_sc2.png", 
       plot = biomass_rda_plot_sc2, 
       scale = 1,
       width = 8,
       height = 6)

save(biomass_pca_sc1, biomass_pca_sc1_drm, biomass_pca_sc1_depth,
     biomass_rda_plot_sc1, stretch,
     biomass_rda_plot_sc1_depth, biomass_rda_plot_sc1_drm, file = "data/biomass_ord_sc1.RData")
