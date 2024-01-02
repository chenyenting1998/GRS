#########################
# Macrofauna composition
#########################
# Description
# Macrofauna assemblage data exploration

# Author: Yen-Ting Chen
# Date of creation: 2023/07/05
# Date of last modification: 2023/10/24

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
               R.squared = c(RsquareAdj(fullmodel)$r.squared, 
                             RsquareAdj(backward_model)$r.squared),
               Adj.R.squared = c(RsquareAdj(fullmodel)$adj.r.squared, 
                                 RsquareAdj(backward_model)$adj.r.squared))
    
  }


# set up left_joining spatial data
env_sp_metadata <- env[c("Month", "Station", "Depth", "DRM")]
point_size_range = c(3, 6.5)
text_size_range = c(1.7, 3.5)
size_depth_breaks <- c(30,50,70,90)
size_drm_breaks = c(15, 20, 25)
stretch = 0.2

####################################################
# 1. Principal component analysis -- Count data ####
####################################################
# Box-Cox-chord transformation
count_chord <- 
  count_wide %>% 
  select(-all_of(c("Month", "Station", "Deployment", "Tube"))) %>% 
  box.cox.chord(count_exp)

# tb-PCA
count_pca <- rda(count_chord)
summary(count_pca)

# get output scaling = 1
count_sc1 <- 
  get_pca_output(count_pca, 
                 metadata = left_join(count_wide[1:4], env_sp_metadata),
                 scaling = 1)
# get output scaling = 2
count_sc2 <- 
  get_pca_output(count_pca, 
                 metadata = left_join(count_wide[1:4], env_sp_metadata),
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

## plot pca plots scaling = 1 ####
count_pca_sc1 <- 
  plot_pca(count_sc1$pca_sites, 
           count_sc1$pca_species, 
           scaling = 1, 
           eig_vector = count_sc1$pca_eig,
           stretch = 0.2)

## plot pca plots scaling = 2 ####
count_pca_sc2 <- 
  plot_pca(count_sc2$pca_sites, 
           count_sc2$pca_species, 
           scaling = 2, 
           eig_vector = count_sc2$pca_eig, 
           stretch = 1)

## plot pca_sc1 depth ####
count_pca_sc1_depth <- 
  ggplot() +
  
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # species
  geom_segment(data = count_sc1$pca_species[count_sc1$pca_species$Show == TRUE,],
               aes(x = 0, 
                   y = 0,
                   xend = PC1 * stretch, 
                   yend = PC2 * stretch),
               size = .4, 
               color = "purple")+
  geom_label(data = count_sc1$pca_species[count_sc1$pca_species$Show == TRUE,],
             aes(x = PC1 * stretch, 
                 y = PC2 * stretch, 
                 label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +
  # sites point
  geom_point(data = count_sc1$pca_sites,
             aes(x = PC1,
                 y = PC2,
                 color = Month,
                 size = Depth)) +
  scale_size_binned("Depth (m)", range = point_size_range, breaks = size_depth_breaks) +
  # sites text
  new_scale("size") +
  geom_text(data = count_sc1$pca_sites,
            aes(x = PC1,
                y = PC2,
                color = Month,
                size = Depth,
                label = Station),
            color = "white") +
  scale_size_binned("Depth (m)", range = text_size_range, breaks = size_depth_breaks) +
  # change axis label
  xlab(paste0("PC1 (", count_sc1$pca_eig[1], "% of the total variance)")) +
  ylab(paste0("PC2 (", count_sc1$pca_eig[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

## plot pca sc1 DRM ####
count_pca_sc1_drm <- 
  ggplot() +
  
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # species
  geom_segment(data = count_sc1$pca_species[count_sc1$pca_species$Show == TRUE,],
               aes(x = 0, 
                   y = 0,
                   xend = PC1 * stretch, 
                   yend = PC2 * stretch),
               size = .4, 
               color = "purple")+
  geom_label(data = count_sc1$pca_species[count_sc1$pca_species$Show == TRUE,],
             aes(x = PC1 * stretch, 
                 y = PC2 * stretch, 
                 label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +
  # sites point
  geom_point(data = count_sc1$pca_sites,
             aes(x = PC1,
                 y = PC2,
                 color = Month,
                 size = DRM)) +
  scale_size_binned("DRM (km)", range = point_size_range, breaks = size_drm_breaks) +
  # site text
  new_scale("size") +
  geom_text(data = count_sc1$pca_sites,
            aes(x = PC1,
                y = PC2,
                label = Station,
                size = DRM),
            color = "white") +
  scale_size_binned("DRM (km)", range = text_size_range, breaks = size_drm_breaks) +
  # change axis label
  xlab(paste0("PC1 (", count_sc1$pca_eig[1], "% of the total variance)")) +
  ylab(paste0("PC2 (", count_sc1$pca_eig[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

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

ggsave(filename = "figure/pca/count_pca_sc1_depth.png", 
       plot = count_pca_sc1_depth, 
       scale = 1,
       width = 8,
       height = 6)

ggsave(filename = "figure/pca/count_pca_sc1_drm.png", 
       plot = count_pca_sc1_drm, 
       scale = 1,
       width = 8,
       height = 6)


#########################################
# 3. PERMANOVA and PERMDISP -- Count ####
#########################################
# sp
env_sp <- env[, c("Month", "Station", env_variables_spatial)]
env_sp$Depth <- scale(env_sp$Depth)
env_sp$DRM <- scale(env_sp$DRM)
data_temp <- left_join(count_wide, env_sp)

## permanova ####
set.seed(14)
count_permanova <- 
  adonis2(count_chord ~ Depth * DRM  + Depth * Month + DRM * Month,
          data = data_temp,
          method = "euclidean",
          permutations = 9999)

write_xlsx(list(PERMANOVA = cbind(" " = rownames(count_permanova), count_permanova)),
           path = "table/permanova/count_permanova.xlsx")

#####################################################
# 4. Canonical redundancy analysis -- Count data ####
#####################################################
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
  geom_point(aes(x = RDA5, 
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
write_xlsx(list(count_rda_statistics = ord.rsquare(count_rda_full, count_rda_back),
                count_rda_axis = cbind(" " = rownames(count_rda_axis), count_rda_axis),
                count_rda_margin = cbind(" " = rownames(count_rda_margin), count_rda_margin)),
           path = "table/rda/count_rda_anova.xlsx")

## scaling = 1 ####
count_rda_output_sc1 <- 
  get_rda_output(count_rda_back, 
                 left_join(count_wide[1:4], env_sp_metadata),
                 env_variables_abbr,
                 scaling = 1)
count_rda_plot_sc1 <- 
  plot_rda(rda_sites   = count_rda_output_sc1$rda_sites, 
           rda_env     = count_rda_output_sc1$rda_env,
           rda_species = count_rda_output_sc1$rda_species,
           rda_result  = count_rda_back,
           scaling = 1,
           stretch = .2)
ggsave("figure/rda/count_rda_plot_sc1.png", 
       plot = count_rda_plot_sc1, 
       scale = 1,
       width = 8,
       height = 6)

## scaling = 2 ####
count_rda_output_sc2 <- 
  get_rda_output(count_rda_back, 
                 left_join(count_wide[1:4], env_sp_metadata),
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

## count_rda_plot_sc1_depth ####
count_rda_plot_sc1_depth <- 
  ggplot() +
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # # plot env
  geom_segment(data = count_rda_output_sc1$rda_env,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 0),
               size = .5, color = "black")+
  geom_text(data = count_rda_output_sc1$rda_env,
            aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = abbr),
            size = 3,
            parse = TRUE) +
  # species vectors
  geom_segment(data = count_rda_output_sc1$rda_species[count_rda_output_sc1$rda_species$Show == TRUE,],
               aes(x = 0, y = 0, xend = RDA1 * stretch, yend = RDA2 * stretch),
               arrow = arrow(angle = 0),
               size = .5, 
               color = "purple")+
  geom_label(data = count_rda_output_sc1$rda_species[count_rda_output_sc1$rda_species$Show == TRUE,], 
             aes(x = RDA1 * stretch, y = RDA2 * stretch, label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +   
  # plot stations
  geom_point(data = count_rda_output_sc1$rda_sites,
             aes(x = RDA1,
                 y = RDA2,
                 color = Month,
                 size = Depth)) +
  scale_size_binned("Depth (m)", range = point_size_range, breaks = size_depth_breaks) +
  new_scale("size") +
  geom_text(data = count_rda_output_sc1$rda_sites,
            aes(x = RDA1,
                y = RDA2,
                label = Station,
                size = Depth),
            color = "white") +
  scale_size_binned("Depth (m)", range = text_size_range, breaks = size_depth_breaks) +
  
  # add R2
  xlab(paste0("RDA1 (", rda_eig_percent(count_rda_back)[1], "% of the total variance)")) +
  ylab(paste0("RDA2 (", rda_eig_percent(count_rda_back)[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

ggsave("figure/rda/count_rda_plot_sc1_depth.png", 
       plot = count_rda_plot_sc1_depth, 
       scale = 1,
       width = 8,
       height = 6)

## count_rda_plot_sc1_drm ####
count_rda_plot_sc1_drm <- 
  ggplot() +
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # # plot env
  geom_segment(data = count_rda_output_sc1$rda_env,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 0),
               size = .5, color = "black")+
  geom_text(data = count_rda_output_sc1$rda_env,
            aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = abbr),
            size = 3,
            parse = TRUE) +
  # species vectors
  geom_segment(data = count_rda_output_sc1$rda_species[count_rda_output_sc1$rda_species$Show == TRUE,],
               aes(x = 0, y = 0, xend = RDA1 * stretch, yend = RDA2 * stretch),
               arrow = arrow(angle = 0),
               size = .5, 
               color = "purple")+
  geom_label(data = count_rda_output_sc1$rda_species[count_rda_output_sc1$rda_species$Show == TRUE,], 
             aes(x = RDA1 * stretch, y = RDA2 * stretch, label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +   
  # plot stations
  geom_point(data = count_rda_output_sc1$rda_sites,
             aes(x = RDA1,
                 y = RDA2,
                 color = Month,
                 size = DRM)) +
  scale_size_binned("DRM (km)", range = point_size_range, breaks = size_drm_breaks) +
  new_scale("size") +
  geom_text(data = count_rda_output_sc1$rda_sites,
            aes(x = RDA1,
                y = RDA2,
                label = Station,
                size = DRM),
            color = "white") +
  scale_size_binned("DRM (km)", range = text_size_range, breaks = size_drm_breaks) +
  
  # add R2
  xlab(paste0("RDA1 (", rda_eig_percent(count_rda_back)[1], "% of the total variance)")) +
  ylab(paste0("RDA2 (", rda_eig_percent(count_rda_back)[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

ggsave("figure/rda/count_rda_plot_sc1_drm.png", 
       plot = count_rda_plot_sc1_drm, 
       scale = 1,
       width = 8,
       height = 6)

# pca
save(count_pca_sc1, count_pca_sc2, 
     count_pca_sc1_depth, count_pca_sc1_drm, 
     # rda
     count_rda_plot_sc1,
     count_rda_plot_sc1_drm, count_rda_plot_sc1_depth,
     file = "data/count_ord_sc1.RData")

