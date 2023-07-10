#########################
# Macrofauna composition
#########################
# Description
# Macrofauna assemblage data exploration

# Author: Yen-Ting Chen
# Date of creation: 2023/07/05
# Date of last modification: 2023/07/08

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
library(GRSmacrofauna)
library(writexl)

# Load data
load("data/taxa_rank.Rdata")
load("data/taxa_color.Rdata")
load("data/cruise_color.Rdata")
load("data/env_variables.Rdata")
# save("data/macrofauna composition.RData")
load("data/wide_data.Rdata")

# Load functions
log10_hell_dist <- function(data){
  dist<- 
    data %>% 
    ungroup() %>%
    select(-Cruise, -Station, -Deployment, -Tube) %>% 
    # log_10(x) + 1 transform
    decostand(method = "log", logbase = 10) %>%
    # Hellinger transformation to deal with double zero
    decostand(method = "hellinger") 
  return(dist)
}

# extract pca data for plotting
get_pca_output <- 
  function(rda_object, 
           metadata, 
           scaling = 1, 
           goodness_threshold = 0.4){
    # get site scores
    pca_sites <-
      scores(rda_object, scaling = scaling)$sites %>%
      as.data.frame() %>% 
      cbind(metadata) 
    
    # get species scores
    pca_species <- 
      scores(rda_object, scaling = scaling)$species %>% 
      as.data.frame()
    pca_species$Taxon <- rownames(pca_species)
    
    # calculate species goodness of fit
    pca_goodness <- 
      goodness(rda_object, choices = 1:2, display = "species", model = "CA") %>% 
      as.data.frame()
    # match goodness to sc2 species scores
    pca_species$goodness <- pca_goodness$PC2[match(rownames(pca_goodness), rownames(pca_species))]
    # use 40% as a cutoff criteria
    pca_species$Show <- pca_species$goodness > goodness_threshold
    
    # extract eigenvalue
    pca_eig <- round(eigenvals(rda_object)/sum(eigenvals(rda_object))*100, 2)
    
    return(list(pca_sites = pca_sites, 
                pca_species = pca_species, 
                pca_goodness = pca_goodness, 
                pca_eig = pca_eig))
  }
# produce screeplot
plot_scree <-
  function(rda_object){
  percent <- round(eigenvals(rda_object)/ sum(eigenvals(rda_object)) * 100, 2) %>% as.vector
  data <- 
    data.frame(PC = factor(names(eigenvals(rda_object)), names(eigenvals(rda_object))),
               percent = percent)
  ggplot(data, aes(x = PC, y = percent)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = mean(percent), color = "red", linetype = 2) +
    xlab("Principal component") +
    ylab("Variance explained (%)")+
    theme_bw()
  }

plot_pca <-
  function(sites, species, eig_vector, scaling){
  if(scaling == 1){
    plot <-
      ggplot()+
      # plot 
      stat_ellipse(data = sites, 
                   aes(x = PC1, y = PC2, color = Cruise, fill = Cruise), 
                   type = "norm", geom = "polygon",
                   size = 1.5,
                   alpha = 0.15,
                   level = 0.95) +
      # plot stations
      geom_point(data = sites, 
                 aes(x = PC1, y = PC2, color = Cruise)) +
      geom_text_repel(data = sites,
                      aes(x = PC1, y = PC2, color = Cruise, label = Station),
                      seed = 1) +
      # change axis label
      xlab(paste0("PC1 (", eig_vector[1], "% of total variance explained)")) +
      ylab(paste0("PC2 (", eig_vector[2], "% of total variance explained)")) +
      scale_color_manual(values = cruise_color) +
      scale_fill_manual(values = cruise_color) +
      coord_fixed() +
      theme_bw()
  }
  if(scaling == 2){
    plot <- 
      ggplot()+
      # plot 
      stat_ellipse(data = sites, 
                   aes(x = PC1, y = PC2, color = Cruise, fill = Cruise), 
                   type = "norm", geom = "polygon",
                   size = 1.5,
                   alpha = 0.15,
                   level = 0.95) +
      # plot stations
      geom_point(data = sites, 
                 aes(x = PC1, y = PC2, color = Cruise),
                 shape = 1) +
      
      # plot sp.
      geom_segment(data = species[species$Show == TRUE,],
                   aes(x = 0, y = 0, xend = PC1, yend = PC2),
                   size = .4, color = "black")+
      geom_label(data = species[species$Show == TRUE,],
                 aes(x = PC1, y = PC2, label = Taxon)) +
      
      # change axis label
      xlab(paste0("PC1 (", eig_vector[1], "% of total variance explained)")) +
      ylab(paste0("PC2 (", eig_vector[2], "% of total variance explained)")) +
      scale_color_manual(values = cruise_color) +
      scale_fill_manual(values = cruise_color) +
      coord_fixed() +
      theme_bw()
  }
  return(plot)
}
# produce pca plot for assemblage ordination
###############################################
# 1. Principal component analysis -- Count data
###############################################
# Plot PCA with count data 
count_dist <- log10_hell_dist(count_wide) 
count_pca <- rda(count_dist)

# get output
count_sc1 <- get_pca_output(count_pca, metadata = count_wide[1:2], scaling = 1)
count_sc2 <- get_pca_output(count_pca, metadata = count_wide[1:2], scaling = 2)
# plot eigenvalues
count_scree <- plot_scree(count_pca)
# plot pca plots
count_pca_sc1 <- plot_pca(count_sc1$pca_sites, count_sc1$pca_species, scaling = 1, eig_vector = count_sc1$pca_eig)
count_pca_sc2 <- plot_pca(count_sc2$pca_sites, count_sc2$pca_species, scaling = 2, eig_vector = count_sc2$pca_eig)

# output
ggsave(filename = "figure/polished/count_pca_sc1.png", plot = count_pca_sc1)
ggsave(filename = "figure/polished/count_pca_sc2.png", plot = count_pca_sc2)

###############################################
# 2. Principal component analysis -- Biomass data
###############################################
# Plot PCA with count data 
biomass_dist <- log10_hell_dist(biomass_wide) 
biomass_pca <- rda(biomass_dist)

# get output
biomass_sc1 <- get_pca_output(biomass_pca, metadata = biomass_wide[1:2], scaling = 1)
biomass_sc2 <- get_pca_output(biomass_pca, metadata = biomass_wide[1:2], scaling = 2)
# plot eigenvalues
count_scree <- plot_scree(count_pca)
# plot pca plots
biomass_pca_sc1 <- plot_pca(biomass_sc1$pca_sites, biomass_sc1$pca_species, scaling = 1, eig_vector = biomass_sc1$pca_eig)
biomass_pca_sc2 <- plot_pca(biomass_sc2$pca_sites, biomass_sc2$pca_species, scaling = 2, eig_vector = biomass_sc2$pca_eig)

# output
ggsave(filename = "figure/polished/biomass_pca_sc1.png", plot = biomass_pca_sc1)
ggsave(filename = "figure/polished/biomass_pca_sc2.png", plot = biomass_pca_sc2)
