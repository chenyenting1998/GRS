##################################
# Gaoping River-shelf environment
##################################
# Description
# Analyze the environmental condition of the GRS

# Author: Yen-Ting Chen
# Date of creation: Unknown
# Date of last modification: 2023/07/06

#######################
# Set up environment
#######################
# clean environment
rm(list = ls())

# Load packages
library(dplyr) # data manipulation
library(tidyr)
library(vegan) # data analysis
library(rstatix) # pipe-friendly stats package
library(ggplot2) # visualization
library(GGally)
library(ggrepel)
library(writexl) # xlsx
library(GRSmacrofauna)  # data package

# load variables
load("data/cruise_color.RData")
load("data/env_variables.RData")
load("data/env_long.RData")

################################
# 1. Select interested variables
################################
# omit variables
omit <- c("Salinity", "Density", "SigmaTheta", "Oxygen", "Transmission",
          "Clay", "Silt", "Sand",
          "delta13C",
          "TN",
          "WC")
env_spatialvariables <- c("Depth", "DRM")
# Reasoning:
#  Salinity: salinity were within the open water range, and snapshots 
#            of salinity change cannot inform us the variability of salinity.
#  Density: not known to be an important factor.
#  Transmission: strong correlation with temperature
#  Oxygen: bottom water is well oxygenated.
#  clay, silt, sand: D50 was used as the sole granulometry proxy 
#  Delta13C: strong correlation with Chla; mixed signals based on geological research
#  TN: strong collinearity with porosity.
#  Water content: strong collinearity with porosity.


# subset env
env_variables_selected <- colnames(env)[!colnames(env) %in% c(env_metadata, env_spatialvariables, omit)]
env_selected <- env[,c(env_variables_selected)]

#############
# 2. Pairplot
#############
env_pairplot <-
  env_selected %>% 
  ggpairs()

ggsave("figure/polished/env_pairplot.png", scale = 1.4, plot = env_pairplot)

###########################
# 3. Wilcoxon rank-sum test 
###########################
# set up an empty object
wilcox_table <- 
  env_long %>% 
  group_by(Variables) %>% 
  wilcox_test(Values ~ Cruise,
              comparisons = list(c("OR1-1219", "OR1-1242")))
# output
write_xlsx(wilcox_table, "table/env_wilcox_table.xlsx")

###########################
# 4. PERMANOVA and PERMDISP
###########################
set.seed(10)
# run permanova
env_permanova <- 
  adonis2(scale(env_selected) ~ Cruise,
  # adonis2(scale(env_selected) ~ Cruise,
          data = env, 
          method = "euclidean",
          permutations = 99999)
# run permdisp
env_disp <-
  betadisper(vegdist(env_selected,
                     method = "euclidean", scale = TRUE), 
             group = env$Cruise)
# env_disp <- betadisper(vegdist(env_selected, method = "euclidean", scale = TRUE), group = env$Cruise)
env_permdisp <- permutest(env_disp, permutations = 99999)

# output
permanova_output <-
  list(PERMANOVA = as.data.frame(env_permanova),
       PERMDISP = as.data.frame(env_permdisp$tab))
write_xlsx(permanova_output,
           "table/env_permanova.xlsx")

##################################
# 5. Principle coordinate analysis
##################################
env_pca <- rda(scale(env_selected))

# extract sites
env_pca_sites <- 
  scores(env_pca, scaling = 1)$sites %>% 
  as.data.frame() %>% 
  cbind(env[c("Cruise", "Station")])

# extract species (env var.)
env_pca_species <- 
  scores(env_pca, scaling = 1)$species %>%
  as.data.frame()

# change names to abbr
env_pca_species$abbr <- env_variables_abbr[match(rownames(env_pca_species), names(env_variables_abbr))]

# extract eig
env_pca_eig <- 
  round(eigenvals(env_pca) / sum(eigenvals(env_pca)) * 100, 2) %>% 
  as.vector()

env_pca_plot <-
  ggplot() +
  # plot variable segments
  geom_segment(data = env_pca_species,
               aes(x = PC1, y = PC2, xend = 0, yend = 0)) +
  # plot variable text
  geom_text(data = env_pca_species,
            aes(x = PC1 * 1.09, 
                y = PC2 * 1.09, 
                label = abbr),
            parse = TRUE) +
  # plot stations
  geom_point(data = env_pca_sites, 
             aes(x = PC1, y = PC2, color = Cruise)) +
  geom_label(data = env_pca_sites,#[-c(2:3),], 
             aes(x = PC1, y = PC2, label = Station, color = Cruise)) +
  # geom_label_repel(data = env_all_pca_sites[c(2:3),],
  #                  aes(x = PC1, y = PC2, label = Station, color = Cruise),
  #                  min.segment.length = 0) +
  # change axis label
  xlab(paste0("PC1 (", env_pca_eig[1], "% of total variance explained)")) +
  ylab(paste0("PC2 (", env_pca_eig[2], "% of total variance explained)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed()+
  theme_bw() +
  theme(legend.position = c(0.1,0.1))

ggsave("figure/polished/env_pca_plot.png", 
       plot = env_pca_plot,
       width = 8,
       height = 6,
       scale = 1)

###########
# 6. Output
###########
# add metadata to env_selected
env_selected <- env[,c(env_metadata, env_variables_selected)]
save(env_variables_selected, 
     env_selected, 
     file = "data/env_selected.RData")
