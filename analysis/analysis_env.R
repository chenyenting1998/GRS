##################################
# Gaoping River-shelf environment
##################################
# Description
# Analyze the environmental condition of the GRS

# Author: Yen-Ting Chen
# Date of creation: Unknown
# Date of last modification: 2023/10/20

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

# add month
add_month <- function(x){
  x$Month <- if_else(x$Cruise == "OR1-1219", "March", "October")
  return(x)
}

################################
# 1. Select interested variables
################################
# omit variables
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

omit <- c("Salinity", "Density", "SigmaTheta", "Oxygen", "Transmission",
          "Clay", "Silt", "Sand",
          "delta13C",
          "TN",
          "WC")

# subset env
env_variables_selected <- colnames(env)[!colnames(env) %in% c(env_metadata, env_variables_spatial, omit)]
env_selected <- env[,c(env_metadata, env_variables_selected)]

#############
# 2. Pairplot
#############
env_pairplot <-
  env %>% 
  add_month() %>% 
  ggpairs(columns = c(env_variables_spatial, env_variables_selected),
          aes(color = Month)) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  theme_bw()

ggsave("figure/pairplot/env_pairplot.png", 
       scale = 1.4,
       h = 5,
       w = 8,
       plot = env_pairplot)

###########################
# 3. Wilcoxon rank-sum test 
###########################
# set up an empty object
wilcox_table <- 
  env %>% 
  pivot_longer(cols = all_of(env_variables), 
               names_to = "Variables", 
               values_to = "Values") %>% 
  add_month() %>% 
  group_by(Variables) %>% 
  wilcox_test(Values ~ Month,
              comparisons = list(c("March", "October")))

# subset
wilcox_table_selected <- wilcox_table[wilcox_table$Variables %in% env_variables_selected,]

# output
write_xlsx(list(wilcox_table = wilcox_table,
                wilcox_table_selected = wilcox_table_selected), 
           "table/wilcox/env_wilcox_table.xlsx")

###########################
# 4. PERMANOVA and PERMDISP
###########################
set.seed(10)
# run permanova
env_sel_month <- add_month(env_selected)
env_permanova <- 
  adonis2(scale(env_sel_month[env_variables_selected]) ~ Month,
  # adonis2(scale(env_sel_month) ~ Month,
          data = env_sel_month, 
          method = "euclidean",
          permutations = 9999)
# run permdisp
env_disp <-
  betadisper(vegdist(env_sel_month[env_variables_selected],
                     method = "euclidean", scale = TRUE), 
             group = env_sel_month$Month)
# env_disp <- betadisper(vegdist(env_selected, method = "euclidean", scale = TRUE), group = env$Cruise)
env_permdisp <- permutest(env_disp, permutations = 9999)

# output
permanova_output <-
  list(PERMANOVA = as.data.frame(env_permanova),
       PERMDISP = as.data.frame(env_permdisp$tab))
write_xlsx(permanova_output,
           "table/permanova/env_permanova.xlsx")

##################################
# 5. Principle coordinate analysis
##################################
env_pca <- rda(scale(env_sel_month[env_variables_selected]))

# extract sites
env_pca_sites <- 
  scores(env_pca, scaling = 1)$sites %>% 
  as.data.frame() %>% 
  cbind(env_sel_month[c("Month", "Station")])

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
  # # plot ellipse
  # geom_polygon(data = env_pca_sites, 
  #              aes(x = PC1, y = PC2, color = Month, fill = Month),
  #              stat = "ellipse",
  #              size = 1.5,
  #              alpha = 0.3) +
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
             aes(x = PC1, y = PC2, color = Month)) +
  geom_label(data = env_pca_sites,#[-c(2:3),], 
             aes(x = PC1, y = PC2, label = Station, color = Month)) +
  # geom_label_repel(data = env_all_pca_sites[c(2:3),],
  #                  aes(x = PC1, y = PC2, label = Station, color = Month),
  #                  min.segment.length = 0) +
  # change axis label
  xlab(paste0("PC1 (", env_pca_eig[1], "% of the total variance)")) +
  ylab(paste0("PC2 (", env_pca_eig[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed()+
  theme_bw() +
  theme(legend.position = c(0.85,
                            0.85))

ggsave("figure/pca/env_pca_plot.png", 
       plot = env_pca_plot,
       width = 8,
       height = 6,
       scale = 1)

###########
# 6. Output
###########
# add metadata to env_selected
save(env_variables_selected, 
     env_selected, 
     file = "data/env_selected.RData")

save(env_pca, file = "data/env_pca.RData")
