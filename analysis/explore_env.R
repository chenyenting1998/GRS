####################################################
# Explore the Gaoping River-shelf environmental data
####################################################
# Description
# Data exploration of the GRS environmental condition 

# Author: Yen-Ting Chen
# Date of creation: 2023/07/08
# Date of last modification: 2023/07/08

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
library(Hmisc) # corrplot essentials
library(corrplot)
library(writexl) # xlsx
library(GRSmacrofauna)  # data package

# load variables
load("data/cruise_color.RData")
load("data/env_variables.RData")

################################
# 1. Convert env to long format
################################
# env long
env_long <- 
  env %>% 
  filter(!Station %in% c("GC1", "GS1")) %>% 
  pivot_longer(cols = all_of(env_variables),
               names_to = "Variables",
               values_to = "Values")

save(env_long, file = "data/env_long.RData")

###########################
# 2. Quantile-Quantile plot
###########################
env_qqplot <-
  ggplot(env_long, aes(sample = Values))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~Variables, scales = "free") +
  scale_color_manual(values = cruise_color) +
  theme_bw()

ggsave("figure/env_explore_qqplot.png", scale = 1.5, plot = env_qqplot)

#######################
# 3. Correlation matrix
####################### 
# make rcorr object
env_corr <-
  env[,env_variables] %>%
  as.matrix() %>% 
  rcorr() 

# plot corrplot
png(filename = "figure/env_explore_corr.png",
    width = 1600,
    height = 800)

corrplot(env_corr$r,
         # corrplot style
         method = "ellipse",
         type = "upper",
         diag = FALSE,
         addCoef.col = TRUE,
         # significance
         p.mat = env_corr$P,
         insig = "blank")
dev.off()

# create corr table
# locate insignificant cells
env_corr_insig <- env_corr$P <= 0.05
# remove insignificant values
env_corr$r[!env_corr_insig] <- NA
# remove diagonal values
diag(env_corr$r) <- NA

# output table
env_corr_table <- as.data.frame(env_corr$r)
write_xlsx(list(corr_table = env_corr_table,
                corr_table_rounded = round(env_corr_table, digits = 2)),
           "table/env_corr.xlsx")


############
# 4. Boxplot
############
env_boxplot <-
  env_long %>% 
  # filter(!Variables %in% c("Sigma-Theta", "WC")) %>% 
  ggplot(aes(x = Cruise, y = Values, color = Cruise))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = "jitter") +
  facet_wrap( ~ factor(Variables, env_variables), 
              scales = "free",
              labeller = as_labeller(env_variables_names, label_parsed)) +
  scale_color_manual(values = cruise_color) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8),
        # legend.position = c(0.91, 0.1),
        strip.text = element_text(size = 8))

ggsave("figure/env_explore_boxplot.png", 
       plot = env_boxplot,
       width = 8,
       height = 6,
       scale = 1.2)

##################################
# 5. Principle coordinate analysis
##################################
env_pca <- rda(scale(env[,env_variables]))

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
  geom_label(data = env_pca_sites, 
             aes(x = PC1, y = PC2, label = Station, color = Cruise)) +
  
  # change axis label
  xlab(paste0("PC1 (", env_pca_eig[1], "% of total variance explained)")) +
  ylab(paste0("PC2 (", env_pca_eig[2], "% of total variance explained)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed()+
  theme_bw() +
  theme(legend.position = c(0.1,0.1))

ggsave("figure/env_explore_pca_plot.png", 
       plot = env_pca_plot,
       width = 8,
       height = 6,
       scale = 1)

