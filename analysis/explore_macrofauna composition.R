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

# Define multiple corer area; diameter = 10 cm
mc_area <- (0.1 / 2) ^ 2 * pi

# Load functions
source("analysis/box.cox.chord.R")

# Define internal function
plot_taxa_boxplot <- function(data, 
                              value = "Value", 
                              exp = 1,
                              filter_zero = TRUE){
  # remove metadata columns
  data_preped <- 
    data %>% 
    ungroup() %>% 
    select(- Cruise, - Station, - Deployment, - Tube)
  # take exp
  data_preped <- data_preped^exp
  # convert to long data
  data_preped <- 
    data_preped %>%
    mutate(Cruise = data$Cruise) %>% 
    pivot_longer(cols = -Cruise,
                 names_to = "Taxon",
                 values_to = "Value")
  # remove zeros
  if(isTRUE(filter_zero)){
    data_preped <- 
      data_preped %>% 
      filter(Value != 0)
  }
  #
  figure <-
    data_preped %>% 
    ggplot() +
    geom_boxplot(aes(x = Taxon, y = Value), outlier.shape = NA) +
    geom_point(aes(x = Taxon,   y = Value, color = Cruise), position = "jitter") +
    scale_color_manual(values = cruise_color) +
    ylab(paste0("log[10] (", value, ")")) +
    coord_flip() +
    theme_bw()
  return(figure)
}

plot_taxa_dist <- function(data, 
                           value = "Value", 
                           exp = 1,
                           filter_zero = TRUE){
  # remove metadata columns
  data_preped <- 
    data %>% 
    ungroup() %>% 
    select(- Cruise, - Station, - Deployment, - Tube) 
  
  # take the exp
  data_preped <- data_preped^exp
  
  # convert to long data
  data_preped <-
    data_preped %>% 
    mutate(Cruise = data$Cruise) %>% 
    pivot_longer(cols = -Cruise,
                 names_to = "Taxon",
                 values_to = "Value")
  # remove zeros
  if(isTRUE(filter_zero)){
  data_preped <- 
    data_preped %>% 
    filter(Value != 0)
  }
  figure <-
    data_preped %>% 
    ggplot() +
    geom_density(aes(y = Value)) +
    scale_color_manual(values = cruise_color) +
    ylab(paste0("log[10] (", value, ")")) +
    facet_wrap(~Taxon, scales = "free") +
    coord_flip() +
    theme_bw()
  return(figure)
}

#############################################
# 1. Calculate unit area density and biomass 
#############################################
composition <- 
  macrofauna_biomass %>%
  # exclude GC1 and GS1
  filter(!Station %in% c("GC1", "GS1")) %>% 
  
  # include only head-intact individuals
  filter(Condition %in% c("C", "FH")) %>%
  
  # remove unknown/terrestrial/pelagic individuals
  filter(!Taxon %in% c("Unknown", "Insecta","Calanoida")) %>%
  
  # remove hydroids noted as stalks since no tissue are attached to the specimens
  filter(!(Taxon %in% "Hydrozoa" & Note %in% "Stalk")) %>%
  
  # remove the stony coral; only one individual and biases the community biomass
  # maybe I should not include the core at all?
  filter(!Taxon %in% c("Scleractinia")) %>% 
  group_by(Cruise, Station, Deployment, Tube, Taxon) %>% 
  summarize(Count = n(),
            Biomass = sum(WM))

##############################
# 2. Plot density composition
##############################
density_composition <-
  add_coarse_taxa(composition, match_file = rank_den, output = "Taxa") %>%
  ggplot(aes(
    x = Station,
    y = Count / mc_area / 3, # unit area; 3 cores
    fill = Taxa
  )) +
  geom_bar(stat = "identity") +
  facet_grid(~ Cruise, scales = "free") +
  scale_fill_manual(values = taxa_den_color) +
  ylab(Density~(ind.~m^{-2})) +
  theme_bw()

density_percentage_composition <-
  add_coarse_taxa(composition, match_file = rank_den, output = "Taxa") %>%
  ggplot(aes(x = Station, y = Count, fill = Taxa)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid( ~ Cruise, scales = "free") +
  scale_fill_manual(values = taxa_den_color) +
  scale_y_continuous(labels = scales::percent) +
  ylab(Density~("%")) +
  theme_bw()

##############################
# 3. Plot biomass composition
##############################
biomass_composition <-
  add_coarse_taxa(composition, match_file = rank_bio, output = "Taxa") %>%
  ggplot(aes(
    x = Station,
    y = Biomass / mc_area / 1000 / 3, # unit area; mg->g; 3 cores
    fill = Taxa
  )) +
  geom_bar(stat = "identity") +
  facet_grid( ~ Cruise, scales = "free") +
  scale_fill_manual(values = taxa_bio_color) +
  ylab(Wet ~ mass ~ (g~m ^ {-2})) +
  theme_bw()

biomass_percentage_composition <- 
  add_coarse_taxa(composition, match_file = rank_bio, output = "Taxa") %>% 
  ggplot(aes(x = Station, y = Biomass, fill = Taxa)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~Cruise, scales = "free") +
  scale_fill_manual(values = taxa_bio_color)+
  scale_y_continuous(labels = scales::percent) +
  ylab(Wet~mass~("%")) +
  theme_bw()

######################
# 4. Render wide data
######################
# calculate density to wide data.frame
count_wide <-
  composition %>% 
  select(-Biomass) %>% 
  pivot_wider(names_from = Taxon, values_from = Count, values_fill = 0) %>% 
  ungroup()

# calculate biomass to wide data.frame
biomass_wide <-
  composition %>% 
  select(-Count) %>% 
  pivot_wider(names_from = Taxon, values_from = Biomass, values_fill = 0) %>% 
  ungroup()

# yield Box-Cox-chord transformation results
# Note that n < 3*p, Dagnelie's Test too liberal 
# if p > 0.05, the Dagnelie's test result is trustworthy
count_BCD <- as.data.frame(BCD(count_wide[-(1:4)]))
biomass_BCD <- as.data.frame(BCD(biomass_wide[-(1:4)]))
write_xlsx(count_BCD, path = "table/count_BCD.xlsx")
write_xlsx(biomass_BCD, path = "table/biomass_BCD.xlsx")

##########
# 5. Plot
##########
# wide data transformed following Anderson et al.(2006)
# plot the frequency of # of individuals by each taxa (i.e., matrix entries)
hist(as.vector(as.matrix(count_wide[-c(1:4)])),
     main = "Histogram of taxa count",
     xlab = "Number of individuals")

# species-specific density plot
count_taxa_distplot <- plot_taxa_dist(count_wide, "Count", exp = 0.5)
biomass_taxa_distplot <- plot_taxa_dist(biomass_wide, "Biomass", exp = 0.5)

# species_specific boxplot
count_taxa_boxplot <- plot_taxa_boxplot(count_wide, "Count", exp = 0.5)
biomass_taxa_boxplot <- plot_taxa_boxplot(biomass_wide, "Biomass", exp = 0.5)

##############
# Data output
##############
save(composition, file = "data/macrofauna composition.RData")
save(count_wide, biomass_wide, file = "data/wide_data.Rdata")

##################
# Figure output
##################
ggsave("figure/density_composition.png", plot = density_composition)
ggsave("figure/biomass_composition.png", plot = biomass_composition)
ggsave("figure/density_percentage_composition.png", plot = density_percentage_composition)
ggsave("figure/biomass_percentage_composition.png", plot = biomass_percentage_composition)
ggsave("figure/count_taxa_boxplot.png", plot = count_taxa_boxplot)
ggsave("figure/biomass_taxa_boxplot.png", plot = biomass_taxa_boxplot)
ggsave("figure/count_taxa_distplot.png", plot = count_taxa_distplot)
ggsave("figure/biomass_taxa_distplot.png", plot = biomass_taxa_distplot)
