#########################
# Macrofauna composition
#########################
# Description
# Assigning ranks and colors to each taxa for consistency in data visualization

# Author: Yen-Ting Chen
# Date of creation: 2023/07/05
# Date of last modification: 2023/07/05

#######################
# Clean up environment
#######################
rm(list = ls())

################
# Load packages
################
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(GGally)
library(GRSmacrofauna)
library(writexl)

############
# Load data
############
load("data/taxa_rank.Rdata")
load("data/taxa_color.Rdata")
load("data/cruise_color.Rdata")

################
# Define values
################
# multiple corer area; diameter = 10 cm
mc_area <- (0.1/2)^2*pi

###################
# Define functions
###################
plot_taxa_boxplot <- function(data, value = "Value"){
  figure <- 
    data %>% 
    ungroup() %>% 
    select(- Cruise, - Station, - Deployment, - Tube) %>% 
    decostand(method = "log", logbase = 10) %>% 
    mutate(Cruise = data$Cruise) %>% 
    pivot_longer(cols = -Cruise,
                 names_to = "Taxon",
                 values_to = "Value") %>%
    filter(Value != 0) %>%
    ggplot() +
    geom_boxplot(aes(x = Taxon, y = Value), outlier.shape = NA) +
    geom_point(aes(x = Taxon,   y = Value, color = Cruise), position = "jitter") +
    scale_color_manual(values = cruise_color) +
    ylab(paste0("log[10] (", value, ")")) +
    coord_flip() +
    theme_bw()
  return(figure)
}

#############################################
# 1. Calculate unit area biomass and density 
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
  summarize(Abundance = n(),
            Biomass = sum(WM))

##############################
# 2. Plot density composition
##############################
density_composition <-
  add_coarse_taxa(composition, match_file = rank_den, output = "Taxa") %>%
  ggplot(aes(
    x = Station,
    y = Abundance / mc_area / 3, # unit area; 3 cores
    fill = Taxa
  )) +
  geom_bar(stat = "identity") +
  facet_grid(~ Cruise, scales = "free") +
  scale_fill_manual(values = taxa_den_color) +
  ylab(Density~(ind.~m^{-2})) +
  theme_bw()

density_percentage_composition <-
  add_coarse_taxa(composition, match_file = rank_den, output = "Taxa") %>%
  ggplot(aes(x = Station, y = Abundance, fill = Taxa)) +
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
  ylab(Wet ~ weight ~ (g~m ^ {-2})) +
  theme_bw()

biomass_percentage_composition <- 
  add_coarse_taxa(composition, match_file = rank_bio, output = "Taxa") %>% 
  ggplot(aes(x = Station, y = Biomass, fill = Taxa)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~Cruise, scales = "free") +
  scale_fill_manual(values = taxa_bio_color)+
  scale_y_continuous(labels = scales::percent) +
  ylab(Wet~weight~("%")) +
  theme_bw()

######################
# 4. Render wide data
######################
# calculate density to wide data.frame
density_wide <-
  composition %>% 
  select(-Biomass) %>% 
  pivot_wider(names_from = Taxon, values_from = Abundance, values_fill = 0)

# calculate biomass to wide data.frame
biomass_wide <-
  composition %>% 
  select(-Abundance) %>% 
  pivot_wider(names_from = Taxon, values_from = Biomass, values_fill = 0)

####################
# 5. 
#######################
# wide data transformed following Anderson et al.(2006)
# plot the frequency of # of individuals by each taxa (i.e., matrix entries)
hist(as.vector(as.matrix(density_wide[-c(1:4)])),
     main = "Histogram of taxa count",
     xlab = "Number of individuals")

# species_specific boxplot
density_taxa_boxplot <- plot_taxa_boxplot(density_wide, "Count")
biomass_taxa_boxplot <- plot_taxa_boxplot(biomass_wide, "Biomass")


#########
# Data output
#########
save(composition, file = "data/macrofauna composition.RData")
save(density_wide, biomass_wide, file = "data/wide_data.Rdata")

##################
# Figure output
##################
ggsave("figure/density_composition.png", plot = density_composition)
ggsave("figure/biomass_composition.png", plot = biomass_composition)
ggsave("figure/density_percentage_composition.png", plot = density_percentage_composition)
ggsave("figure/biomass_percentage_composition.png", plot = biomass_percentage_composition)
ggsave("figure/density_taxa_boxplot.png", plot = density_taxa_boxplot)
ggsave("figure/biomass_taxa_boxplot.png", plot = biomass_taxa_boxplot)
