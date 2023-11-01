##############################
# Setup - Taxa rank and color
##############################
# Description
# Assigning ranks and colors to each taxa for consistency in data visualization

# Author: Yen-Ting Chen
# Date of creation: Unknown
# Date of last modification: 2023/10/24

#####################
# Set up environment
#####################
rm(list = ls())

################
# Load packages
################
library(dplyr)
library(Polychrome)
library(GRSmacrofauna)

##################
# Define function 
##################
# could only be used in this script
reduce_taxa_length <- function(data, percent_threshold = 1) {
  # label rare taxa as others
  data$Taxa <-
    if_else(data$percent < percent_threshold, "Others", "Dominant")
  # Relabel taxa > 1% with their original name
  data[data$Taxa == "Dominant",]$Taxa <-
    data[data$Taxa == "Dominant",]$Taxon
  
  # factorize $Taxa
  rank_order <- data$Taxa[length(unique(data$Taxa)):1]
  data$Taxa <- factor(data$Taxa, rank_order)
  return(data)
}

###################
# 1. Data cleaning
###################
data <-
  macrofauna_biomass %>%
  # exclude GC1 and GS1
  filter(!Station %in% c("GC1", "GS1")) %>% 
  # include only head-intact individuals
  filter(Condition %in% c("C", "FH")) %>%
  # remove unknown/terrestrial/pelagic individuals
  filter(!Taxon %in% c("Unknown", "Insecta", "Calanoida")) %>%
  # remove hydroids noted as stalks since no tissue are attached to the specimens
  filter(!(Taxon %in% "Hydrozoa" & Note %in% "Stalk"))# %>%
  # remove the stony coral; only one individual and biases the community biomass
  # maybe I should not include the core at all?
  # filter(!Taxon %in% c("Scleractinia"))

#################################
# 2. Define taxa rank by density
#################################
# Rank taxa
rank_den <-
  data %>%
  group_by(Taxon) %>%
  summarize(Count = n()) %>%
  mutate(percent = round(Count / sum(Count) * 100, 2)) %>%
  arrange(desc(percent))

rank_den <- reduce_taxa_length(rank_den)

#################################
# 3. Define taxa rank by biomass
#################################
rank_bio <-
  data %>%
  group_by(Taxon) %>%
  summarize(Biomass = sum(Size)) %>%
  mutate(percent = round(Biomass / sum(Biomass) * 100, 2)) %>%
  arrange(desc(percent))

rank_bio <- reduce_taxa_length(rank_bio, .5)
View(rank_bio)
##################################
# 4. Assign color by reduced taxa
##################################
# all taxa list
all_taxa <-
  as.character(rank_bio$Taxa) %>%
  c(as.character(rank_den$Taxa)) %>%
  unique()
# grab colors
taxa_color <- kelly.colors(length(all_taxa) + 1)[-1]
# assign colors to taxa
names(taxa_color) <- all_taxa

####################
# 5. den_taxa_color
####################
# assign fixed color 
taxa_den_color <- taxa_color[names(taxa_color) %in% rank_den$Taxa] 
# reorder color with density rank
taxa_den_color <- taxa_den_color[match(levels(rank_den$Taxa), names(taxa_den_color))]

####################
# 6. bio_taxa_color
####################
# assign fixed color 
taxa_bio_color <- taxa_color[names(taxa_color) %in% rank_bio$Taxa]
# reorder color with biomass rank
taxa_bio_color <- taxa_bio_color[match(levels(rank_bio$Taxa), names(taxa_bio_color))]

############
# 5. Output
############
save(rank_den, rank_bio,
     file = "data/taxa_rank.Rdata")
save(taxa_den_color, taxa_bio_color,
     file = "data/taxa_color.Rdata")
