##############################
# Setup - Taxa rank and color
##############################
# Description
# Assigning ranks and colors to each taxa for consistency in data visualization

# Author: Yen-Ting Chen
# Date of creation: Unknown
# Date of last modification: 2023/07/05

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
reduce_taxa_length <- function(data) {
  # label rare taxa as others
  data$Dominant <-
    if_else(data$percent < 2, "Others", "Dominant")
  # Relabel Dominant > 1% with their original name
  data[data$Dominant == "Dominant",]$Dominant <-
    data[data$Dominant == "Dominant",]$Taxa
  
  # factorize $Dominant_family
  rank_order <- data$Dominant[length(unique(data$Dominant)):1]
  data$Dominant <- factor(data$Dominant, rank_order)
  return(data)
}

###################
# 1. Data cleaning
###################
data <- 
  macrofauna_biomass %>% 
  # cruise to month
  filter(Taxon == "Polychaeta") %>% 
  # remove NAs and unknowns
  filter(!is.na(Family)) %>% 
  filter(!Family %in% c("Unknown", "unknown")) %>% 
  # filter head intact specimens
  filter(Condition %in% c("C", "FH")) %>% 
  # use one column to describe taxonomic resolution
  mutate(Taxa = if_else(!is.na(Genus), paste0("G. ", Genus), paste0("F. ", Family)))

# Note: there is a possibility that a fragmented specimen belongs to 
# a family/genus that does not appear in complete/head-intact specimen.

#################################
# 2. Define family rank by density
#################################
# Rank taxa
rank_den <-
  data %>%
  group_by(Taxa) %>%
  summarize(Count = n()) %>%
  mutate(percent = round(Count / sum(Count) * 100, 2)) %>%
  arrange(desc(percent))

rank_den <- reduce_taxa_length(rank_den)

#################################
# 3. Define taxa rank by biomass
#################################
# rank_bio <-
#   data %>%
#   group_by(Taxa) %>%
#   summarize(WM = sum(WM)) %>%
#   mutate(percent = round(WM / sum(WM) * 100, 2)) %>%
#   arrange(desc(percent))
# 
# rank_bio <- reduce_taxa_length(rank_bio)

##################################
# 4. Assign color by reduced taxa
##################################
# all taxa list
all_taxa <- as.character(rank_den$Dominant) %>% unique()

# how many colors do we need
all_taxa

# grab colors
taxa_color <- kelly.colors(length(all_taxa) + 1)[-1]
# assign colors to taxa
names(taxa_color) <- all_taxa

####################
# 5. den_family_color
####################
# assign fixed color 
taxa_den_color <- taxa_color[names(taxa_color) %in% rank_den$Dominant] 
# reorder color with density rank
taxa_den_color <- taxa_den_color[match(levels(rank_den$Dominant), names(taxa_den_color))]

####################
# 6. bio_family_color
####################
# assign fixed color 
# family_bio_color <- family_color[names(family_color) %in% rank_bio$Dominant_family]
# # reorder color with biomass rank
# family_bio_color <- family_bio_color[match(levels(rank_bio$Dominant_family), names(family_bio_color))]

############
# 5. Output
############
save(rank_den, #rank_bio,
     file = "data/polychaete_rank.Rdata")
save(taxa_den_color, #family_bio_color,
     file = "data/polychaete_family_color.Rdata")
