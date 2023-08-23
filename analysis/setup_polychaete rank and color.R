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
  data$Dominant_family <-
    if_else(data$percent < 1, "Others", "Dominant")
  # Relabel Dominant_family > 1% with their original name
  data[data$Dominant_family == "Dominant",]$Dominant_family <-
    data[data$Dominant_family == "Dominant",]$Family
  
  # factorize $Dominant_family
  rank_order <- data$Dominant_family[length(unique(data$Dominant_family)):1]
  data$Dominant_family <- factor(data$Dominant_family, rank_order)
  return(data)
}

###################
# 1. Data cleaning
###################
data <- 
  macrofauna_biomass %>% 
  filter(Taxon == "Polychaeta") %>% 
  # remove families entries
  filter(!is.na(Family)) %>% 
  filter(!Family %in% c("Unknown", "unknown")) %>% 
  # include specimen with good condition
  filter(Condition %in% c("C", "FH"))

# Note: there is a possibility that a fragmented specimen belongs to 
# a family/genus that does not appear in complete/head-intact specimen.

#################################
# 2. Define family rank by density
#################################
# Rank taxa
rank_den <-
  data %>%
  group_by(Family) %>%
  summarize(Count = n()) %>%
  mutate(percent = round(Count / sum(Count) * 100, 2)) %>%
  arrange(desc(percent))

rank_den <- reduce_taxa_length(rank_den)

#################################
# 3. Define taxa rank by biomass
#################################
rank_bio <-
  data %>%
  group_by(Family) %>%
  summarize(WM = sum(WM)) %>%
  mutate(percent = round(WM / sum(WM) * 100, 2)) %>%
  arrange(desc(percent))

rank_bio <- reduce_taxa_length(rank_bio)

##################################
# 4. Assign color by reduced taxa
##################################
# all taxa list
all_family <-
  as.character(rank_bio$Dominant_family) %>%
  c(as.character(rank_den$Dominant_family)) %>%
  unique()

# grab colors
family_color <- kelly.colors(length(all_family) + 1)[-1]
# assign colors to taxa
names(family_color) <- all_family

####################
# 5. den_family_color
####################
# assign fixed color 
family_den_color <- family_color[names(family_color) %in% rank_den$Dominant_family] 
# reorder color with density rank
family_den_color <- family_den_color[match(levels(rank_den$Dominant_family), names(family_den_color))]

####################
# 6. bio_family_color
####################
# assign fixed color 
family_bio_color <- family_color[names(family_color) %in% rank_bio$Dominant_family]
# reorder color with biomass rank
family_bio_color <- family_bio_color[match(levels(rank_bio$Dominant_family), names(family_bio_color))]

############
# 5. Output
############
save(rank_den, rank_bio,
     file = "data/polychaete_rank.Rdata")
save(family_den_color, family_bio_color,
     file = "data/polychaete_family_color.Rdata")
