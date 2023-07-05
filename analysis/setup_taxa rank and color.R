###############################
# Setup - Taxa rank and color
##############################
# Description
# Assigning ranks and colors to each taxa for consistency in data visualization

# Author: Yen-Ting Chen
# Date of creation: Unknown
# Date of last modification: 2023/07/05

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
  rank_den$Taxa <-
    if_else(rank_den$percent < 1, "Others", "Dominant")
  # Relabel taxa > 1% with their original name
  rank_den[rank_den$Taxa == "Dominant",]$Taxa <-
    rank_den[rank_den$Taxa == "Dominant",]$Taxon
  
  # factorize $Taxa
  rank_den_order <- rank_den$Taxa[length(unique(rank_den$Taxa)):1]
  rank_den$Taxa <- factor(rank_den$Taxa, rank_den_order)
  return(rank_den)
}

###################
# 1. Data cleaning
###################
data <-
  macrofauna_biomass %>%
  # include only head-intact individuals
  filter(Condition %in% c("C", "FH")) %>%
  # remove unclear/ pelagic individuals
  filter(!Taxon %in% c("Unknown", "Calanoida")) %>%
  # remove hydroids noted as stalks since no tissue are attached to the specimens
  filter(!(Taxon %in% "Hydrozoa" & Note %in% "Stalk")) %>%
  # remove the stony coral; only one individual and biases the community biomass
  # maybe I should not include the core at all?
  filter(!Taxon %in% c("Scleractinia"))

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

rank_bio <- reduce_taxa_length(rank_bio)

##################################
# 4. Assign color by reduced taxa
##################################
# all taxa list
all_taxa <- unique(c(rank_den$Taxa, rank_bio$Taxa))
# grab colors
taxa_color <- kelly.colors(length(all_taxa) + 1)[-1]
# assign colors to taxa
names(taxa_color) <- all_taxa

# den_taxa_color
taxa_den_color <- taxa_color[names(taxa_color) %in% rank_den$Taxa]

# bio_taxa_color
taxa_bio_color <- taxa_color[names(taxa_color) %in% rank_bio$Taxa]

############
# 5. Output
############
save(rank_den, rank_bio,
     file = "data/taxa_rank.Rdata")
save(taxa_den_color, taxa_bio_color,
     file = "data/taxa_color.Rdata")