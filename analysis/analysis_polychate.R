########################
# Polychaete composition
########################
# Description
# Polychate assemblage data exploration

# Author: Yen-Ting Chen
# Date of creation: 2023/08/15
# Date of last modification: 2023/08/15

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
library(ggrepel)
library(Polychrome)
library(GRSmacrofauna)
library(writexl)

# Load data
load("data/polychaete_rank.Rdata")
load("data/polychaete_family_color.Rdata")

# Load functions
source("source/box.cox.chord.R")

pol <- 
  macrofauna_biomass %>% 
  filter(Taxon == "Polychaeta")

#####################
# Inspect condition #
#####################
ggplot(pol, aes(x = Family, fill = Condition)) +
  geom_bar(position = "stack") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0.99))

ggplot(pol, aes(x = Family, y = WM, fill = Condition)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0.99))

# update pol
pol <- 
  macrofauna_biomass %>% 
  filter(Taxon == "Polychaeta") %>% 
  # remove families entries
  filter(!is.na(Family)) %>% 
  filter(!Family %in% c("Unknown", "unknown")) %>% 
  # include specimen with good condition
  filter(Condition %in% c("C", "FH"))

# Note: there is a possibility that a fragmented specimen belongs to 
# a family/genus that does not appear in complete/head-intact specimen.

############
# Barplots #
############
# family-level long dataframe
pol_family_long <- 
  pol %>% 
  group_by(Cruise, Station, Deployment, Tube, Family) %>% 
  summarise(Count = n(),
            WM = sum(WM)) %>% 
  ungroup()

# plotting
add_coarse_family <- function(data, 
                              match_file, 
                              output = "Family"){
  # collecting the two columns needed for this function
  # Family and the columne named "output"
  extract <- 
    data.frame(Family = match_file$Family, 
               match_file[,colnames(match_file) == output])
  # matching rows of two dataframes
  comparing <- match(data$Family, extract$Family)
  # 
  result <- extract[comparing, 2]
  
  # attach the coase family as a new column
  data$output <- result
  # change the new column's name
  colnames(data)[ncol(data)] <- output
  return(data)
}

# set up ggplot function
plot_comp <- function(data,
                      rank,
                      xx,
                      yy, 
                      position,
                      color){
  data_rank <-
    data %>% 
    add_coarse_family(rank, output = "Dominant_family")
    
  plot <- 
    ggplot(data_rank, 
           aes(x = .data[[xx]], y = .data[[yy]], fill = Dominant_family)) +
    geom_bar(stat = "identity", position = position) +
    facet_grid(~Cruise, scales = "free") +
    scale_fill_manual(values = color) +
    theme_bw()
  return(plot)
}

# den
plot_comp(pol_family_long, rank_den, "Station", "Count", "stack", family_den_color)
plot_comp(pol_family_long, rank_den, "Station", "Count", "fill", family_den_color)
# bio
plot_comp(pol_family_long, rank_bio, "Station", "WM", "stack", family_bio_color)
plot_comp(pol_family_long, rank_bio, "Station", "WM", "fill", family_bio_color)

###########
# Heatmap #
###########
pol_count_wide <-
  pol_family_long %>% 
  select(-WM) %>% 
  pivot_wider(names_from = Family, 
              values_from = Count,
              values_fill = 0)

pol_count_long <-
  pol_count_wide %>% 
  pivot_longer(cols = -(1:4), 
               names_to = "Family", 
               values_to = "Count")

# heatmap of all station
pol_count_long %>% 
  ggplot() +
  geom_tile(aes(x = Station, 
                y = factor(Family, levels = rank_den$Family), 
                fill = log10(Count + 1))) +
  facet_grid(~Cruise, scales = "free")+
  viridis::scale_fill_viridis() +
  ylab("Family")+
  theme_bw()

# heatmap of duplicated stations
pol_count_long %>%
  filter(Station %in% c("GC1", "GS1",
                        "S3", "S5", "S6", "S7")) %>% 
  ggplot() +
  geom_tile(aes(x = Cruise, 
                y = factor(Family, levels = rank_den$Family), 
                fill = log10(Count + 1))) +
  facet_grid(~Station, scales = "free")+
  viridis::scale_fill_viridis() +
  ylab('Family') +
  theme_bw()

###########
# iNEXT #
###########
library(iNEXT)
setColnames <- function(data, name){
  colnames(data)
}
pol_wide_station <- 
  pol_count_long %>%
  group_by(Cruise, Station, Family) %>% 
  summarise(Count = sum(Count)) %>% 
  pivot_wider(names_from = "Family", 
              values_from = "Count",
              values_fill = 0) %>% 
  ungroup()

pol_wide_station_t <-
  pol_wide_station %>%
  select(-Cruise, - Station) %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame()
# set colnames
colnames(pol_wide_station_t) <- paste0(pol_wide_station$Cruise,
                                       "_",
                                       pol_wide_station$Station)
# run iNEXT  
x <- iNEXT(pol_wide_station_t)

# output
x$DataInfo
x$AsyEst
ggiNEXT(x, type = 1)
ggiNEXT(x, type = 2)
ggiNEXT(x, type = 3)
