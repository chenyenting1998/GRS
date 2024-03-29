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
source("source/box.cox.chord.R")

# Define internal function
plot_taxa_boxplot <- function(data, 
                              value = "Value", 
                              exp = 1,
                              filter_zero = TRUE){
  # remove metadata columns
  data_preped <- 
    data %>% 
    ungroup() %>% 
    select(- Month, - Station, - Deployment, - Tube)
  # take exp
  data_preped <- data_preped^exp
  # convert to long data
  data_preped <- 
    data_preped %>%
    mutate(Month = data$Month) %>% 
    pivot_longer(cols = -Month,
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
    geom_point(aes(x = Taxon,   y = Value, color = Month), position = "jitter") +
    scale_color_manual(values = cruise_color) +
    ylab(paste0("log[10] (", value, ")")) +
    coord_flip() +
    theme_bw()
  return(figure)
}

plot_taxa_den <- function(data, 
                           value = "Value", 
                           exp = 1,
                           filter_zero = TRUE){
  # remove metadata columns
  data_preped <- 
    data %>% 
    ungroup() %>% 
    select(- Month, - Station, - Deployment, - Tube)
  
  # take the exp
  data_preped <- data_preped^exp
  
  # convert to long data
  data_preped <-
    data_preped %>% 
    mutate(Month = data$Month) %>% 
    pivot_longer(cols = -Month,
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

add_month <- function(x){
  x$Month <- if_else(x$Cruise == "OR1-1219", "March", "October")
  return(x)
}

#############################################
# 1. Calculate unit area density and biomass 
#############################################
# remove unwanted specimen
x <- 
  macrofauna_biomass %>%
  # exclude GC1 and GS1
  filter(!Station %in% c("GC1", "GS1")) %>% 
  
  # include only head-intact individuals
  filter(Condition %in% c("C", "FH")) %>%
  
  # remove unknown/terrestrial/pelagic individuals
  filter(!Taxon %in% c("Unknown", "Insecta","Calanoida")) %>%
  
  # remove hydroids noted as stalks since no tissue are attached to the specimens
  filter(!(Taxon %in% "Hydrozoa" & Note %in% "Stalk"))

# total taxa
unique(x$Taxon)  # 27 taxa
# percentage
count_per_table <-
  x %>% 
  add_month() %>% 
  group_by(Month, Station, Taxon) %>% 
  summarise(Count_sum = n()) %>% 
  mutate(Count_per = round(Count_sum/ sum(Count_sum) * 100, 2)) %>% 
  select(-Count_sum) %>%
  pivot_wider(names_from = "Taxon", 
              values_from = "Count_per",
              values_fill = 0)

WM_per_table <-
  x %>% 
  add_month() %>% 
  group_by(Month, Station, Taxon) %>% 
  summarise(WM_sum = sum(WM)) %>% 
  mutate(WM_per = round(WM_sum / sum(WM_sum) * 100, 2)) %>% 
  select(-WM_sum) %>%
  pivot_wider(names_from = "Taxon", 
              values_from = "WM_per",
              values_fill = 0)

# composition
composition <- 
  x %>%
  add_month() %>% 
  group_by(Month, Station, Deployment, Tube, Taxon) %>% 
  summarise(Count = n(),
            Biomass = sum(WM))

##############################
# 2. Plot density composition
##############################
# essential data for highlighting reivsited stations
revisited <- 
  data.frame(Station = c("S3", "S5", "S6", "S7"))

# actual density
density_composition <-
  add_coarse_taxa(composition, match_file = rank_den, output = "Taxa") %>%
  ggplot(aes(
    x = Station,
    y = Count / mc_area / 3, # unit area; 3 cores
    fill = Taxa
  )) +
  geom_bar(stat = "identity") +
  facet_grid(~ Month, scales = "free") +
  scale_fill_manual(values = taxa_den_color) +
  ylab(Density~(ind.~m^{-2})) +
  theme_bw()

# percentage
density_percentage_composition <-
  add_coarse_taxa(composition, match_file = rank_den, output = "Taxa") %>%
  ggplot() +
  geom_bar(aes(x = Station, y = Count, fill = Taxa), position = "fill", stat = "identity") +
  # adding red brackets to revisited stations
  # geom_tile(data = revisited, aes(x = Station, y = .5), color = "red", alpha = 0, linewidth = 1.5) +
  facet_grid(Month ~ ., scales = "free", switch = "y") +
  scale_fill_manual(values = taxa_den_color) +
  scale_y_continuous(labels = scales::percent, position = "right") +
  ylab(Density~("%")) + 
  theme_bw() 

##############################
# 3. Plot biomass composition
##############################
# actual wet mass
biomass_composition <-
  add_coarse_taxa(composition, match_file = rank_bio, output = "Taxa") %>%
  ggplot(aes(
    x = Station,
    y = Biomass / mc_area / 1000 / 3, # unit area; mg->g; 3 cores
    fill = Taxa
  )) +
  geom_bar(stat = "identity") +
  facet_grid( ~ Month, scales = "free") +
  scale_fill_manual(values = taxa_bio_color) +
  ylab(Wet ~ mass ~ (g~m ^ {-2})) +
  theme_bw()

# percentage wet mass
biomass_percentage_composition <- 
  add_coarse_taxa(composition, match_file = rank_bio, output = "Taxa") %>% 
  ggplot() +
  geom_bar(aes(x = Station, y = Biomass, fill = Taxa), position = "fill", stat = "identity") +
  # adding red brackets to revisited stations
  # geom_tile(data = revisited, aes(x = Station, y = .5), color = "red", alpha = 0, linewidth = 1.5) +
  facet_grid(Month ~ ., scales = "free", switch = "y") +
  scale_fill_manual(values = taxa_bio_color)+
  scale_y_continuous(labels = scales::percent, position = "right") +
  ylab(Wet~mass~("%")) +
  theme_bw()

# density_percentage_composition - biomass_percentage_composition

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

#############
# 5. Heat map
#############
count_heatmap <-
  count_wide %>% 
  pivot_longer(cols = -(1:4), names_to = "Taxon", values_to = "Count") %>% 
  ggplot() +
  geom_tile(aes(x = Station, y = factor(Taxon, rev(rank_den$Taxon)), fill = log10(Count + 1))) +
  viridis::scale_fill_viridis() +
  facet_grid(~Month, scales = "free") +
  ylab("Taxon") +
  theme_bw()

biomass_heatmap <- 
  biomass_wide %>% 
  pivot_longer(cols = -(1:4), names_to = "Taxon", values_to = "Biomass") %>% 
  ggplot() +
  geom_tile(aes(x = Station, y = factor(Taxon, rev(rank_bio$Taxon)), fill = log10(Biomass + 1))) +
  viridis::scale_fill_viridis() +
  facet_grid(~Month, scales = "free") +
  ylab("Taxon") +
  theme_bw()

########################
# 6. Data transformation
########################
# yield Box-Cox-chord transformation results
# Note that n < 3*p, Dagnelie's Test too liberal 
# if p > 0.05, the Dagnelie's test result is trustworthy
count_BCD <- as.data.frame(BCD(count_wide[-(1:4)])) # exp = 0.3
biomass_BCD <- as.data.frame(BCD(biomass_wide[-(1:4)])) # exp = 0.1

# plot
plot(count_BCD$BC.exp, count_BCD$`BC.chord_p-val`)
plot(biomass_BCD$BC.exp, biomass_BCD$`BC.chord_p-val`)
# write xl
write_xlsx(count_BCD, path = "table/BCc/count_BCD.xlsx")
write_xlsx(biomass_BCD, path = "table/BCc/biomass_BCD.xlsx")

# set exponent 
count_exp <- 0.3
biomass_exp <- 0.1

##########
# 6. Plot
##########
# wide data transformed following Anderson et al.(2006)
# plot the frequency of # of individuals by each taxa (i.e., matrix entries)
hist(as.vector(as.matrix(count_wide[-c(1:4)])),
     main = "Histogram of taxa count",
     xlab = "Number of individuals")

# species-specific density plot
count_taxa_denplot <- plot_taxa_den(count_wide, "Count", exp = count_exp)
biomass_taxa_denplot <- plot_taxa_den(biomass_wide, "Biomass", exp = biomass_exp)

# species_specific boxplot
count_taxa_boxplot <- plot_taxa_boxplot(count_wide, "Count", exp = count_exp)
biomass_taxa_boxplot <- plot_taxa_boxplot(biomass_wide, "Biomass", exp = biomass_exp)

##############
# Data output
##############
save(composition, file = "data/macrofauna composition.RData")
save(density_percentage_composition, biomass_percentage_composition, file = "data/composition_figure.RData")
save(count_wide, biomass_wide, file = "data/wide_data.Rdata")
save(count_exp, biomass_exp, file = "data/BC_exp.Rdata")

##################
# Figure output
##################
ggsave("figure/barplot/density_composition.png", plot = density_composition)
ggsave("figure/barplot/biomass_composition.png", plot = biomass_composition)
ggsave("figure/barplot/density_percentage_composition.png", plot = density_percentage_composition)
ggsave("figure/barplot/biomass_percentage_composition.png", plot = biomass_percentage_composition)
ggsave("figure/boxplot/count_taxa_boxplot.png", plot = count_taxa_boxplot)
ggsave("figure/boxplot/biomass_taxa_boxplot.png", plot = biomass_taxa_boxplot)
ggsave("figure/denplot/count_taxa_denplot.png", plot = count_taxa_denplot)
ggsave("figure/denplot/biomass_taxa_denplot.png", plot = biomass_taxa_denplot)
ggsave("figure/heatmap/count_heatmap.png", plot = count_heatmap)
ggsave("figure/heatmap/biomass_heatmap.png", plot = biomass_heatmap)

