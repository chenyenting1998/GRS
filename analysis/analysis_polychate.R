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
library(iNEXT)
library(vegan)
library(GGally)
library(ggrepel)
library(patchwork)
library(Polychrome)
library(GRSmacrofauna)
library(writexl)

# Load data
load("data/cruise_color.Rdata")
load("data/polychaete_rank.Rdata")
load("data/polychaete_family_color.Rdata")

# Load functions
source("source/box.cox.chord.R")

pol <- 
  macrofauna_biomass %>% 
  filter(Taxon == "Polychaeta")

pol %>%
  select(Family, Genus) %>% 
  mutate(FG = paste0(pol$Family, "_", pol$Genus)) %>% 
  group_by(Family, Genus, FG) %>% 
  summarize(Count = n()) %>%
  View

#######################
# 1. Specimen condition
#######################
ind_cond. <- 
  ggplot(pol, aes(x = Family, fill = Condition)) +
  geom_bar(position = "stack") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0.99))

wm_cond. <- 
  ggplot(pol, aes(x = Family, y = WM, fill = Condition)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0.99))

ggsave("figure/barplot/pol_condition.png",
       plot = ind_cond./wm_cond.,
       width = 8,
       height = 8)

# only include C and FH 
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

#############
# 2. Barplots
#############
# family-level long dataframe
pol_family_long <- 
  pol %>% 
  group_by(Cruise, Station, Deployment, Tube, Family) %>% 
  summarise(Count = n(),
            WM = sum(WM)) %>% 
  ungroup()

# set up functions
add_coarse_family <- function(data, match_file, output = "Family"){
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
plot_comp <- function(data, rank, xx, yy, position, color){
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
pol_den_identity <- plot_comp(pol_family_long, rank_den, "Station", "Count", "stack", family_den_color)
pol_den_stack    <- plot_comp(pol_family_long, rank_den, "Station", "Count", "fill", family_den_color)
# bio
pol_bio_identity <- plot_comp(pol_family_long, rank_bio, "Station", "WM", "stack", family_bio_color)
pol_bio_stack    <- plot_comp(pol_family_long, rank_bio, "Station", "WM", "fill", family_bio_color)

ggsave("figure/barplot/pol_den_identity.png", plot = pol_den_identity, scale = 1.2, width = 6, height = 5)
ggsave("figure/barplot/pol_den_stack.png", plot = pol_den_stack, scale = 1.2, width = 6, height = 5)
ggsave("figure/barplot/pol_bio_identity.png", plot = pol_bio_identity, scale = 1.2, width = 6, height = 5)
ggsave("figure/barplot/pol_bio_stack.png", plot = pol_bio_stack, scale = 1.2, width = 6, height = 5)

############
# 3. Heatmap
############
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

# heatmap by actual counts
heatmap_count <- 
  pol_count_long %>% 
  ggplot() +
  geom_tile(aes(x = Station, 
                y = factor(Family, levels = rank_den$Family), 
                fill = Count)) +
  facet_grid(~Cruise, scales = "free")+
  scale_fill_gradient(low = "white", high = "red")+
  ylab("Family")+
  theme_bw()

# heatmap by percentage
heatmap_percentage <-
  pol_count_long %>%
  group_by(Cruise, Station,) %>% 
  summarise(Family = Family,
            Percentage = Count / sum(Count) * 100) %>% 
  ggplot() +
  geom_tile(aes(x = Station, 
                y = factor(Family, levels = rank_den$Family), 
                fill = Percentage)) +
  facet_grid(~Cruise, scales = "free")+
  scale_fill_gradient(low = "white", high = "red")+
  ylab('Family') +
  theme_bw()

ggsave("figure/heatmap/heatmap_count.png", plot = heatmap_count, width = 7, height = 8)
ggsave("figure/heatmap/heatmap_percentage.png", plot = heatmap_percentage, width = 7, height = 8)

##########
# 4. iNEXT
##########
setColnames <- function(data, name){
  colnames(data)
}
# data formatting for iNEXT
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
set.seed(1)
x <- iNEXT(pol_wide_station_t, datatype = "abundance", q = c(0, 1, 2))

# output
write_xlsx(list(`Raw data` = x$DataInfo, `Asymptote estimation`= x$AsyEst),
           path = "table/data/pol_iNEXT.xlsx")

# 
x$AsyEst
x$DataInfo$Cruise <- gsub("_.*", "", x$DataInfo$site) 
x$DataInfo$Station <- gsub(".*_", "", x$DataInfo$site) 

# abundance vs. number of families
pol_n_spp <- 
  ggplot(x$DataInfo, aes(x = n, y = S.obs, color = Cruise, label = Station)) +
  geom_point() +
  geom_text_repel(seed = 1, size = 3) +
  scale_color_manual(values = cruise_color) +
  scale_x_continuous(limits = c(0, 185)) +
  scale_y_continuous(limits = c(0,  25)) +
  xlab(Polychaete~sample~size) +
  ylab(Number~of~families~observed) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.1))

ggsave("figure/scatterplot/pol_n_spp.png", plot = pol_n_spp, width = 6, height = 5)

# Number of rare families vs. number of families
pol_f_spp <- 
  ggplot(x$DataInfo, aes(x = f1 + f2 + f3 , y = S.obs, color = Cruise, label = Station)) +
  geom_point() +
  geom_text_repel(seed = 1, size = 3) +
  scale_color_manual(values = cruise_color) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0, 25)) +
  xlab(Rare~species~(f1+f2+f3)) +
  ylab(Number~of~families~observed) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.1))

ggsave("figure/scatterplot/pol_f_spp.png", plot = pol_f_spp, width = 6, height = 5)

pol_n_f <- 
  ggplot(x$DataInfo, aes(x = n, y = f1 + f2 + f3, color = Cruise, label = Station)) +
  geom_point() +
  geom_text_repel(seed = 1, size = 3) +
  scale_color_manual(values = cruise_color) +
  scale_x_continuous(limits = c(0, 185)) +
  scale_y_continuous(limits = c(0, 20)) +
  xlab(Polychaete~sample~size) +
  ylab(Rare~species~(f1+f2+f3)) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.1))

ggsave("figure/scatterplot/pol_n_f.png", plot = pol_n_f, width = 6, height = 5)

# Rarefaction curve
station_color <- kelly.colors(10)[-1]
names(station_color) <- unique(macrofauna_biomass$Station)
samplesize_rarefaction <- 
  fortify(x) %>% 
  mutate(Station = gsub(".*_", "", site),
         Cruise = gsub("_.*", "", site)) %>% 
  ggplot()+
  geom_line(aes(x = x, y = y, color = Station, linetype = Cruise), size = 1.1) +
  facet_wrap(~order, scales = "free") +
  xlab(Number~of~individuals) +
  ylab(Effective~number~of~species) +
  scale_color_manual(values = station_color) +
  theme_bw()

ggsave("figure/pol_samplesize_rarefaction.png",
       plot = samplesize_rarefaction)

coverage_rarefaction <- 
  fortify(x, type = 3) %>% 
  mutate(Station = gsub(".*_", "", site),
         Cruise = gsub("_.*", "", site)) %>% 
  ggplot()+
  geom_line(aes(x = x, y = y, color = Station, linetype = Cruise), size = 1.1) +
  facet_wrap(~order, scales = "free") +
  xlab(Sample~coverage) +
  ylab(Effective~number~of~species) +
  scale_color_manual(values = station_color) +
  theme_bw()
ggsave("figure/pol_coverage_rarefaction.png",
       plot = coverage_rarefaction)

#####################
# 5. 
#####################
