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
library(ggnewscale)
library(patchwork)
library(GRSmacrofauna)
library(Polychrome)
library(writexl)

# Load data
load("data/env.Rdata")
load("data/env_variables.Rdata")
load("data/env_selected.Rdata")
load("data/cruise_color.Rdata")
load("data/polychaete_rank.Rdata")
load("data/polychaete_family_color.Rdata")

# Load functions
source("source/box.cox.chord.R")
source("source/plot_pca.R")
source("source/plot_rda.R")
source("source/plot_scree.R")
source("source/extract_goodness.R")

ord.rsquare <- 
  function(fullmodel, backward_model){
    data.frame(Model = c("Full model", "Backward selected"),
               R.squared = c(RsquareAdj(fullmodel)$r.squared, 
                             RsquareAdj(backward_model)$r.squared),
               Adj.R.squared = c(RsquareAdj(fullmodel)$adj.r.squared, 
                                 RsquareAdj(backward_model)$adj.r.squared))
    
  }

add_month <- function(x){
  x$Month <- if_else(x$Cruise == "OR1-1219", "March", "October")
  return(x)
}

# subset
pol <- 
  macrofauna_biomass %>% 
  filter(Taxon == "Polychaeta")
# attach NA with unknown
pol$Family[is.na(pol$Family)] <- "Unknown"

# how many specimen can / cannot be classified to family level
sum(pol$Family != "Unknown"); sum(pol$Family == "Unknown")

# How many species by percentage are identified to family
sum(pol$Family != "Unknown") / nrow(pol) * 100

# how many specimens are / are not classified to genus level
sum(!is.na(pol$Genus)) ; sum(is.na(pol$Genus))

# subset all speciemens that can be classified to family level
pol %>% filter(Family != "Unknown") %>% nrow()

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
  # remove NAs and unknowns
  filter(!is.na(Family)) %>% 
  filter(!Family %in% c("Unknown", "unknown"))
  
  
pol_c <- pol %>% filter(Condition %in% c("C", "FH"))
pol_f <- pol %>% filter(!Condition %in% c("C", "FH"))

# what fragmented taxa are unique to cephalic intact ones?
unique(pol_f$Family)[!unique(pol_f$Family) %in% unique(pol_c$Family)] 
# How abundant is fabriciidae?
pol %>% filter(Family == "Fabriciidae") %>% nrow()

# What fragmented taxa are unique from the cephalic intact ones?
unique(pol_f$Genus)[!unique(pol_f$Genus) %in% unique(pol_c$Genus)] 

# Only one family (which has one individual) is omitted from the dataset due to lack of cephalic intact specimen a 
# No genuses were omitted due to the cephalic intact regime
# We therefore consider that our data filtering criteria can represent the polychaete composition.

# Note: there is a possibility that a fragmented specimen belongs to 
# a family/genus that does not appear in complete/head-intact specimen.

pol <-
  macrofauna_biomass %>% 
  # cruise to month
  filter(Taxon == "Polychaeta") %>% 
  # remove canyon
  filter(!Station %in% c("GC1",  "GS1")) %>% 
  # remove NAs and unknowns
  filter(!is.na(Family)) %>% 
  filter(!Family %in% c("Unknown", "unknown")) %>% 
  # filter head intact specimens
  filter(Condition %in% c("C", "FH")) %>% 
  # use one column to describe taxonomic resolution
  mutate(Taxa = if_else(!is.na(Genus), paste0("G. ", Genus), paste0("F. ", Family))) %>% 
  # data cleaning
  add_month() %>% 
  relocate(Month, .before = Cruise) %>% 
  select(-Family, -Genus, -Cruise) 

#############
# 2. Barplots
#############
# family-level long dataframe
pol_long <- 
  pol %>% 
  group_by(Month, Station, Deployment, Tube, Taxa) %>% 
  summarise(Count = n()) %>% 
  ungroup()

# set up functions
add_coarse_taxa <- function(data, match_file, output = "Taxa"){
  # collecting the two columns needed for this function
  # Family and the columne named "output"
  extract <- 
    data.frame(Taxa = match_file$Taxa, 
               match_file[,colnames(match_file) == output])
  # matching rows of two dataframes
  comparing <- match(data$Taxa, extract$Taxa)
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
    add_coarse_taxa(rank, output = "Dominant")
    
  plot <- 
    ggplot(data_rank, 
           aes(x = .data[[xx]], y = .data[[yy]], fill = Dominant)) +
    geom_bar(stat = "identity", position = position) +
    facet_grid(~Month, scales = "free") +
    scale_fill_manual(values = color) +
    theme_bw()
  return(plot)
}

# den
pol_den_identity <- plot_comp(pol_long, rank_den, "Station", "Count", "stack", taxa_den_color)
pol_den_stack    <- plot_comp(pol_long, rank_den, "Station", "Count", "fill", taxa_den_color)

ggsave("figure/barplot/pol_den_identity.png", plot = pol_den_identity, scale = 1.2, width = 6, height = 5)
ggsave("figure/barplot/pol_den_stack.png", plot = pol_den_stack, scale = 1.2, width = 6, height = 5)

############
# 3. Heatmap
############
pol_count_long <-
  pol_long %>% 
  pivot_wider(names_from = Taxa, 
              values_from = Count,
              values_fill = 0) %>% 
  pivot_longer(cols = -(1:4), 
               names_to = "Taxa", 
               values_to = "Count")


# heatmap by actual counts
heatmap_count <- 
  pol_count_long %>% 
  ggplot() +
  geom_tile(aes(x = Station, 
                y = factor(Taxa, levels = rank_den$Taxa), 
                fill = Count)) +
  facet_grid(~Month, scales = "free")+
  scale_fill_gradient(low = "white", high = "red")+
  ylab("Taxa")+
  theme_bw()

# heatmap by percentage
heatmap_percentage <-
  pol_count_long %>%
  group_by(Month, Station,) %>% 
  summarise(Taxa = Taxa,
            Percentage = Count / sum(Count) * 100) %>% 
  ggplot() +
  geom_tile(aes(x = Station, 
                y = factor(Taxa, levels = rank_den$Taxa), 
                fill = Percentage)) +
  facet_grid(~Month, scales = "free")+
  scale_fill_gradient(low = "white", high = "red")+
  ylab('Family') +
  theme_bw()

ggsave("figure/heatmap/heatmap_count.png", plot = heatmap_count, width = 7, height = 8)
ggsave("figure/heatmap/heatmap_percentage.png", plot = heatmap_percentage, width = 7, height = 8)

# ##########
# # 4. iNEXT
# ##########
# setColnames <- function(data, name){
#   colnames(data)
# }
# # data formatting for iNEXT
# pol_wide_station <- 
#   pol_count_long %>%
#   group_by(Cruise, Station, Family) %>% 
#   summarise(Count = sum(Count)) %>% 
#   pivot_wider(names_from = "Family", 
#               values_from = "Count",
#               values_fill = 0) %>% 
#   ungroup()
# 
# pol_wide_station_t <-
#   pol_wide_station %>%
#   select(-Cruise, - Station) %>% 
#   as.matrix() %>% 
#   t() %>% 
#   as.data.frame()
# # set colnames
# colnames(pol_wide_station_t) <- paste0(pol_wide_station$Cruise,
#                                        "_",
#                                        pol_wide_station$Station)
# # run iNEXT  
# set.seed(1)
# x <- iNEXT(pol_wide_station_t, datatype = "abundance", q = c(0, 1, 2))
# 
# # output
# write_xlsx(list(`Raw data` = x$DataInfo, `Asymptote estimation`= x$AsyEst),
#            path = "table/data/pol_iNEXT.xlsx")
# 
# # 
# x$AsyEst
# x$DataInfo$Cruise <- gsub("_.*", "", x$DataInfo$site) 
# x$DataInfo$Station <- gsub(".*_", "", x$DataInfo$site) 
# 
# # abundance vs. number of families
# pol_n_spp <- 
#   ggplot(x$DataInfo, aes(x = n, y = S.obs, color = Cruise, label = Station)) +
#   geom_point() +
#   geom_text_repel(seed = 1, size = 3) +
#   scale_color_manual(values = cruise_color) +
#   scale_x_continuous(limits = c(0, 185)) +
#   scale_y_continuous(limits = c(0,  25)) +
#   xlab(Polychaete~sample~size) +
#   ylab(Number~of~families~observed) +
#   theme_bw() +
#   theme(legend.position = c(0.9, 0.1))
# 
# ggsave("figure/scatterplot/pol_n_spp.png", plot = pol_n_spp, width = 6, height = 5)
# 
# # Number of rare families vs. number of families
# pol_f_spp <- 
#   ggplot(x$DataInfo, aes(x = f1 + f2 + f3 , y = S.obs, color = Cruise, label = Station)) +
#   geom_point() +
#   geom_text_repel(seed = 1, size = 3) +
#   scale_color_manual(values = cruise_color) +
#   scale_x_continuous(limits = c(0, 20)) +
#   scale_y_continuous(limits = c(0, 25)) +
#   xlab(Rare~species~(f1+f2+f3)) +
#   ylab(Number~of~families~observed) +
#   theme_bw() +
#   theme(legend.position = c(0.9, 0.1))
# 
# ggsave("figure/scatterplot/pol_f_spp.png", plot = pol_f_spp, width = 6, height = 5)
# 
# pol_n_f <- 
#   ggplot(x$DataInfo, aes(x = n, y = f1 + f2 + f3, color = Cruise, label = Station)) +
#   geom_point() +
#   geom_text_repel(seed = 1, size = 3) +
#   scale_color_manual(values = cruise_color) +
#   scale_x_continuous(limits = c(0, 185)) +
#   scale_y_continuous(limits = c(0, 20)) +
#   xlab(Polychaete~sample~size) +
#   ylab(Rare~species~(f1+f2+f3)) +
#   theme_bw() +
#   theme(legend.position = c(0.9, 0.1))
# 
# ggsave("figure/scatterplot/pol_n_f.png", plot = pol_n_f, width = 6, height = 5)
# 
# # Rarefaction curve
# station_color <- kelly.colors(10)[-1]
# names(station_color) <- unique(macrofauna_biomass$Station)
# samplesize_rarefaction <- 
#   fortify(x) %>% 
#   mutate(Station = gsub(".*_", "", site),
#          Cruise = gsub("_.*", "", site)) %>% 
#   ggplot()+
#   geom_line(aes(x = x, y = y, color = Station, linetype = Cruise), size = 1.1) +
#   facet_wrap(~order, scales = "free") +
#   xlab(Number~of~individuals) +
#   ylab(Effective~number~of~species) +
#   scale_color_manual(values = station_color) +
#   theme_bw()
# 
# ggsave("figure/pol_samplesize_rarefaction.png",
#        plot = samplesize_rarefaction)
# 
# coverage_rarefaction <- 
#   fortify(x, type = 3) %>% 
#   mutate(Station = gsub(".*_", "", site),
#          Cruise = gsub("_.*", "", site)) %>% 
#   ggplot()+
#   geom_line(aes(x = x, y = y, color = Station, linetype = Cruise), size = 1.1) +
#   facet_wrap(~order, scales = "free") +
#   xlab(Sample~coverage) +
#   ylab(Effective~number~of~species) +
#   scale_color_manual(values = station_color) +
#   theme_bw()
# ggsave("figure/pol_coverage_rarefaction.png",
#        plot = coverage_rarefaction)

#####################
# 5. Composition ####
#####################
# p >> n
# extract taxa > 1 percent
taxa_subset <- rank_den[rank_den$percent > 1,]$Taxa

# wide data
pol_count_wide <-
  pol_long %>% 
  filter(Taxa %in% taxa_subset) %>% 
  pivot_wider(names_from = Taxa, 
            values_from = Count,
            values_fill = 0)

# find best BC.chord value
pol_bcd <- BCD(pol_count_wide[-(1:4)])
write_xlsx(list(pol_BCD = as.data.frame(pol_bcd)), path = "table/BCc/pol_count_BCD.xlsx")
pol_exp <- 0.4

# box cox chord transformation
pol_chord <- 
  pol_count_wide %>% 
  select(-all_of(c("Station", "Deployment", "Tube", "Month"))) %>% 
  box.cox.chord(pol_exp)

################
# PERMANOVA ####
################
# sp
env_sp <- env[, c("Month", "Station", env_variables_spatial)]
env_sp$Depth <- scale(env_sp$Depth)
env_sp$DRM <- scale(env_sp$DRM)
data_temp <- left_join(pol_count_wide, env_sp)

## permanova
set.seed(100)
pol_permanova <- 
  adonis2(pol_chord ~ Depth * Month + DRM * Month + Depth * DRM,
          data = data_temp,
          method = "euclidean",
          permutations = 9999)

write_xlsx(list(PERMANOVA = cbind(" " = rownames(pol_permanova), pol_permanova)),
           path = "table/permanova/pol_permanova.xlsx")

# extract env_sp_metadata
env_sp_metadata <- env[c("Month", "Station", "Depth", "DRM")]
size_range = c(2, 5)
size_depth_breaks <- c(30,50,70,90)
size_drm_breaks = c(15, 20, 25)
stretch = 0.2

# tb-PCA
pol_pca <- rda(pol_chord)

screeplot(pol_pca)

pol_pca_goodness <- extract_goodness(pol_pca, "CA")
pol_pca_goodness_plot <-
  ggplot(pol_pca_goodness) +
  geom_point(aes(x = PC2, 
                 y = Taxon)) +
  xlab("Cummulative variance explained to PC2") +
  theme_bw()


## get output scaling = 1 ####
pol_sc1 <- 
  get_pca_output(pol_pca, 
                 metadata = left_join(pol_count_wide[1:4], env_sp_metadata),
                 scaling = 1)

## get output scaling = 2 ####
pol_sc2 <- 
  get_pca_output(pol_pca, 
                 metadata = left_join(pol_count_wide[1:4], env_sp_metadata),
                 scaling = 2)

## plot pca plots scaling = 1 ####
pol_pca_sc1 <- 
  plot_pca(pol_sc1$pca_sites, 
           pol_sc1$pca_species, 
           scaling = 1, 
           eig_vector = pol_sc1$pca_eig,
           stretch = .2)

ggsave("figure/rda/pol_pca_sc1.png", 
       plot = pol_pca_sc1, 
       scale = 1,
       width = 8,
       height = 6)

## plot pca plots scaling = 2 ####
pol_pca_sc2 <- 
  plot_pca(pol_sc2$pca_sites, 
           pol_sc2$pca_species, 
           scaling = 2, 
           eig_vector = pol_sc2$pca_eig, 
           stretch = 1)

ggsave("figure/rda/pol_pca_sc2.png", 
       plot = pol_pca_sc2, 
       scale = 1,
       width = 8,
       height = 6)

## depth pca ####
pol_pca_sc1_depth <- 
  ggplot() +
  
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # species
  geom_segment(data = pol_sc1$pca_species[pol_sc1$pca_species$Show == TRUE,],
               aes(x = 0, 
                   y = 0,
                   xend = PC1 * stretch, 
                   yend = PC2 * stretch),
               size = .4, 
               color = "purple")+
  geom_label(data = pol_sc1$pca_species[pol_sc1$pca_species$Show == TRUE,],
             aes(x = PC1 * stretch, 
                 y = PC2 * stretch, 
                 label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +
  # sites
  geom_point(data = pol_sc1$pca_sites,
             aes(x = PC1,
                 y = PC2,
                 color = Month,
                 size = Depth)) +
  scale_size_binned("Depth (m)", range = size_range + 1.5, breaks = size_depth_breaks) +
  new_scale("size") +
  geom_text(data = pol_sc1$pca_sites,
            aes(x = PC1,
                y = PC2,
                color = Month,
                size = Depth,
                label = Station),
            color = "white") +
  scale_size_binned("Depth (m)", range = size_range - 1.5, breaks = size_depth_breaks) +
  # change axis label
  xlab(paste0("PC1 (", pol_sc1$pca_eig[1], "% of the total variance)")) +
  ylab(paste0("PC2 (", pol_sc1$pca_eig[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

ggsave("figure/rda/pol_pca_sc1_depth.png", 
       plot = pol_pca_sc1_depth, 
       scale = 1,
       width = 8,
       height = 6)


## pca drm ####
pol_pca_sc1_drm <- 
  ggplot() +
  
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # species
  geom_segment(data = pol_sc1$pca_species[pol_sc1$pca_species$Show == TRUE,],
               aes(x = 0, 
                   y = 0,
                   xend = PC1 * stretch, 
                   yend = PC2 * stretch),
               size = .4, 
               color = "purple")+
  geom_label(data = pol_sc1$pca_species[pol_sc1$pca_species$Show == TRUE,],
             aes(x = PC1 * stretch, 
                 y = PC2 * stretch, 
                 label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +
  # sites
  geom_point(data = pol_sc1$pca_sites,
             aes(x = PC1,
                 y = PC2,
                 color = Month,
                 size = DRM)) +
  scale_size_binned("DRM (km)", range = size_range + 1.5, breaks = size_drm_breaks) +
  new_scale("size") +
  geom_text(data = pol_sc1$pca_sites,
            aes(x = PC1,
                y = PC2,
                color = Month,
                size = DRM,
                label = Station),
            color = "white") +
  scale_size_binned("DRM (km)", range = size_range - 1.5, breaks = size_drm_breaks) +
  # change axis label
  xlab(paste0("PC1 (", pol_sc1$pca_eig[1], "% of the total variance)")) +
  ylab(paste0("PC2 (", pol_sc1$pca_eig[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

ggsave("figure/rda/pol_pca_sc1_drm.png", 
       plot = pol_pca_sc1_drm, 
       scale = 1,
       width = 8,
       height = 6)

# low PCA axis representative yet good PERMANOVA results
# more species can explain the spatiotemporal trend


####################################
# Canonical redundancy analysis ####
####################################
# match env_selected data.frame with biomass_wide
env_selected_expand <- 
  left_join(pol_count_wide, env_selected) %>% 
  select(all_of(env_variables_selected))

# full RDA model
pol_rda_full <- rda(pol_chord ~ ., data = as.data.frame(scale(env_selected_expand)))

# backward selection
pol_rda_back <- ordistep(pol_rda_full, method = "backward", permutations = 9999)

# compare full and reduced model
anova(pol_rda_full) # significant
anova(pol_rda_back) # significant
anova(pol_rda_full, pol_rda_back) # no sig. diff. btw the full and reduced

# rsquared
RsquareAdj(pol_rda_full)
RsquareAdj(pol_rda_back) # adj.r2 increased

# vif
vif.cca(pol_rda_full) # TOC and porosity have high vif
vif.cca(pol_rda_back)

# residual plots
ordiresids(pol_rda_full)
ordiresids(pol_rda_back)


summary(pol_rda_back)

pol_rda_for_goodness <- extract_goodness(pol_rda_back, "CCA")
pol_rda_for_goodness_plot <-
  ggplot(pol_rda_for_goodness) +
  geom_point(aes(x = RDA6, 
                 y = Taxon)) +
  xlab("Cummulative variance explained to RDA6") +
  theme_bw()

## scaling = 1 ####
pol_rda_output_sc1 <- 
  get_rda_output(pol_rda_back, 
                 # left_join(biomass_wide[1:4], env_sp),
                 left_join(pol_count_wide[1:4], env_sp_metadata),
                 env_variables_abbr,
                 scaling = 1)

pol_rda_plot_sc1 <- 
  plot_rda(rda_sites = pol_rda_output_sc1$rda_sites, 
           rda_env = pol_rda_output_sc1$rda_env,
           rda_species = pol_rda_output_sc1$rda_species,
           rda_result = pol_rda_back,
           scaling = 1,
           stretch = 0.2)

ggsave("figure/rda/pol_rda_plot_sc1.png", 
       plot = pol_rda_plot_sc1, 
       scale = 1,
       width = 8,
       height = 6)

# extract reduced model statistics
set.seed(10)
pol_rda_axis <- anova.cca(pol_rda_back, by = "axis", permutations = 9999)
# first two rda axises are sig.
pol_rda_margin <- anova.cca(pol_rda_back, by = "margin", permutations = 9999)
# all variables are sig.
write_xlsx(list(pol_rda_statistics = ord.rsquare(pol_rda_full, pol_rda_back),
                pol_rda_axis = cbind(" " = rownames(pol_rda_axis), pol_rda_axis),
                pol_rda_margin = cbind(" " = rownames(pol_rda_margin), pol_rda_margin)),
           path = "table/rda/pol_rda_anova.xlsx")

## pol_rda_plot_sc1_depth ####
pol_rda_plot_sc1_depth <-
  ggplot() +
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # # plot env
  geom_segment(data = pol_rda_output_sc1$rda_env,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 0),
               size = .5, color = "black")+
  geom_text(data = pol_rda_output_sc1$rda_env,
            aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = abbr),
            size = 3,
            parse = TRUE) +
  # # species vectors
  geom_segment(data = pol_rda_output_sc1$rda_species[pol_rda_output_sc1$rda_species$Show == TRUE,],
               aes(x = 0, y = 0, xend = RDA1 * stretch, yend = RDA2 * stretch),
               arrow = arrow(angle = 0),
               size = .5,
               color = "purple")+
  geom_label(data = pol_rda_output_sc1$rda_species[pol_rda_output_sc1$rda_species$Show == TRUE,],
             aes(x = RDA1 * stretch, y = RDA2 * stretch, label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +
  # plot stations
  geom_point(data = pol_rda_output_sc1$rda_sites,
             aes(x = RDA1,
                 y = RDA2,
                 color = Month,
                 size = Depth)) +
  scale_size_binned("Depth (m)", range = size_range + 1.5, breaks = size_depth_breaks) +
  new_scale("size") +
  geom_text(data = pol_rda_output_sc1$rda_sites,
            aes(x = RDA1,
                y = RDA2,
                label = Station,
                size = Depth),
            color = "white") +
  scale_size_binned("Depth (m)", range = size_range - 1.5, breaks = size_depth_breaks) +
  
  # add R2
  xlab(paste0("RDA1 (", rda_eig_percent(pol_rda_back)[1], "% of the total variance)")) +
  ylab(paste0("RDA2 (", rda_eig_percent(pol_rda_back)[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

ggsave("figure/rda/pol_rda_plot_sc1_depth.png", 
       plot = pol_rda_plot_sc1_depth, 
       scale = 1,
       width = 8,
       height = 6)

## pol_rda_plot_sc1_drm ####
pol_rda_plot_sc1_drm <-
  ggplot() +
  # add v and hline
  geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
  geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
  
  # # plot env
  geom_segment(data = pol_rda_output_sc1$rda_env,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 0),
               size = .5, color = "black")+
  geom_text(data = pol_rda_output_sc1$rda_env,
            aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = abbr),
            size = 3,
            parse = TRUE) +
  # # species vectors
  geom_segment(data = pol_rda_output_sc1$rda_species[pol_rda_output_sc1$rda_species$Show == TRUE,],
               aes(x = 0, y = 0, xend = RDA1 * stretch, yend = RDA2 * stretch),
               arrow = arrow(angle = 0),
               size = .5,
               color = "purple")+
  geom_label(data = pol_rda_output_sc1$rda_species[pol_rda_output_sc1$rda_species$Show == TRUE,],
             aes(x = RDA1 * stretch, y = RDA2 * stretch, label = Taxon),
             size = 3,
             color = "purple",
             label.padding = unit(0.15, "lines")) +
  # plot stations
  geom_point(data = pol_rda_output_sc1$rda_sites,
             aes(x = RDA1,
                 y = RDA2,
                 color = Month,
                 size = DRM)) +
  scale_size_binned("DRM (km)", range = size_range + 1.5, breaks = size_drm_breaks) +
  new_scale("size") +
  geom_text(data = pol_rda_output_sc1$rda_sites,
            aes(x = RDA1,
                y = RDA2,
                label = Station,
                size = DRM),
            color = "white") +
  scale_size_binned("DRM (km)", range = size_range - 1.5, breaks = size_drm_breaks) +
  
  # add R2
  xlab(paste0("RDA1 (", rda_eig_percent(pol_rda_back)[1], "% of the total variance)")) +
  ylab(paste0("RDA2 (", rda_eig_percent(pol_rda_back)[2], "% of the total variance)")) +
  scale_color_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

ggsave("figure/rda/pol_rda_plot_sc1_drm.png", 
       plot = pol_rda_plot_sc1_drm, 
       scale = 1,
       width = 8,
       height = 6)