###################################
###### Macrofauna composition #####
###################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(GGally)
library(GRSmacrofauna)
library(writexl)

# load
load("data/taxa_rank.Rdata")
load("data/taxa_color.Rdata")
load("data/cruise_color.Rdata")

# variables
sw <- 1.13 # specific weight
mc_area <- (0.1/2)^2*pi # multiple corer area; diameter = 10 cm

# calculate unit area biomass and density 
comp <- 
  size %>% 
  filter(Condition %in% c("C", "FH")) %>% 
  filter(!Taxon %in% c("Unknown", "Calanoida")) %>% 
  filter(!(Taxon %in% "Hydrozoa" & Note %in% "Stalk")) %>% 
  group_by(Cruise, Station, Deployment, Tube, Taxon) %>% 
  summarize(Abundance = n(),
            Biomass = sum(Size * sw))
save(comp, file = "data/composition.RData")

# plot composition ----
density_comp <- 
  add_coarse_taxon(comp, match_file = rank_den, output = "Taxa") %>% 
  ggplot(aes(x = Station, y = Abundance, fill = Taxa)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~Cruise, scales = "free") +
  scale_fill_manual(values = taxa_den_color)+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()

biomass_comp <- 
  add_coarse_taxon(comp, match_file = rank_bio, output = "Taxa") %>% 
  ggplot(aes(x = Station, y = Biomass, fill = Taxa)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~Cruise, scales = "free") +
  scale_fill_manual(values = taxa_bio_color)+
  scale_y_continuous(labels = scales::percent) +
  theme_bw()

ggsave("figure/density_composition.png", 
       plot = density_comp)
ggsave("figure/biomass_composition.png", 
       plot = biomass_comp)

# composition data
# set up wide data ----
# calculate density to wide dataframe
den_wide <-
  comp %>% 
  select(-Biomass) %>% 
  pivot_wider(names_from = Taxon, values_from = Abundance, values_fill = 0)

# calculate biomass to wide dataframe
bio_wide <-
  comp %>% 
  select(-Abundance) %>% 
  pivot_wider(names_from = Taxon, values_from = Biomass, values_fill = 0)

save(den_wide, bio_wide, file = "data/wide_data.Rdata")
save(den_wide, bio_wide, file = "data/wide_data.Rdata")

# plot nMDS ----
# plot taxon density distribution
# all species histogram
den_wide[-c(1:4)] %>% as.matrix() %>% as.vector %>% hist()
# species_specific boxplot
den_wide %>% 
  pivot_longer(col = -(1:4),
               names_to = "Taxon",
               values_to = "Abundance") %>% 
  ggplot() +
  geom_boxplot(aes(x = Taxon, y = Abundance^ .25), outlier.shape = NA) +
  geom_point(aes(x = Taxon,   y = Abundance^ .25, color = Cruise)) +
  scale_color_manual(values = cruise_color) +
  coord_flip() +
  theme_bw()

# calculate density dist
den_dist <- 
  den_wide[-c(1:4)]^ .25 %>% 
  vegdist(method = "bray")

# calculate density nmds
den_nmds <- metaMDS(den_dist, autotransform = FALSE)
stressplot(den_nmds)

# extract score and stress
den_nmds_site <- den_nmds$points %>% data.frame %>% cbind(den_wide[1:4])
den_nmds_stress <- den_nmds$stress %>% round(2)

# plot nmds plot
den_nmds_plot <-
  ggplot(den_nmds_site, aes(x = MDS1, y = MDS2, color = Cruise, label = Station)) +
  stat_ellipse(aes(fill = Cruise), 
               type = "norm", geom = "polygon",
               size = 1.5,
               alpha = 0.15,
               level = 0.95) +
  geom_label() +
  annotate(geom = "text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.5, label = paste0("Stress: ", den_nmds_stress)) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

ggsave("figure/density_nmds_plot.png", plot = den_nmds_plot)

# plot biomass
bio_wide[-c(1:4)]^0.5 %>% as.matrix() %>% as.vector %>% hist()

# species_specific boxplot
bio_wide %>% 
  pivot_longer(col = -(1:4),
               names_to = "Taxon",
               values_to = "Biomass") %>% 
  ggplot() +
  geom_boxplot(aes(x = Taxon, y = Biomass^0.25), outlier.shape = NA) +
  geom_point(aes(x = Taxon, y = Biomass^0.25, color = Cruise), position = "jitter") +
  scale_color_manual(values = cruise_color) +
  # facet_wrap(~Cruise) +
  coord_flip() +
  theme_bw()

# calculate biomass dist
bio_dist <- 
  bio_wide[-c(1:4)]^0.5 %>%
  vegdist(method = "bray")

# calculate biomass nmds
bio_nmds <- metaMDS(bio_dist, autotransform = FALSE)
stressplot(bio_nmds)

# extract score and stress
bio_nmds_site <- bio_nmds$points %>% data.frame %>% cbind(bio_wide[1:4])
bio_nmds_stress <- bio_nmds$stress %>% round(2)

# plot nmds plot
bio_nmds_plot <- 
  ggplot(bio_nmds_site, aes(x = MDS1, y = MDS2, color = Cruise, label = Station)) +
  stat_ellipse(aes(fill = Cruise), 
               type = "norm", geom = "polygon",
               size = 1.5,
               alpha = 0.15,
               level = 0.95) +
  geom_label() +
  annotate(geom = "text", x = Inf, y = Inf, hjust = 1.2, vjust = 1.5, label = paste0("Stress: ", bio_nmds_stress)) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

ggsave("figure/biomass_nmds_plot.png",plot = bio_nmds_plot)

save(den_dist, bio_dist, file = "data/dist_data.Rdata")