###################################
##### Environmental condition #####
###################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(GGally)
library(GRSmacrofauna)

# load variables
load("data/cruise_color.RData")

# define env variables
env_var <- 
  c("Depth", "DRM", "Temperature", "Salinity",
  "Sigma-Theta", "Density", "Oxygen", "Fluorescence",
  "Transmission", "d50", "Clay", "Silt", "Sand",
  "CN", "TOC", "TN", "WC", "TOCd13", "Chla", "Porosity")

# env long
env_l <- 
  env %>% 
  pivot_longer(cols = all_of(env_var),
               names_to = "Variables",
               values_to = "Values")
# define env variable order
env_var_ft <- 
  c("Depth", "DRM", "Temperature", "Salinity",
    "Oxygen", "Fluorescence", "Transmission", 
    "Sigma-Theta", "Density", "d50", "Clay", "Silt", "Sand",
    "Porosity", "TOC", "TN", "CN", "WC", "TOCd13", "Chla")

# define env full names and units
env_var_fn <- 
  c("Depth" = "Depth~(m)", 
  "DRM" = "Distance~to~river~mouth~(km)", 
  "Temperature" = "Temperature~({}^o*C)", 
  "Salinity" = "Salinity~(PSU)",
  # "Sigma-Theta", 
  "Density" = "Density~(kg/m^3)",
  "Oxygen" = "Oxygen~(mg/L)", 
  "Fluorescence" = "Fluorescence~(mg/L)",
  "Transmission" = "Transmission~('%')", 
  "d50" = "Median~grain~size~(mu*m)", 
  "Clay" = "Clay~('%')", 
  "Silt" = "Silt~('%')", 
  "Sand" = "Sand~('%')",
  "CN" = "CN~(ratio)", 
  "TOC" = "TOC~('%')", 
  "TN" = "TN~('%')",
  "WC" = "Water~content~('%')",
  "TOCd13" = "delta^13*C~({}^o/{}[oo])", 
  "Chla" = "Chlorophyll~a~(ng/g)", 
  "Porosity" = "Porosity~('%')")

# define env abbr
env_var_abbr <- 
  c("Depth" = "Depth", 
    "DRM" = "DRM", 
    "Temperature" = "Temp", 
    "Salinity" = "Sal",
    # "Sigma-Theta", 
    "Density" = "Den",
    "Oxygen" = "Oxy", 
    "Fluorescence" = "Fluo",
    "Transmission" = "Trans", 
    "d50" = "D50", 
    "Clay" = "Clay", 
    "Silt" = "Silt", 
    "Sand" = "Sand",
    "CN" = "CN", 
    "TOC" = "TOC", 
    "TN" = "TN",
    "WC" = "WC",
    "TOCd13" = "delta^13*C", 
    "Chla" = "Chla", 
    "Porosity" = "Por")

# plot qqplot ----
env_qqplot <-
  ggplot(env_l, aes(sample = Values))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~Variables, scales = "free") +
  scale_color_manual(values = cruise_color) +
  theme_bw()
  
ggsave("figure/env_qqplot.png", plot = env_qqplot)

# plot boxplot ----
env_boxplot <-
  env_l %>% 
  filter(Variables != "Sigma-Theta") %>% 
  ggplot(aes(x = Cruise, y = Values, color = Cruise))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = "jitter") +
  facet_wrap(~factor(Variables, env_var_ft), scales = "free", labeller = as_labeller(env_var_fn, label_parsed)) +
  scale_color_manual(values = cruise_color) +
  theme_bw() +
  theme(legend.position = c(0.91, 0.1),
        strip.text = element_text(size = 8))

ggsave("figure/env_boxplot.png", plot = env_boxplot,
       width = 8,
       height = 6,
       scale = 1.2)

# plot correlation matrix
png(filename = "figure/env_corr_all.png",
    width = 1600,
    height = 800)
corrplot::corrplot(cor(env[,env_var]),
                   addCoef.col = TRUE)
dev.off()

# Variables were selected via prior knowledge and collinearity of |r| > 0.7
# Depth and DRM are dropped since I am interested in analyzing 
#   the relationship between in situ environment and the biota.
# Oxygen was removed since no bottom water in the study area is hypoxic
# Salinity was removed as salinity levels at a certain timepoint does not affect macrofauna.
# Sigma-theta and density were removed as density is not important to macrofauna
# I use d50 as the sole predictor for sediment granuolometry; 
#   clay, silt, sand were all removed.
# Water content was removed due to strong collinearity with porosity.
# Transmission was omitted due to strong collinearity with temperature.
# TN is removed due to strong collineraity with porosity.

# remaining variables are following
# select specific variables ------
env_sel <- 
  c("Temperature", "Fluorescence","d50", 
    "CN", "TOC", "Porosity", "Chla")

# plot corrplot
png(filename = "figure/env_corr_sel.png",
    width = 1600,
    height = 800)
corrplot::corrplot(cor(env[,env_sel]),
                   addCoef.col = TRUE)
dev.off()

# plot PCA ----
# scale data to zero mean and unit standard deviation
env_sel_scale <- cbind(env[1:4], scale(env[,env_sel]))

# plot scaled env distribution
env_sel_scale %>% 
  pivot_longer(cols = -(1:4),
               names_to = "Variables",
               values_to = "Values") %>% 
  ggplot(aes(x = Variables, y = Values)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Cruise), position = "jitter") +
  scale_color_manual(values = cruise_color) +
  theme_bw()

# calculate dist
env_dist <- vegdist(env_sel_scale[-c(1:4)], method = "euclidean")

# permanova
set.seed(1)
env_disp <- betadisper(env_dist, env_sel_scale$Cruise)
env_permdisp <- permutest(env_disp, permutations = 9999)

env_perm <- 
  adonis2(env_dist ~ Cruise,
          data = env_sel_scale,
          permutations = 9999)
env_perm
# similar distribution between cruises
# significant differences between cruises

# calculate PCA
env_pca <- rda(scale(env[,env_sel]))

# get scores
# attach stations 
env_pca_sites <- 
  scores(env_pca, scaling = "sites")$sites %>% 
  cbind(env[c("Cruise", "Station")]) %>% 
  as.data.frame()

env_pca_species <- 
  scores(env_pca, scaling = "sites")$species %>% 
  as.data.frame()

env_pca_species$abbr <- env_var_abbr[match(rownames(env_pca_species), names(env_var_abbr))]

# extract eigenvalue
env_pca_eig <- round((env_pca$CA$eig/sum(env_pca$CA$eig))*100, 2)

# plot screeplot
data_scree_plot <- 
  data.frame(pc = names(env_pca$CA$eig),
             eigen = env_pca$CA$eig,
             var= env_pca_eig)

ggplot(data_scree_plot) +
  geom_bar(aes(x = pc, y = var), stat = "identity") +
  xlab("Principal component") +
  ylab(Variance~(`%`)) +
  theme_bw()

# plot PCA
env_pca_plot <- 
  ggplot()+
  # plot 
  stat_ellipse(data = env_pca_sites, 
               aes(x = PC1, y = PC2, color = Cruise, fill = Cruise), 
               type = "norm", geom = "polygon",
               size = 1.5,
               alpha = 0.15,
               level = 0.95) +
  
  # plot env variables
  geom_segment(data = env_pca_species,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               size = .8, color = "blue")+
  geom_label(data = env_pca_species, 
             aes(x = PC1, y = PC2, label = abbr)) +

  # plot stations
  geom_label(data = env_pca_sites, 
             aes(x = PC1, y = PC2, color = Cruise, label = Station)) +

  # change axis label
  xlab(paste0("PC1 (", env_pca_eig[1], "% of total variance explained)")) +
  ylab(paste0("PC2 (", env_pca_eig[2], "% of total variance explained)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()

ggsave("figure/env_pca_plot.png", plot = env_pca_plot,
       width = 8,
       height = 6,
       scale = 0.8)

save(env_sel, env_sel_scale, file = "data/env_sel.RData")
