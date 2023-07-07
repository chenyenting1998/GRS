##################################
# Gaoping River-shelf environment
##################################
# Description
# Analyze the environmental condition of the GRS

# Author: Yen-Ting Chen
# Date of creation: Unknown
# Date of last modification: 2023/07/06

#######################
# Set up environment
#######################
# clean environment
rm(list = ls())

# Load packages
library(dplyr) # data manipulation
library(tidyr)
library(vegan) # data analysis
library(rstatix) # pipe-friendly stats package
library(ggplot2) # visualization
library(GGally)
library(ggrepel)
library(writexl) # xlsx
library(GRSmacrofauna)  # data package


# load variables
load("data/cruise_color.RData")
load("data/env_variables.RData")

################################
# 1. Convert env to long format
################################
# env long
env_long <- 
  env %>% 
  filter(!Station %in% c("GC1", "GS1")) %>% 
  pivot_longer(cols = all_of(env_variables),
               names_to = "Variables",
               values_to = "Values")

#####################
# 2. Data exploration
#####################
# qqplot
env_qqplot <-
  ggplot(env_long, aes(sample = Values))+
  geom_qq()+
  geom_qq_line()+
  facet_wrap(~Variables, scales = "free") +
  scale_color_manual(values = cruise_color) +
  theme_bw()
  
ggsave("figure/env_qqplot.png", scale = 1.5, plot = env_qqplot)

# Correlation matrix
png(filename = "figure/env_corr_all.png",
    width = 1600,
    height = 800)

env_corr <- Hmisc::rcorr(as.matrix(env[,env_variables]))
corrplot::corrplot(env_corr$r,
                   # corrplot style
                   method = "ellipse",
                   type = "upper",
                   diag = FALSE,
                   addCoef.col = TRUE,
                   # significance
                   p.mat = env_corr$P,
                   insig = "blank")
dev.off()

# boxplot
env_boxplot <-
  env_long %>% 
  filter(!Variables %in% c("Sigma-Theta", "WC")) %>% 
  ggplot(aes(x = Cruise, y = Values, color = Cruise))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = "jitter") +
  facet_wrap( ~ factor(Variables, env_variables), 
             scales = "free") +
             # labeller = as_labeller(env_var_fn, label_parsed)) +
  scale_color_manual(values = cruise_color) +
  theme_bw() +
  theme(legend.position = c(0.91, 0.1),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 8))

ggsave("figure/env_boxplot.png", 
       plot = env_boxplot,
       width = 8,
       height = 6,
       scale = 1.2)

###################
# 3. Kruskal-Wallis 
###################
# set up an empty object
kruskal_table <- NULL
?t_test
x <- 
  env_long %>% 
  group_by(Variables) %>% 
  wilcox_test(Values ~ Cruise,
              comparisons = list(c("OR1-1219", "OR1-1242")))
View(x)
# for-loop kruskal.test() to all the variables
for(variable in env_variables){
  # write a formula for kruskal.test()
  kruskal_formula <- as.formula(paste0(variable,"~", "Cruise"))
  
  # run test
  output <- kruskal.test(kruskal_formula, data = env)
  
  # extract test result
  result <- 
    data.frame("Variable" = gsub(" .*", "", output$data.name),
               "df" = output$parameter,
               "Chi-squared" = output$statistic,
               "p value" = output$p.value,
               row.names = which(env_variables == variable))

  # stitch the result with the previous ones
  kruskal_table <- rbind(kruskal_table, result)
}

# output
write_xlsx(kruskal_table, "table/env_kruskal_table.xlsx")

###########################
# 6. PERMANOVA and PERMDISP
###########################
set.seed(10)
# run permanova
env_permanova <- 
  adonis2(scale(env[,colnames(env) %in% env_variables]) ~ Cruise,
  # adonis2(scale(env_selected) ~ Cruise,
          data = env, 
          method = "euclidean",
          permutations = 9999)
# run permdisp
env_disp <-
  betadisper(vegdist(env[, colnames(env) %in% env_variables],
                     method = "euclidean", scale = TRUE), 
             group = env$Cruise)
# env_disp <- betadisper(vegdist(env_selected, method = "euclidean", scale = TRUE), group = env$Cruise)
env_permdisp <- permutest(env_disp, permutations = 9999)

# check output
env_permanova
env_permdisp

##################################
# 5. Principle coordinate analysis
##################################
env_pca <- rda(scale(env[, colnames(env) %in% env_variables]))

# extract sites
env_pca_sites <- 
  scores(env_pca, scaling = 1)$sites %>% 
  as.data.frame() %>% 
  cbind(env[c("Cruise", "Station")])

# extract species (env var.)
env_pca_species <- 
  scores(env_pca, scaling = 1)$species %>%
  as.data.frame()

# change names to abbr
env_pca_species$abbr <- env_variables_abbr[match(rownames(env_pca_species), names(env_variables_abbr))]

# extract eig
env_pca_eig <- 
  round(eigenvals(env_pca) / sum(eigenvals(env_pca)) * 100, 2) %>% 
  as.vector()

env_pca_plot <-
  ggplot() +
  # plot variable segments
  geom_segment(data = env_pca_species,
               aes(x = PC1, y = PC2, xend = 0, yend = 0)) +
  # plot variable text
  geom_text(data = env_pca_species,
            aes(x = PC1 * 1.09, 
                y = PC2 * 1.09, 
                label = abbr),
            parse = TRUE) +
  # plot stations
  geom_point(data = env_pca_sites, 
             aes(x = PC1, y = PC2, color = Cruise)) +
  geom_label(data = env_pca_sites,#[-c(2:3),], 
             aes(x = PC1, y = PC2, label = Station, color = Cruise)) +
  # geom_label_repel(data = env_all_pca_sites[c(2:3),],
  #                  aes(x = PC1, y = PC2, label = Station, color = Cruise),
  #                  min.segment.length = 0) +
  # change axis label
  xlab(paste0("PC1 (", env_pca_eig[1], "% of total variance explained)")) +
  ylab(paste0("PC2 (", env_pca_eig[2], "% of total variance explained)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed()+
  theme_bw() +
  theme(legend.position = c(0.1,0.1))

ggsave("figure/env_pca_plot.png", 
       plot = env_pca_plot,
       width = 8,
       height = 6,
       scale = 1)

###############################
# 6. Select interested variables
###############################
# omit variables
omit <- c("Salinity", "Density", "SigmaTheta", "Oxygen",
          "Clay", "Silt", "Sand", "WC")
omit_space <- c("Depth", "DRM")
# Reasoning:
#  Salinity: snapshots of salinity change cannot inform us the variability of salinity.
#  Density: might not be an important factor
#  Oxygen: bottom water is well oxygenated.
#  clay, silt, sand: D50 was used as the sole granulometry proxy 
#  Water content: strong collinearity with porosity.

# subset env
env_selected <- env[,!colnames(env) %in% c(env_metadata, omit_space, omit)]


# omit some clearly correlated variables
corrplot::corrplot(cor(env[,env_var[!env_var %in% c(omit, omit_space)]]),
                   addCoef.col = TRUE)
# pair plot (with reduced variables)
pairplot_ggplot <-
  env %>% 
  ggpairs(columns = c(env_var[!env_var %in% c(omit)]),
          legend = 1,
          aes(color = Cruise, fill = Cruise)) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  theme_bw()

ggsave("figure/env_pairs.png", 
       pairplot_ggplot, 
       width = 16,
       height = 9,
       scale = 1.2)

# Depth and DRM are dropped since I am interested in analyzing 
#   the relationship between in situ environment and the biota.

save(env_selected, file = "data/env_sel.RData")
