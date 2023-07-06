#############################################
# Macrofauna abundance, biomass, and drivers 
#############################################
# Description
# Create plots and tables for macrofauna density and biomass 

# Author: Yen-Ting Chen
# Date of creation: 2023/07/05
# Date of last modification: 2023/07/05

#######################
# Set up environment
#######################
# clean environment
rm(list = ls())

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(GGally)
library(ggrepel)
library(GRSmacrofauna)
library(writexl)

# Load data
load("data/macrofauna composition.Rdata")
load("data/cruise_color.Rdata")

# define multiple corer area; diameter = 10 cm
mc_area <- (0.1 / 2) ^ 2 * pi

###############################################
# 1. Calculate unit area abundance and biomass
###############################################
standingstock <-
  composition %>%
  group_by(Cruise, Station, Deployment, Tube) %>%
  summarize(
    Abundance = sum(Abundance) / mc_area,
    Biomass = sum(Biomass) / 1000 / mc_area # mg to g
  )
# output .xlsx table
standingstock_table <-
  standingstock %>% 
  ungroup() %>% 
  group_by(Cruise, Station) %>% 
  summarize(Abundance_mean = mean(Abundance),
            Abundance_sd = sd(Abundance),
            Biomass_mean = mean(Biomass),
            Biomass_sd = sd(Biomass)) 

write_xlsx(standingstock_table, "table/abundance_biomass.xlsx")

############################################################
# 2. Plot standing stock relationships with station and DRM
############################################################
# set up data.frame
standingstock_env <-
  standingstock %>% 
  pivot_longer(cols = c("Abundance", "Biomass"),
              names_to = "Variable",
              values_to = "Value") %>% 
  group_by(Cruise, Station, Variable) %>% 
  summarize(mean = mean(Value),
            sd = sd(Value)) %>% 
  left_join(env)
  

# vs station (discret var.)
standingstock_vs_Station <- 
  standingstock_env %>% 
  ggplot(aes(x = Station, y = mean, ymin = mean - sd, ymax = mean + sd, fill = Cruise, color = Cruise)) +
  geom_pointrange(position = position_dodge(width = .5)) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  facet_grid(Variable~., scales = "free") +
  ylab("Value") +
  theme_bw()

# vs DRM (continuous var.)
standingstock_vs_DRM <-
  standingstock_env %>% 
  ggplot(aes(x = DRM, y = mean, ymin = mean - sd, ymax = mean + sd, label = Station, fill = Cruise, color = Cruise)) +
  geom_smooth(method = "lm", se = TRUE)+
  geom_pointrange() +
  geom_text_repel(seed = 1) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  facet_grid(Variable~Cruise, scales = "free") +
  ylab("Value") +
  theme_bw()

# fit abundance and biomass to DRM ----
standing_stock_spatial <- 
  standing_stock %>% 
  left_join(ou) %>% 
  full_join(env[c("Cruise", "Station", "Depth", "DRM", "Porosity")]) %>% 
  pivot_longer(cols = c("Abundance", "Biomass", "TOU"),
               names_to = "Variables",
               values_to = "Values")

ab_fn <- 
  c("OR1-1219" = "OR1-1219",
    "OR1-1242" = "OR1-1242", 
    "Abundance" = "Abundance~(ind./m^2)",
    "Biomass" = "Biomass~(g/m^2)",
    "TOU" = "Oxygeon~utilization")

ab_drm_plot <- 
  ggplot(ab_spatial, aes(x = DRM, y = Values)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(Variables~Cruise, scales = "free", labeller = as_labeller(ab_fn, label_parsed)) +
  xlab(Distance~to~river~mouth~(km)) +
  theme_bw()

ggsave("figure/abundance_biomass_DRM_plot.png",
       plot = ab_drm_plot)

ab_porosity_plot <-
  ggplot(ab_spatial, aes(x = Porosity, y = Values)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = Cruise)) +
  facet_grid(Variables~., scales = "free", labeller = as_labeller(ab_fn, label_parsed)) +
  xlab(Porosity~('%')) +
  scale_color_manual(values = cruise_color) +
  theme_bw()

ggsave("figure/abundance_biomass_porosity_plot.png",
       plot = ab_porosity_plot)

# spatial-temporal variation of abundance and biomass
lm(Values~DRM * Cruise, data = ab_spatial[ab_spatial$Variables == "Abundance",]) %>% summary
lm(Values~DRM * Cruise, data = ab_spatial[ab_spatial$Variables == "Abundance",]) %>% plot

lm(Values~DRM * Cruise, data = ab_spatial[ab_spatial$Variables == "Biomass",]) %>% summary
lm(Values~DRM * Cruise, data = ab_spatial[ab_spatial$Variables == "Biomass",]) %>% plot

# model averaging ----
# set up data
load("data/env_sel.RData")
ab_env <- 
  ab %>% 
  full_join(env_sel_scale) %>% 
  pivot_longer(cols = c("Abundance", "Biomass"),
               names_to = "Variables",
               values_to = "Values") %>% 
  select(-all_of(c("Date", "Latitude")))

# set up formula
ab_formula <- as.formula(paste("Values ~", paste0(env_sel, collapse = " + ")))

# set up density global model
den_env <- ab_env[ab_env$Variables == "Abundance",]
den_globalmodel <- lm(ab_formula, data = den_env, na.action = "na.fail")

# dredge
den_dredge <- dredge(den_globalmodel, rank = "AICc")

# model.avg 95% confint 
den_modavg <- 
  model.avg(den_dredge, 
            revised.var = T)
summary(den_modavg)
plot(den_modavg)


# set up biomass global model
bio_env <- ab_env[ab_env$Variables == "Biomass" & ab_env$Cruise == "OR1-1219",]
bio_globalmodel <- lm(ab_formula, data = bio_env, na.action = "na.fail")

# dredge
bio_dredge <- dredge(bio_globalmodel, rank = "AICc")

# model.avg 95% confint 
bio_modavg <- 
  model.avg(bio_dredge, 
            revised.var = T)
summary(bio_modavg)
plot(bio_modavg)



# testing: OU modelling ----
ou_env <- ou %>% left_join(env_sel_scale)
ou_formula <- as.formula(paste0("BOU ~ ", paste(env_sel[-1], collapse = " + ")))
ou_globalmodel <- lm(ou_formula, data = ou_env, na.action = "na.fail")
ou_dredge <- dredge(ou_globalmodel, rank = "AICc")
ou_modavg <- model.avg(ou_dredge, revised.var = T)
summary(ou_modavg)
plot(ou_modavg)


# testing: ou relationship with ab
ab$Deployment <- ab$Deployment %>% as.character()
ou$Tube <- ou$Tube %>% as.double()
ou_ab <- ou %>% left_join(ab) %>% pivot_longer(cols = c("Abundance", "Biomass"), names_to = "Variables",  values_to = "Values", )
ggplot(ou_ab, aes(x = Values, y = BOU, color = Cruise)) +
  geom_smooth(method = "lm") +
  geom_point()+
  facet_wrap(~Variables, scales = "free_x") +
  scale_color_manual(values = cruise_color) +
  theme_bw()

# fit lm
lm(BOU~Values * Cruise, data = ou_ab[ou_ab$Variables == "Abundance",]) %>% summary
lm(BOU~Values * Cruise, data = ou_ab[ou_ab$Variables == "Biomass",]) %>% summary

################
# Figure output
################
ggsave("figure/standingstock_vs_DRM.png", plot = standingstock_vs_DRM)
ggsave("figure/standingstock_vs_Station.png", plot = standingstock_vs_Station)
