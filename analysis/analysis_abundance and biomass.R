######################################################
##### Macrofauna abundance, biomass, and drivers #####
######################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(GGally)
library(GRSmacrofauna)
library(writexl)
library(MuMIn)

# load
load("data/composition.Rdata")

# variables
sw <- 1.13 # specific weight
mc_area <- (0.1/2)^2*pi # multiple corer area; diameter = 10 cm

# plot abundance and biomass ----
ab <- 
  comp %>% 
  group_by(Cruise, Station, Deployment, Tube) %>% 
  summarize(Abundance = sum(Abundance)/mc_area,
            Biomass = sum(Biomass)/1000/mc_area)
# table
ab %>% 
  ungroup() %>% 
  group_by(Cruise, Station) %>% 
  summarize(mean_a = mean(Abundance),
            sd_a = sd(Abundance),
            mean_b = mean(Biomass),
            sd_b = sd(Biomass)) %>% 
  write_xlsx("xlsx/abundance_biomass.xlsx")

# fit abundance and biomass to DRM ----
ab_spatial <- 
  ab %>% 
  full_join(env[c("Cruise", "Station", "Depth", "DRM", "Porosity")]) %>% 
  pivot_longer(cols = c("Abundance", "Biomass"),
               names_to = "Variables",
               values_to = "Values")

ab_fn <- 
  c("OR1-1219" = "OR1-1219",
    "OR1-1242" = "OR1-1242", 
    "Abundance" = "Abundance~(ind./m^2)",
    "Biomass" = "Biomass~(g/m^2)")

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
  ggplot(ab_spatial, aes(x = Porosity, y = Values, color = Cruise)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_grid(Variables~., scales = "free", labeller = as_labeller(ab_fn, label_parsed)) +
  xlab(Porosity~('%')) +
  theme_bw()

ggsave("figure/abundance_biomass_porosity_plot.png",
       plot = ab_porosity_plot)

# spatial-temporal variation of abundance and biomass
lm(Values~DRM * Cruise, data = ab_spatial[ab_spatial$Variables == "Abundance",]) %>% summary
lm(Values~DRM * Cruise, data = ab_spatial[ab_spatial$Variables == "Biomass",]) %>% summary


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
ou_ab <- ou %>% left_join(ab) %>% pivot_longer(cols = c("Abundance", "Biomass"), names_to = "Variables",  values_to = "Values", )
ggplot(ou_ab, aes(x = Values, y = BOU, color = Cruise)) +
  geom_smooth(method = "lm") +
  geom_point()+
  facet_wrap(~Variables, scales = "free_x") +
  theme_bw()
