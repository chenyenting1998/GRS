##############################################
# Macrofauna abundance, biomass, and diveristy 
##############################################
# Description
# Create plots and tables for macrofauna density and biomass 

# Author: Yen-Ting Chen
# Date of creation: Unknown
# Date of last modification: 2023/07/06

#######################
# Set up environment
#######################
# clean environment
rm(list = ls())

# Load packages
library(dplyr)
library(tidyr)
library(vegan)
library(Hmisc)
library(MuMIn)
library(corrplot)
library(ggplot2)
library(GGally)
library(ggrepel)
library(patchwork)
library(GRSmacrofauna)
library(writexl)

# Load data
load("data/macrofauna composition.Rdata")
load("data/wide_data.Rdata")
load("data/cruise_color.Rdata") 
load("data/env_variables.Rdata")
load("data/env_selected.Rdata")
load("data/ou.Rdata")

# define function
get_model.avg_results <- function(sum_object, path){

  # extract objects
  coefficients <- as.data.frame(sum_object$coefficients)
  msTable <- as.data.frame(sum_object$msTable)
  sw <- as.data.frame(sum_object$sw)
  coefmat.full <- as.data.frame(sum_object$coefmat.full)
  coefmat.subset <- as.data.frame(sum_object$coefmat.subset)
  
  # add rownames
  coefficients <- cbind(" " = rownames(coefficients), coefficients)
  msTable <- cbind(" " = rownames(msTable), msTable)
  sw <- cbind(" " = rownames(sw), sw)
  coefmat.full <- cbind(" " = rownames(coefmat.full), coefmat.full)
  coefmat.subset <- cbind(" " = rownames(coefmat.subset), coefmat.subset)
  
  # write xlsx
  write_xlsx(list(coefficients = coefficients, 
                  msTable = msTable, 
                  sw = sw, 
                  coefmat.full = coefmat.full, 
                  coefmat.subset = coefmat.subset),
             path = path)
}

# define multiple corer area; diameter = 10 cm
mc_area <- (0.1 / 2) ^ 2 * pi

###############################################
# 1. Calculate unit area abundance and biomass
###############################################
ss <-
  composition %>%
  group_by(Cruise, Station, Deployment, Tube) %>%
  summarise(
    Abundance = sum(Count) / mc_area,
    Biomass = sum(Biomass) / 1000 / mc_area) %>%  # mg to g
  ungroup()

# output .xlsx table
ss_table <-
  ss %>% 
  group_by(Cruise, Station) %>% 
  summarise(Abundance_mean = mean(Abundance),
            Abundance_se = sd(Abundance) / (3^0.5),
            Biomass_mean = mean(Biomass),
            Biomass_se = sd(Biomass) / (3^0.5)) 

write_xlsx(ss_table, "table/abundance_biomass.xlsx")

############################################################
# 2. Plot standing stock relationships with station and DRM
############################################################
# set up data.frame
ss_plot_data <-
  ss %>% 
  pivot_longer(cols = c("Abundance", "Biomass"),
              names_to = "Variable",
              values_to = "Value") %>% 
  group_by(Cruise, Station, Variable) %>% 
  summarise(mean = mean(Value),
            sd = sd(Value)) %>% 
  left_join(env)
  
plot_ss <- function(object, x, y, sd, text_repel = FALSE){
  plot <- 
    object %>% 
    ggplot(aes(x = .data[[x]], 
               y = .data[[y]], 
               ymin = .data[[y]] - .data[[sd]], 
               ymax = .data[[y]] + .data[[sd]],
               color = Cruise)) +
    geom_pointrange(position = position_dodge(width = .5))
  
    if(text_repel == TRUE){
    plot <-
      plot + 
      geom_text_repel(data = object,
                      aes(x = .data[[x]],
                          y = .data[[y]],
                          label = Station,
                          color = Cruise),
                      seed = 10)
    }
  plot <- 
    plot +
    scale_color_manual(values = cruise_color) +
    facet_wrap(~Variable, scales = "free") +
    ylab("Value") +
    theme_bw()
  
  return(plot)
}

# vs station (discret var.)
ss_vs_Station <- plot_ss(ss_plot_data, "Station", "mean", "sd")
ss_vs_DRM <- plot_ss(ss_plot_data, "DRM", "mean", "sd", TRUE) + xlab("Distance to river mouth (km)")
ss_vs_depth <- plot_ss(ss_plot_data, "Depth", "mean", "sd", TRUE) + xlab("Depth (m)")

ggsave("figure/standing_stock_station.png",
       plot = ss_vs_Station,
       scale = 1.2,
       width = 6,
       height = 3)

ggsave("figure/standing_stock_drm_depth.png",
       plot = ss_vs_depth/ss_vs_DRM,
       scale = 1.2,
       width = 8,
       height = 6)

#######################################################
# 3. Statistical relationships between SS and DRM/Depth
#######################################################
# extract spatiotemporal data.frame
sp <- 
  env %>% 
  select(Cruise, Station, Depth, DRM) 
# scale Depth and DRM
sp$Depth <- (sp$Depth - mean(sp$Depth)) / sd(sp$Depth)
sp$DRM <- (sp$DRM - mean(sp$DRM)) / sd(sp$DRM)

ss_sp <- left_join(ss, sp)

# generate a fullmodel regarding spatiotemporal relationships
# Develop hypotheses for variable inputs:
# Single variables:
#    DRM    : Areas closer to the river is prone to sedimentation
#    Depth  : deeper waters can act as a refuge for the benthos
#    Cruise : before and after flood season + a typhoon

# Two-way Interaction:
#    DRM:Depth    : further/deeper waters seem to be less prone to sedimentation
#    Cruise:Depth : shallower waters might be more abu. before flood season.
#    Cruise:DRM   : closer waters might be more abu. before flood season.

# Three-way interaction:
#    Cruise:DRM:Depth : The interaction between water depth and distance varies at flood season.

#############
# Abundance #
#############
abu_sp_fullmodel <- 
  lm(log10(Abundance) ~
       DRM + 
       Depth + 
       Cruise +
       DRM:Depth +
       DRM:Cruise +
       Depth:Cruise +
       Depth:DRM:Cruise,
     data = ss_sp,
     na.action = "na.fail")

abu_sp_full_sum <- summary(abu_sp_fullmodel)
abu_sp_full_sum
# there is some heteroscedasticity
# Imma just pretend this does not matter
plot(abu_sp_fullmodel)

# dredge 
abu_sp_d <- dredge(abu_sp_fullmodel)
nrow(abu_sp_d) # 19 models in total
View(abu_sp_d)
head(abu_sp_d)

# extract models
abu_sp_d_subst <- get.models(abu_sp_d, subset = delta<6)
# model averaging
abu_sp_avg <- model.avg(abu_sp_d_subst)
# store averaging result
get_model.avg_results(abu_sp_avg, "table/abu_sp_model.avg.xlsx")

###########
# Biomass #
###########
bio_sp_fullmodel <- 
  lm(log10(Biomass) ~ 
       # single variable
       DRM + Depth + Cruise +
       # two-variable interaction
       DRM:Depth + DRM:Cruise + Depth:Cruise +
       # three variable interaction
       Depth:DRM:Cruise, 
     data = ss_sp,
     na.action = "na.fail")

bio_sp_sum <- summary(bio_sp_fullmodel)
bio_sp_sum
plot(bio_sp_fullmodel) # good residual patterns

# dredge 
bio_sp_d <- dredge(bio_sp_fullmodel)
View(bio_sp_d)
head(bio_sp_d) # the first model has the best AICc

# extract model with delta < 6
bio_sp_subset <- get.models(bio_sp_d, subset = delta < 6)
bio_sp_avg <- model.avg(bio_sp_subset)
get_model.avg_results(bio_sp_avg_sum, "table/bio_sp_model.avg.xlsx")


#######################
# 4. Correlation matrix
#######################
# using station-level resolution doing correlation
ss_ou_env <- 
  ss_table %>% 
  select(Cruise, Station, Abundance_mean, Biomass_mean) %>% 
  mutate(Log10Abundance = log10(Abundance_mean),
         Log10Biomass = log10(Biomass_mean)) %>% 
  left_join(ou_station[c("Cruise", 
                         "Station", 
                         "In_situ_TOU_mean",
                         "In_situ_DOU_mean",
                         "OPD_mean")]) %>% 
  mutate(In_situ_BOU = In_situ_TOU_mean - In_situ_DOU_mean) %>% 
  left_join(env[c("Cruise", "Station",env_variables)])

ss_ou_env_cor <- 
  ss_ou_env[-(1:2)] %>% 
  as.matrix %>% 
  rcorr()

# plot corrplot
png(filename = "figure/ss_ou_env_cor.png",
    width = 1600,
    height = 800)

corrplot(ss_ou_env_cor$r,
         # corrplot style
         method = "ellipse",
         type = "upper",
         diag = FALSE,
         addCoef.col = TRUE,
         # significance
         p.mat = ss_ou_env_cor$P,
         insig = "blank")
dev.off()

#########################################################
# 4. Predictors of Abundance, biomass, and o2 utilization
#########################################################
# scaling selected env variables
env_scaled <- env_selected
env_scaled[env_variables_selected] <- scale(env_scaled[env_variables_selected])
# matching ou, ss, env at core-level resolution
ss_ou_env_core <- 
  ss %>%
  # normalize data
  mutate(Log10Abundance = log10(Abundance),
         Log10Biomass = log10(Biomass)) %>% 
  select(-Abundance, - Biomass) %>% 
  left_join(ou_core[c("Cruise",
                      "Station",
                      "Deployment",
                      "Tube",
                      "In_situ_TOU",
                      "In_situ_DOU_mean",
                      "OPD_mean")], 
            by = c("Cruise", "Station", "Deployment", "Tube")) %>% 
  left_join(env_scaled) %>% 
  # remove excess metadata
  select(-all_of(c("Habitat", "Latitude", "Longitude", "Date")))

# scatter plot
ggplot(ss_ou_env_core) +
  geom_point(aes(x = Log10Biomass, y = In_situ_TOU, 
                 color = Cruise)) +
  scale_color_manual(values = cruise_color) +
  theme_bw()

# pairplot
ss_ou_pairplot <- 
  ggpairs(ss_ou_env_core,
          aes(color = Cruise, fill = Cruise),
          columns = c("Log10Abundance",
                      "Log10Biomass",
                      env_variables_selected)) +
  theme_bw()
ggsave("figure/ss_ou_pairplot.png", 
       ss_ou_pairplot,
       width = 10,
       height = 8)

#
plot_ou <- function(variable){
  plot <-
    ou_station %>% 
    left_join(env, by = c("Cruise", "Station")) %>% 
    left_join(env, by = c("Cruise", "Station")) %>% 
    ggplot(aes(x = .[[variable]], 
               y = In_situ_TOU_mean,
               color = Cruise, 
               label = Station)) +
    geom_pointrange(aes(ymin = In_situ_TOU_mean - In_situ_TOU_sd,
                        ymax = In_situ_TOU_mean + In_situ_TOU_sd)) +
    geom_text_repel() +
    scale_color_manual(values = cruise_color) +
    xlab(variable) +
    theme_bw()
  return(plot)
}

plot_ou("Temperature") + facet_wrap(~Cruise)
plot_ou("Porosity")
plot_ou("DRM")

################
# MMI: Abundance
################
# abundance full model
abu_fm <- 
  lm(Log10Abundance ~
       Temperature +
       Fluorescence +
       D50 + 
       TOC +
       CN +
       Chla +
       Porosity,
     data = ss_ou_env_core,
     na.action = "na.fail")

# full model statistics
summary(abu_fm)
# full model residual variations
plot(abu_fm)

# dredge all combinations of parameter inputs
abu_d <- dredge(abu_fm)
# View model set
View(abu_d)
# number of models
nrow(abu_d)

# extract top models 
abu_topmodels <- get.models(abu_d, subset = delta <4)

# model average
abu_avg <- model.avg(abu_topmodels)
abu_sum <- summary(abu_avg)

get_model.avg_results(abu_sum, "table/abu_model.avg.xlsx")

coefTable(abu_avg)
plot(abu_avg)

##############
# MMI: Biomass
##############
# Biomass full model
bio_fm <- 
  lm(Log10Biomass ~
       Temperature +
       Fluorescence +
       D50 + 
       TOC +
       CN +
       Chla +
       Porosity,
     data = ss_ou_env_core,
     na.action = "na.fail")

# full model statistics
summary(bio_fm)
# full model residual variations
plot(bio_fm)

# dredge all combinations of parameter inputs
bio_d <- dredge(bio_fm)
# View model set
View(bio_d)
# number of models
nrow(bio_d)

# extract top models
bio_topmodels <- get.models(bio_d, subset = delta <4)

# model average
bio_avg <- model.avg(bio_topmodels)
bio_sum <- summary(bio_avg)

get_model.avg_results(bio_sum, "table/bio_model.avg.xlsx")

#######################################
# MMI: In Situ Total Oxygen Utilization
#######################################
# in situ TOU full model
TOU_fullmodel <-
  lm(In_situ_TOU ~
       # scale(Log10Abundance) +
       # scale(Log10Biomass) +
       # Temperature +
       Fluorescence +
       D50 +
       TOC +
       CN +
       Chla +
       Porosity,
     data = ss_ou_env_core,
     na.action = "na.fail")

# fullmodel evaluation
summary(TOU_fullmodel)
plot(TOU_fullmodel)

# multimodel inference
# dredge all model combinations
TOU_d <- dredge(TOU_fullmodel) 
View(TOU_d)

# extract top models
TOU_topmodels <- get.models(TOU_d, subset = delta < 4)

# average
TOU_avg <- model.avg(TOU_topmodels)
TOU_sum <- summary(TOU_avg)

get_model.avg_results(TOU_sum, "table/TOU_model.avg.xlsx")


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
ggsave("figure/standing_stock_vs_DRM.png", plot = ss_vs_DRM)
ggsave("figure/standing_stock_vs_Station.png", plot = ss_vs_Station)
