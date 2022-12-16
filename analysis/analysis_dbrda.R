library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(GGally)
library(GRSmacrofauna)
library(writexl)

# load
load("data/wide_data.Rdata")
load("data/dist_data.Rdata")
load("data/env_sel.Rdata")

# plot density dbrda ----
den_env <- 
  full_join(den_wide, env_sel_scale) %>% # match data
  select(all_of(env_sel))

# calculate den dbrda
den_dbrda <- 
  dbrda(den_dist ~ .,
        data = den_env[-c(1:4)])

# backward selection
set.seed(1)
den_dbrda_bs <- ordistep(den_dbrda, direction = "backward", permutaions = 9999)
anova(den_dbrda_bs)

# test variable and axis power
den_dbrda_terms <- anova(den_dbrda_bs, by = 'terms') # type I anova;  sequential
den_dbrda_margin <- anova(den_dbrda_bs, by = "margin") # type III anova; influence of one effect when other effects are present
den_dbrda_axis <- anova(den_dbrda_bs, by = "axis") 

# extract summary
den_dbrda_su <- summary(den_dbrda_bs)

# extract scores
den_dbrda_sites <- 
  scores(den_dbrda_bs, scaling = "sites")$sites %>% 
  as.data.frame() %>% 
  cbind(den_env[1:4])

den_dbrda_species <- 
  scores(den_dbrda_bs, scaling = "sites")$biplot %>% 
  as.data.frame()

# plot dbrda
den_dbrda_plot <-
  ggplot(data = den_dbrda_sites, 
         aes(x = dbRDA1, y = dbRDA2)) +
  stat_ellipse(aes(color = Cruise, fill = Cruise), 
               type = "norm", geom = "polygon",
               size = 1.5,
               alpha = 0.15,
               level = 0.95) +
  # station
  geom_label(aes(color = Cruise, label = Station), alpha = 0.5) +
  geom_segment(data = den_dbrda_species, 
               aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               color = "blue", 
               alpha = 0.6, 
               arrow = arrow(length = unit(0.2, "cm"))) +
  # env
  geom_label(data = den_dbrda_species, 
             aes(x = dbRDA1, y = dbRDA2, label = rownames(den_dbrda_species)),
             parse = TRUE,
             alpha = 0.5) +
  annotate(geom = "text", 
           x = Inf, y = Inf, 
           hjust = 1,
           vjust = 1.5,
           label = paste0("Total constrained: ", 
                          round(den_dbrda_su$constr.chi/den_dbrda_su$tot.chi * 100, 1), 
                          "%"))+
  xlab(paste0("dbRDA1 (", round(den_dbrda_su$cont$importance[2,1] * 100, 1),"% constrained)"))+
  ylab(paste0("dbRDA2 (", round(den_dbrda_su$cont$importance[2,2] * 100, 1),"% constrained)"))+
  # scale_color_manual(values = cruise_color) +
  # scale_fill_manual(values = cruise_color) +
  theme_bw()

ggsave("figure/den_dbrda_plot.png",
       device = "png",
       plot = den_dbrda_plot)


# plot biomass dbrda ----
bio_env <- 
  full_join(bio_wide, env_sel_scale) %>% # match data
  select(all_of(env_sel))

# calculate den dbrda
bio_dbrda <- 
  dbrda(bio_dist ~ .,
        data = bio_env[-c(1:4)])

# backward selection
set.seed(1)
bio_dbrda_bs <- ordistep(bio_dbrda, direction = "backward", permutaions = 9999)
anova(bio_dbrda_bs)

# test variable and axis power
bio_dbrda_terms <- anova(bio_dbrda_bs, by = 'terms') # type I anova;  sequential
bio_dbrda_margin <- anova(bio_dbrda_bs, by = "margin") # type III anova; influence of one effect when other effects are present
bio_dbrda_axis <- anova(bio_dbrda_bs, by = "axis") 

# extract summary
bio_dbrda_su <- summary(bio_dbrda_bs)

# extract scores
bio_dbrda_sites <- 
  scores(bio_dbrda_bs, scaling = "sites")$sites %>% 
  as.data.frame() %>% 
  cbind(bio_env[1:4])

bio_dbrda_species <- 
  scores(bio_dbrda_bs, scaling = "sites")$biplot %>% 
  as.data.frame()

# plot dbrda
bio_dbrda_plot <-
  ggplot(data = bio_dbrda_sites, 
         aes(x = dbRDA1, y = dbRDA2)) +
  stat_ellipse(aes(color = Cruise, fill = Cruise), 
               type = "norm", geom = "polygon",
               size = 1.5,
               alpha = 0.15,
               level = 0.95) +
  # station
  geom_label(aes(color = Cruise, label = Station), alpha = 0.5) +
  geom_segment(data = bio_dbrda_species, 
               aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2), 
               color = "blue", 
               alpha = 0.6, 
               arrow = arrow(length = unit(0.2, "cm"))) +
  # env
  geom_label(data = bio_dbrda_species, 
             aes(x = dbRDA1, y = dbRDA2, label = rownames(bio_dbrda_species)),
             parse = TRUE,
             alpha = 0.5) +
  annotate(geom = "text", 
           x = Inf, y = Inf, 
           hjust = 1,
           vjust = 1.5,
           label = paste0("Total constrained: ", 
                          round(bio_dbrda_su$constr.chi/bio_dbrda_su$tot.chi * 100, 1), 
                          "%"))+
  xlab(paste0("dbRDA1 (", round(bio_dbrda_su$cont$importance[2,1] * 100, 1),"% constrained)"))+
  ylab(paste0("dbRDA2 (", round(bio_dbrda_su$cont$importance[2,2] * 100, 1),"% constrained)"))+
  # scale_color_manual(values = cruise_color) +
  # scale_fill_manual(values = cruise_color) +
  theme_bw()

ggsave("figure/bio_dbrda_plot.png",
       device = "png",
       plot = bio_dbrda_plot)


