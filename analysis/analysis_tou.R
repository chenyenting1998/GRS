######################################
##### Sediment oxygen consumption ####
######################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(GGally)
library(GRSmacrofauna)

# BOU
ou_stat <- 
  ou %>% 
  group_by(Cruise, Station) %>% 
  summarize(BOU_mean = -mean(BOU),
            BOU_sd = sd(BOU)) %>% 
  full_join(env)

ggplot(ou_stat, aes(x = Station, y = BOU_mean))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = BOU_mean - BOU_sd, ymax = BOU_mean + BOU_sd),
                width = 0.5) +
  ylab(Benthic-meidated~oxygen~utilization~(mmol/m^2/d)) +
  facet_grid(~Cruise, scales = "free") +
  theme_bw()

# OU vs. DRM
ou_env <- 
  ou %>% 
  full_join(env, by = c("Cruise", "Station")) 

ggplot(ou_env, aes(x = DRM, y = BOU))+
  geom_point()+
  geom_smooth(method = "lm") +
  xlab(Distance~to~river~mouth~(km))+
  ylab(Benthic-meidated~oxygen~utilization~(mmol/m^2/d)) +
  facet_grid(~Cruise, scales = "free") +
  theme_bw()

# how does macrofauna composition affects sediment oxygen utilization?
# how does sedimentary environment affects sediment oxygen utilization?
# Which is more important in mediating sedimetn OU, macrofauna or environment?

# sediment macrofauna assemblage is shaped by organic matter quality, quantitiy and many more factors

