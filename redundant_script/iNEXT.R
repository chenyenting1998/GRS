# Calculating effective number of species

# You should know about
# Hill numbers
# sample completeness

# set up library
library(dplyr)
library(tidyr)
library(iNEXT)
library(ggplot2)

# set up environment 
load("data/macrofauna composition.Rdata")

#######################
# 1. set up iNEXT input
#######################
# pool data by station
count_station <-
  composition %>% 
  group_by(Cruise, Station, Taxon) %>% 
  summarize(Count = sum(Count)) %>% 
  pivot_wider(names_from = "Taxon", values_from = "Count", values_fill = 0)

# as.matrix, transpose, as.data.frame
x_station <- 
  count_station[-(1:2)] %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame()
# fill column names
colnames(x_station) <- paste0(count_station$Cruise, "_", count_station$Station)

##############
# 2. Run iNEXT
##############
iNEXT_output <- iNEXT(x_station, q = 0, datatype = "abundance", endpoint = 500)

# the $DataInfo output
iNEXT_output$DataInfo
# Sites with more rare species have higher richness
# f1: singleton; f2: doubleton; and so on
ggplot(iNEXT_output$DataInfo, aes(x = f1+f2+f3, y = S.obs)) +
  geom_point() +
  stat_smooth(method = "lm") +
  xlab("# of spp. f1, f2, and f3 summed together") +
  ylab("Observed richness") +
  theme_bw()

iNEXT_output$iNextEst

# The $AsyEst output
# Yields the estimation at the asymptote
iNEXT_output$AsyEst
ggplot(iNEXT_output$AsyEst, aes(y = Site)) +
  geom_errorbar(aes(xmin = Estimator -s.e., xmax = Estimator + s.e.)) +
  geom_point(aes(x = Estimator)) +
  facet_wrap(~Diversity, scales = "free")+
  theme_bw()

# sample-size-based R/E curve
ggiNEXT(x, type = 1)
# sample completeness curve
ggiNEXT(x, type = 2)
# Coverage-based R/E curve
ggiNEXT(x, type = 3)
