#################################
# Mapping the Gaoping River-shelf
#################################
# Description
# Plot the map of the nSCS and the Gaoping continental shelf

# Author: Yen-Ting Chen
# Date of creation: 2023/09/05
# Date of last modification: 2023/09/05

#######################
# Set up environment
#######################
# clean environment
rm(list = ls())

# Load packages
library(dplyr) # data manipulation
library(tidyr)
library(ggplot2) # visualization
library(GGally)
library(ggrepel)
library(patchwork)
library(pals)
library(ggsn)
library(GRSmacrofauna)  # data package

# load object
load("data/cruise_color.Rdata")

# extract station lat and long
st <- 
  env %>% 
  select(Station, Latitude, Longitude) %>% 
  group_by(Station) %>% 
  summarise(Latitude = mean(Latitude),
            Longitude = mean(Longitude))

# subset map ----
bathy_map_sub <- 
  bathy_map %>% 
  filter(Longitude > 118.5 & Longitude < 123) %>% 
  filter(Latitude > 21 & Latitude < 26.5)

grs_map_sub <- 
  grs_map %>% 
  filter(Longitude > 120.15 & Longitude < 120.7) %>% 
  filter(Latitude > 22 & Latitude < 22.8)

############
# bathy map 
############
p_fill <- 
  ggplot()+
  # sea
  geom_raster(data = bathy_map_sub[bathy_map_sub$Elevation < 0,], 
              aes(x = Longitude, y = Latitude, fill = Elevation)) +
  #land
  geom_raster(data = bathy_map_sub[bathy_map_sub$Elevation >= 0,], 
              aes(x = Longitude, y = Latitude), fill = "black") +
  # color
  scale_fill_gradientn(colors = kovesi.linear_blue_5_95_c73(100),
                       breaks = c(0, 1:8*(-1000)),
                       labels = abs)+
  xlab(Longitude~(degree*E))+
  ylab(Latitude~(degree*N))+
  coord_fixed(expand = FALSE) +
  labs(fill = "Depth (m)")+
  theme_bw()

# setting up objects
tw <- data.frame(Lo = 120.9, La = 23.6)
ts <- data.frame(Lo = 119.4, La = 24.3)
po <- data.frame(Lo = 122.3, La = 23)
ss <- data.frame(Lo = 119.3, La = 21.5)

# 
p_fill <-
  p_fill +
  # study area
  annotate(geom = "rect", xmin = 120.15, xmax = 120.7, ymin = 22, ymax = 22.8, color = "red", alpha = 0, size = 1.5) +
  # add location names
  annotate(geom = "text", x = tw$Lo, y = tw$La, color = "white", label = "Taiwan") +
  annotate(geom = "text", x = ts$Lo, y = ts$La, color = "black", label = "Taiwan Strait") +
  annotate(geom = "text", x = po$Lo, y = po$La, color = "white", label = "Pacific Ocean") +
  annotate(geom = "text", x = ss$Lo, y = ss$La, color = "black", label = "South China Sea") 

#########
# grs map
#########
# set up bathymetric color
p1_fill <- 
  ggplot()+
  # land
  geom_tile(data = grs_map_sub[grs_map_sub$Elevation >= 0,], 
            aes(x = Longitude, y = Latitude), fill = "black") +
  # sea
  geom_tile(data = grs_map_sub[grs_map_sub$Elevation < 0,],
            aes(x = Longitude, y = Latitude, fill = Elevation)) +
  # 200m-isobath
  geom_contour(data = grs_map_sub,
               aes(x = Longitude, y = Latitude, z = Elevation),
               breaks = seq(-200, -1400, -200),
               size = 0.3,
               color = "black")+
  # color
  scale_fill_gradientn(colors = kovesi.linear_blue_5_95_c73(100),
                       breaks = seq(0,-1400, -200),
                       labels = abs,
                       limits = c(-1400,0))+
  xlab(Longitude~(degree*E))+
  ylab(Latitude~(degree*N))+
  coord_fixed(expand = FALSE) +
  labs(fill = "Depth (m)")+
  theme_bw()

# add station
p1_fill <- 
  p1_fill +
  geom_point(data = st, 
             aes(x = Longitude, 
                 y = Latitude),
             shape = 16,
             stroke = 1.1) +
  geom_text_repel(data = st, 
                  aes(x = Longitude, 
                      y = Latitude,
                      label = Station),
                  seed = 3)
# add names
gr <- data.frame(Lo = 120.423960, La = 22.470504)
gc <- data.frame(Lo = 120.22, La = 22.3)
fc <- data.frame(Lo = 120.5, La = 22.18)

p1_fill <-
  p1_fill +
  # river mouth
  annotate(geom = "point", x = gr$Lo, y = gr$La, color = 'white', stroke = 1.2, shape = 1, size = 3) +
  annotate(geom = "segment", x = gr$Lo, y = gr$La, color = 'white',
           xend = gr$Lo + 0.05, yend = gr$La + 0.07, size = 1) +
  annotate(geom = "text", x = gr$Lo + 0.05, y = gr$La + 0.11, color = 'white', label = "Gaoping\nRiver mouth") +
  # canyons
  annotate(geom = "text", x = gc$Lo, y = gc$La, color = "white", label = "Gaoping\nCanyon") +  
  annotate(geom = "text", x = fc$Lo, y = fc$La, color = "white", label = "Fangliao\nCanyon")  

##############
# Store output
##############
map <- 
  (p_fill + theme(legend.position = c(0.99,  0.99), legend.justification=c(1, 1))) +
  (p1_fill + theme(legend.position = c(0.99, 0.99), legend.justification=c(1, 1))) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")

ggsave("figure/publish/map.png",
       plot = map,
       width = 12,
       height = 6)
