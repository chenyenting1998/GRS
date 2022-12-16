###############################################
###### Macrofauna composition - Hellinger #####
###############################################

den_dist <- 
  den_wide[-c(1:4)] %>% 
  decostand(method = "hellinger")


den_rda <- rda(den_dist)
screeplot(den_rda)

# get scores
# attach stations 
den_rda_sites <- 
  scores(den_rda, scaling = "sites")$sites %>% 
  cbind(env[c("Cruise", "Station")]) %>% 
  as.data.frame()
den_rda_species <- 
  scores(den_rda, scaling = "sites")$species %>% 
  as.data.frame()

# extract eigenvalue
den_rda_eig <- round((den_rda$CA$eig/sum(den_rda$CA$eig))*100, 2)

# plot screeplot
screeplot(den_rda)

# plot PCA
den_rda_plot <- 
  ggplot()+
  # plot 
  stat_ellipse(data = den_rda_sites, 
               aes(x = PC1, y = PC2, color = Cruise, fill = Cruise), 
               type = "norm", geom = "polygon",
               size = 1.5,
               alpha = 0.15,
               level = 0.95) +
  
  # # plot env variables
  # geom_segment(data = den_rda_species,
  #              aes(x = 0, y = 0, xend = PC1, yend = PC2),
  #              size = .8, color = "blue")+
  # geom_label(data = den_rda_species, 
  #            aes(x = PC1, y = PC2, label = rownames(den_rda_species))) +
  # plot stations
  geom_label(data = den_rda_sites, 
             aes(x = PC1, y = PC2, color = Cruise, label = Station)) +
  # change axis label
  xlab(paste0("PC1 (", env_pca_eig[1], "% of total variance explained)")) +
  ylab(paste0("PC2 (", env_pca_eig[2], "% of total variance explained)")) +
  scale_color_manual(values = cruise_color) +
  scale_fill_manual(values = cruise_color) +
  coord_fixed() +
  theme_bw()
