# extract pca output
get_pca_output <- 
  function(rda_object, 
           metadata, 
           scaling = 1, 
           goodness_threshold = 0.4){
    # get site scores
    pca_sites <-
      scores(rda_object, scaling = scaling)$sites %>%
      as.data.frame() %>% 
      cbind(metadata) 
    
    # get species scores
    pca_species <- 
      scores(rda_object, scaling = scaling)$species %>% 
      as.data.frame()
    pca_species$Taxon <- rownames(pca_species)
    
    # calculate species goodness of fit
    pca_goodness <- 
      goodness(rda_object, choices = 1:2, display = "species", model = "CA") %>% 
      as.data.frame()
    # match goodness to sc2 species scores
    pca_species$goodness <- pca_goodness$PC2[match(rownames(pca_goodness), rownames(pca_species))]
    # use 40% as a cutoff criteria
    pca_species$Show <- pca_species$goodness > goodness_threshold
    
    # extract eigenvalue
    pca_eig <- round(eigenvals(rda_object)/sum(eigenvals(rda_object))*100, 2)
    
    return(list(pca_sites = pca_sites, 
                pca_species = pca_species, 
                pca_goodness = pca_goodness, 
                pca_eig = pca_eig))
  }


# produce pca plot for assemblage ordination
plot_pca <-
  function(sites, species, eig_vector, scaling, stretch = 1){
    if(scaling == 1){
      plot <-
        ggplot()+
        
        # add v and hline
        geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
        geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
        
        # plot sp.
        geom_segment(data = species[species$Show == TRUE,],
                     aes(x = 0, 
                         y = 0,
                         xend = PC1 * stretch, 
                         yend = PC2 * stretch),
                     size = .4, 
                     color = "purple")+
        geom_label(data = species[species$Show == TRUE,],
                   aes(x = PC1 * stretch, 
                       y = PC2 * stretch, 
                       label = Taxon),
                   size = 3,
                   color = "purple") +
        # plot stations
        # geom_point(data = sites, 
        #            aes(x = PC1, 
        #                y = PC2, 
        #                color = Cruise)) +
        geom_text(data = sites,
                  aes(x = PC1,
                      y = PC2,
                      color = Cruise,
                      label = Station),
                  size = 2.5) +
        # change axis label
        xlab(paste0("PC1 (", eig_vector[1], "% of the total variance)")) +
        ylab(paste0("PC2 (", eig_vector[2], "% of the total variance)")) +
        scale_color_manual(values = cruise_color) +
        scale_fill_manual(values = cruise_color) +
        coord_fixed() +
        theme_bw()
    }
    if(scaling == 2){
      plot <- 
        ggplot()+
        # plot 
        # add v and hline
        geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
        geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
        
        # plot sp.
        geom_segment(data = species[species$Show == TRUE,],
                     aes(x = 0, 
                         y = 0,
                         xend = PC1 * stretch, 
                         yend = PC2 * stretch),
                     size = .4, 
                     color = "purple")+
        geom_label(data = species[species$Show == TRUE,],
                   aes(x = PC1 * stretch, 
                       y = PC2 * stretch, 
                       label = Taxon),
                   size = 3,
                   color = "purple") +
        
        # plot stations
        geom_text(data = sites,
                  aes(x = PC1,
                      y = PC2,
                      color = Cruise,
                      label = Station),
                  size = 2.5) +        
        # change axis label
        xlab(paste0("PC1 (", eig_vector[1], "% of the total variance)")) +
        ylab(paste0("PC2 (", eig_vector[2], "% of the total variance)")) +
        scale_color_manual(values = cruise_color) +
        scale_fill_manual(values = cruise_color) +
        coord_fixed() +
        theme_bw()
    }
    return(plot)
  }
