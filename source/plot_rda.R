# get rda_output
get_rda_output <- 
  function(rda_result,
           community_matrix_metadata,
           env_abbr,
           scaling = 1,
           goodness_threshold = 0.4){
    
    # extract site scores
    rda_sites <-
      scores(rda_result, scaling = scaling)$sites %>%
      as.data.frame() %>%
      cbind(community_matrix_metadata)
    
    # extract env scores
    rda_env <-
      scores(rda_result, scaling = scaling)$biplot %>%
      as.data.frame()
    rda_env$abbr <- env_abbr[match(rownames(rda_env), names(env_abbr))]
    
    # extract species scores
    rda_species <-
      scores(rda_result, scaling = scaling)$species %>%
      as.data.frame()
    rda_species$Taxon <- rownames(rda_species)
    
    # species goodness
    rda_goodness <-
      goodness(rda_result,
               model = "CCA") %>%
      as.data.frame()
    
    # attach goodness to species
    rda_species$Goodness <- rda_goodness[length(rda_goodness)]
    rda_species$Show <- rda_species$Goodness> goodness_threshold
    
    # combine results
    results <- list(rda_sites = rda_sites, 
                    rda_env = rda_env,
                    rda_species = rda_species)
    return(results)
  }
rda_eig_percent<- function(x){
  # extract rda axis; returns a numeric vector
  rda_axis <- grep("RDA", names(eigenvals(x)))
  # extract eig vals
  rda_eig <- eigenvals(x)[rda_axis]
  # get r.squared in %
  eig_percent <- eigenvals(x)[rda_axis] / sum(eigenvals(x)) * 100
  # round them as percent
  eig_percent <- as.vector(round(eig_percent, 2))
  return(eig_percent)
}
# plot RDA triplot
plot_rda <- 
  function(rda_sites,
           rda_env,
           rda_species,
           rda_result,
           stretch = 1,
           scaling = 1){
    # extract r2 and adj.R2
    # R2 <- round(RsquareAdj(rda_result)$r.squared* 100, 2)
    # adjR2 <- round(RsquareAdj(rda_result)$adj.r.squared* 100, 2)
    
    if(scaling == 1){
      # plot
      rda_plot <- 
        ggplot() +
        # add v and hline
        geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
        geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
        
        # # plot env
        geom_segment(data = rda_env,
                     aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                     arrow = arrow(angle = 0),
                     size = .5, color = "black")+
        geom_text(data = rda_env,
                   aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = abbr),
                   size = 3,
                   parse = TRUE) +
        # species vectors
        geom_segment(data = rda_species[rda_species$Show == TRUE,],
                     aes(x = 0, y = 0, xend = RDA1 * stretch, yend = RDA2 * stretch),
                     arrow = arrow(angle = 0),
                     size = .5, 
                     color = "purple")+
        geom_label(data = rda_species[rda_species$Show == TRUE,], 
                   aes(x = RDA1 * stretch, y = RDA2 * stretch, label = Taxon),
                   size = 3,
                   color = "purple",
                   label.padding = unit(0.15, "lines")) +   
        scale_color_manual(values = cruise_color) +
        # plot stations
        geom_text(data = rda_sites,
                  aes(x = RDA1,
                      y = RDA2,
                      label = Station,
                      color = Month),
                  size = 3) +

        # # plot species
        scale_color_manual(values = cruise_color) +
        # add R2
        xlab(paste0("RDA1 (", rda_eig_percent(rda_result)[1], "% of the total variance)")) +
        ylab(paste0("RDA2 (", rda_eig_percent(rda_result)[2], "% of the total variance)")) +
        coord_fixed() +
        theme_bw()
    }
    if(scaling == 2){
      # plot
      rda_plot <- 
        ggplot() +
        # add v and hline
        geom_vline(xintercept = 0, color = "black", size = .5, linetype = 2) +
        geom_hline(yintercept = 0, color = "black", size = .5, linetype = 2) +
        
        # env vectors
        geom_segment(data = rda_env,
                     aes(x = 0, 
                         y = 0, 
                         xend = RDA1 * stretch, 
                         yend = RDA2 * stretch),
                     arrow = arrow(angle = 0),
                     size = .5, 
                     color = "black")+
        geom_label(data = rda_env,
                   aes(x = RDA1 * stretch, 
                       y = RDA2 * stretch, 
                       label = abbr),
                   color = "black",
                   size = 3,
                   parse = TRUE) +
        # species vectors
        geom_segment(data = rda_species[rda_species$Show == TRUE,],
                     aes(x = 0, 
                         y = 0, 
                         xend = RDA1 * stretch, 
                         yend = RDA2 * stretch),
                     arrow = arrow(angle = 0),
                     size = .5, 
                     color = "purple")+
        geom_label(data = rda_species[rda_species$Show == TRUE,], 
                   aes(x = RDA1 * stretch,
                       y = RDA2 * stretch,
                       label = Taxon),
                   size = 3,
                   color = "purple",
                   label.padding = unit(0.15, "lines")) +
        # plot stations
        geom_text(data = rda_sites[rda_sites$Station %in% c("S1", "S2", "S4"),],
                  aes(x = RDA1, 
                      y = RDA2, 
                      label = Station, 
                      color = Month),
                  size = 3) +
        geom_label(data = rda_sites[rda_sites$Station %in% c("S3", "S5", "S6", "S7"),],
                  aes(x = RDA1, 
                      y = RDA2, 
                      label = Station, 
                      color = Month),
                  size = 3,
                  fontface = "bold",
                  label.padding = unit(0.15, "lines")) +
        scale_color_manual(values = cruise_color) +
        scale_fill_manual(values = cruise_color) +
        # annotate(geom = "text", 
        #          x = Inf,
        #          y = Inf,
        #          hjust = 1.01,
        #          vjust = 1.5,
        #          label = paste0("R squared: ", R2, "%")) +
        # annotate(geom = "text", 
        #          x = Inf,
        #          y = Inf,
        #          hjust = 1.01,
        #          vjust = 2.8,
        #          label = paste0("Adjusted R squared: ", adjR2, "%")) +
        xlab(paste0("RDA1 (", rda_eig_percent(rda_result)[1], "% of the total variance)")) +
        ylab(paste0("RDA2 (", rda_eig_percent(rda_result)[2], "% of the total variance)")) +
        coord_fixed() +
        theme_bw()
      
    }
      return(rda_plot)
  }
