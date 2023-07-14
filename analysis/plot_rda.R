# get rda_output
get_rda_output <- 
  function(rda_result,
           community_matrix_metadata,
           env_abbr,
           goodness_threshold = 0.4){
    # extract site scores
    rda_sites <-
      scores(rda_result, scaling = 1)$sites %>%
      as.data.frame() %>%
      cbind(community_matrix_metadata)
    # extract env scores
    rda_env <-
      scores(rda_result, scaling = 1)$biplot %>%
      as.data.frame
    rda_env$abbr <- env_abbr[match(rownames(rda_env), names(env_abbr))]
    # extract species scores
    rda_species <-
      scores(rda_result, scaling = 2)$species %>%
      as.data.frame()
    # species goodness
    rda_goodness <-
      goodness(rda_result,
               choices = 1:2,
               display = "species",
               model = "CCA") %>%
      as.data.frame()
    # attach goodness to species
    rda_species$Goodness <- rda_goodness$RDA2[match(rownames(rda_goodness), rownames(rda_species))]
    rda_species$Show <- rda_species$Goodness> goodness_threshold
    
    # combine results
    results <- list(rda_sites = rda_sites, 
                    rda_env = rda_env,
                    rda_species = rda_species)
    return(results)
  }
# plot RDA biplot scaling = 1
plot_rda_sc1 <- 
  function(rda_sites,
           rda_env,
           rda_result){
    # internal function
    rda_eig_percent<- function(x){
      # extract rda axis
      rda_axis <- grep("RDA", names(eigenvals(x)))
      rda_eig <- eigenvals(x)[rda_axis]
      as.vector(round((rda_eig / sum(rda_eig) * 100), 2))
    }
    # plot
    rda_plot <- 
      ggplot() +
      stat_ellipse(data = rda_sites,
                   aes(x = RDA1, y = RDA2, color = Cruise, fill = Cruise),
                   geom = "polygon",
                   type = "norm",
                   size = 1.5,
                   level = .95,
                   alpha = .1) +
      # plot stations
      geom_point(data = rda_sites,
                 aes(x = RDA1, y = RDA2, color = Cruise)) +
      geom_text_repel(data = rda_sites,
                      aes(x = RDA1, y = RDA2, label = Station, color = Cruise),
                      seed = 1) +
      # plot species
      # geom_label(data = den_rda_species,
      #            aes(x = RDA1, y = RDA2, label = rownames(den_rda_species))) +
      # plot env
      geom_segment(data = rda_env,
                   aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                   arrow = arrow(),
                   size = .3, color = "black")+
      geom_text(data = rda_env,
                aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = abbr),
                parse = TRUE) +
      scale_color_manual(values = cruise_color) +
      scale_fill_manual(values = cruise_color) +
      annotate(geom = "text", 
               x = Inf,
               y = Inf,
               hjust = 1.01,
               vjust = 1.1,
               label = paste0("Total explained variance: ", 
                              round(RsquareAdj(rda_result)$r.squared* 100, 2),
                              "%")) +
      xlab(paste0("RDA1 (", rda_eig_percent(rda_result)[1], "% of the total explained variance)")) +
      ylab(paste0("RDA2 (", rda_eig_percent(rda_result)[2], "% of the total explained variance)")) +
      coord_fixed() +
      theme_bw()
    return(rda_plot)
  }