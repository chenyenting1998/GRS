# produce screeplot
plot_scree <-
  function(rda_object){
    percent <- round(eigenvals(rda_object)/ sum(eigenvals(rda_object)) * 100, 2) %>% as.vector
    data <- 
      data.frame(PC = factor(names(eigenvals(rda_object)), names(eigenvals(rda_object))),
                 percent = percent)
    ggplot(data, aes(x = PC, y = percent)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = mean(percent), color = "red", linetype = 2) +
      xlab("Principal component") +
      ylab("Variance explained (%)")+
      theme_bw()
  }
