extract_goodness <- function(ord, model){
  goodness_result <- 
    goodness(ord, display = "species", model = model) %>% 
    as.data.frame()
  goodness_result$Taxon <- rownames(goodness_result)
  return(goodness_result)
}
