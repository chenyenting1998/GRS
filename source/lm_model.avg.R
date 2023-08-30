get_lm_summary <- function(x, path){
  # extract coefs ====
  coefficients <- as.data.frame(x$coefficients)
  coefficients <- cbind(" " = row.names(coefficients), coefficients)
  
  # extract full model statistics ====
  r.squared <- x$r.squared
  adj.r.squared <- x$adj.r.squared
  fstatistic <- x$fstatistic
  df <- x$df
  
  ## calculate p-value ====
  p_value <- pf(fstatistic[1], fstatistic[2], fstatistic[3], lower.tail = FALSE)
  
  # set up a dataframe for full model statistics ====
  fullmodel_statistcs <- 
    data.frame("r squared" = r.squared,
               "adj. r squared" = adj.r.squared,
               "f statistic" = fstatistic[1], 
               "df" = df[2], # get n-p only
               "p value" = p_value)
  write_xlsx(list(Coefficients = coefficients,
                  Fullmodel_statistics = fullmodel_statistcs),
             path = path)
}
calculate_adjR2 <- function(fullmodel, modelset){
  # set model 
  n <- nrow(fullmodel$model)
  k <- modelset$df - 2 # minus parameters of intercept and standard error
  # calculate adj.R2
  modelset$`adj.R^2` <- 1 - ((1 - modelset$`R^2`) * (n - 1) / (n - k - 1))
  return(modelset)
}
get_model.avg_results <- function(avg_model, path){
  # summary() model.avg object
  if(!class(avg_model) %in% "summary.averaging"){
    sum_object <- summary(avg_model)
  }
  
  # extract objects
  coefficients <- as.data.frame(sum_object$coefficients)
  msTable <- as.data.frame(sum_object$msTable)
  sw <- as.data.frame(sum_object$sw)
  coefmat.full <- as.data.frame(sum_object$coefmat.full)
  coefmat.subset <- as.data.frame(sum_object$coefmat.subset)
  
  # add rownames
  coefficients <- cbind(" " = rownames(coefficients), coefficients)
  msTable <- cbind(" " = rownames(msTable), msTable)
  sw <- cbind(" " = rownames(sw), sw)
  coefmat.full <- cbind(" " = rownames(coefmat.full), coefmat.full)
  coefmat.subset <- cbind(" " = rownames(coefmat.subset), coefmat.subset)
  
  # write xlsx
  write_xlsx(list(coefficients = coefficients, 
                  msTable = msTable, 
                  sw = sw, 
                  coefmat.full = coefmat.full, 
                  coefmat.subset = coefmat.subset),
             path = path)
}
