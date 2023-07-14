# Box-Cox+Chord transformation
box.cox.chord <-
  function(mat, bc.exp=0)
    {
    # Internal function
    vec.norm <- function(vec) sqrt(sum(vec^2))
    #
    chck <- apply(mat, 1, sum)
    if(any(chck == 0)) stop("Rows",which(chck==0)," of the data matrix sum to 0")
    # 
    # Apply the user-selected Box-Cox exponent (bc.exp) to the frequency data
    if(bc.exp==0) {
      tmp <- log(mat+1)
      } else {
      tmp <- mat^bc.exp
      } 
    row.norms <- apply(tmp, 1, vec.norm)
    #
    # Apply the chord transformation to matrix "tmp" before returning it
    res <- sweep(tmp, 1, row.norms, "/")
  }

    
    

# Box-Cox-Dagnelie
BCD <-
  function(mat,
           bc.exp=c(0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
           chord=TRUE)
    {
    # Internal function
    vec.norm <- function(vec) sqrt(sum(vec^2))
    #
    require(ade4)
    epsilon <- sqrt(.Machine$double.eps)
    mat <- as.matrix(mat)
    n <- nrow(mat)
    p <- ncol(mat)
    n.exp <- length(bc.exp)
    #
    if(chord) {
      res <- matrix(NA,n.exp,5)
      colnames(res) <- c("BC.exp","BC_W","BC_p-val","BC.chord_W","BC.chord_p-val")
      } else {
      res <- matrix(NA,n.exp,3)
      colnames(res) <- c("BC.exp","BC_W","BC_p-val")
      }
    res[,1] <- bc.exp
    #
    if(any(mat < 0)) stop("Negative values not allowed in community data",
                          call.=FALSE)
    chck1 <- apply(mat, 1, sum)
    if(any(chck1 == 0)) stop("One or several rows of 'mat' sum to 0", 
                             call.=FALSE)
    chck2 <- apply(mat, 2, var)
    keep.spec <- which(chck2 > epsilon)
    if(length(keep.spec) < p) {
      cat(length(keep.spec),"Species have variances > 0 and were kept\n") 
      cat("Species",which(chck2 <= epsilon),"were excluded\n")
      mat2 <- mat[,keep.spec]
      } else { mat2 <- mat } 
    #
    for(k in 1:n.exp) {
      if(bc.exp[k]==0) { 
        # If BC exponent = 0, compute log(x+1)
        # Add 1 to the data before log transformation
        tmp <- log(mat2 + 1) 
        # Add 1 to the data before applying a negative exponent 
        } else if(bc.exp[k] < 0) { tmp <- (mat2 + 1)^bc.exp[k]
        # No transformation when bc.exp= 1 
        } else if(bc.exp[k] == 1){ tmp <- mat2 
         #Apply the exponent to the data
        } else { tmp <- mat2^bc.exp[k] } 
      #
      tmp2 <- dagnelie.test(tmp)
      if((max(tmp2$D)-min(tmp2$D))<epsilon)
        stop("All D values are equal, Dagnelie's test cannot be computed.",
             "Check the data.", call.=FALSE)
      res[k,2] <- tmp2$Shapiro.Wilk$statistic
      res[k,3] <- tmp2$Shapiro.Wilk$p.value
      if(chord) {
        # Apply the chord transformation to matrix "tmp"
        row.norms <- apply(tmp, 1, vec.norm)
        mat3 <- sweep(tmp, 1, row.norms, "/")
        tmp2 <- dagnelie.test(mat3)
        res[k,4] <- tmp2$Shapiro.Wilk$statistic
        res[k,5] <- tmp2$Shapiro.Wilk$p.value
        }
      } 
    res
  }
    