filter_PCA_data <- function(D, impute = FALSE, impute_method = "mean", propNA = 0){
  
  # proportion if missing values per protein
  mean_NA <- apply(D, 1, function(x) mean(is.na(x)))

  ### remove rows with too many missing values
  D <- D[mean_NA <= propNA, ]
  
  if (nrow(D) == 0){
    return(NULL)
  }
  
  ## perform imputation
  if (impute) {
    D <- as.data.frame(t(apply(D, 1, function(x) {
      if(anyNA(x)) {
        x[is.na(x)] <- switch(impute_method, mean = mean(x, na.rm = TRUE),
                              median = stats::median(x, na.rm = TRUE))
      }
      return(x)
    })))
    
    
  } else {
    D <- stats::na.omit(D)
    if (nrow(D) == 0){
      return(NULL)
    }
  }
  
  ### remove proteins/peptides with an (almost) constant value (variance near zero)
  v <- matrixStats::rowVars(as.matrix(D))
  ind_zeroVar <- which(v < 1e-25)
  if (length(ind_zeroVar) > 0) D <- D[-ind_zeroVar,]
  
  return(D)
}



