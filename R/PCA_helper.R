#' A method for filtering the data for PCA.
#' 
#' 
#' @param D                \strong{data.frame} \cr
#'                         The data set containing intensities of the sample.
#' @param impute           \strong{logical} \cr
#'                         If \code{TRUE}, missing values will be imputed.
#' @param impute_method    \strong{character} \cr
#'                         The imputation method. Options are "mean" or "median".
#' @param propNA           \strong{numeric} \cr
#'                         The proportion of allowed missing NAs for a protein, before it is discarded.
#' 
#' 
#' @return The filtered data.
#' 

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



