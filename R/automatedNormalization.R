#' Automated normalization of proteomics data
#'
#' @param DATA A data.frame containing the data.
#' @param method A character containing the method of normalization. The possible methods are no normalization "nonorm" or "median", "loess", "quantile" or "lts" normalization.
#' @param id_columns An integer vector containing the columns of the data of the peptide/protein IDs.
#' @param log_transformed If \code{TRUE}, the data is log-transformed.
#' @param log_base A numeric containing the base, in case the data was log-transformed.
#' @param lts.quantile A numeric containing the quantile for the lts normalization.
#'
#' @return The normalized data as well as a message.
#' @export
#'
#' @examples
#' 

automatedNormalization <- function(DATA, 
                                   method = "loess", 
                                   id_columns = NULL,
                                   log_transformed = TRUE,
                                   log_base = 2,
                                   lts.quantile = 0.8){
  
  mess <- ""
  
  if(method == "loess" | method == "quantile" | method == "median"){
    
    #### choose normalization function
    fun <- switch(method,
                  "loess" = limma::normalizeBetweenArrays,
                  "quantile" = limma::normalizeBetweenArrays,
                  "median" = limma::normalizeBetweenArrays)
    
    ### choose arguments for normalization function
    args <- switch(method,
                   "loess" = list(object = DATA, method = "cyclicloess"),
                   "quantile" = list(object = DATA, method = "quantile"),
                   "median" = list(object = DATA, method = "scale"))
    
 
    DATA_norm <- do.call(fun, args)
    DATA_norm <- as.data.frame(DATA_norm)
    
    DATA_norm_2 <- data.frame(id_columns, DATA_norm)
    mess <- paste0(mess, "Normalized data successfully saved. \n")
    
  }
  
  if (method == "lts") {
    
    ### reverse log-transformation, if data is log-transformed
    if(log_transformed){     
      DATA <- log_base^DATA
    }
    
    DATA_norm <- vsn::vsn2(as.matrix(DATA), lts.quantile = lts.quantile)
    DATA_norm <- DATA_norm@hx
    DATA_norm <- as.data.frame(DATA_norm)
    DATA_norm_2 <- cbind(id_columns, DATA_norm)
    mess <- paste0(mess, "Normalized data successfully saved. \n")
  }
  
  if (method == "nonorm") {
    DATA_norm <- DATA
    DATA_norm_2 <- data.frame(id_columns, DATA_norm)
    mess <- paste0(mess, "Normalization method was no normalization. \n")
  }
  
  return(list("data" = DATA_norm, "message" = mess))
  
}