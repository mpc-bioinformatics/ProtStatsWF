#' Automated normalization of proteomics data
#'
#' @param DATA A data.frame containing the data.
#' @param method A character containing the method of normalization. The possible methods are no normalization "nonorm" or "median", "loess", "quantile" or "lts" normalization.
#' @param id_columns An integer vector containing the columns of the data of the peptide/protein IDs.
#' @param groupwise If \code{TRUE}, the data contains groups.
#' @param group The groups of the data, if it has any.
#' @param lts.quantile A numeric containing the quantile for the lts normalization.
#'
#' @return The normalized data as well as a message.
#' @export
#'
#' @examples
#' 

automatedNormalization <- function(DATA, method = "loess", id_columns = NULL,
                                   groupwise = FALSE, group = NULL, lts.quantile = 0.8){
  
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
    
    if (!groupwise) {
      DATA_norm <- do.call(fun, args)
      DATA_norm <- as.data.frame(DATA_norm)
    } else {
      DATA_split <- split.default(DATA, group)
      DATA_split_norm <- lapply(DATA_split, limma::normalizeBetweenArrays, method = args$method)
      DATA_norm <- do.call(cbind, DATA_split_norm)
    }
    
    DATA_norm_2 <- data.frame(id_columns, DATA_norm)
    mess <- paste0(mess, "Normalized data successfully saved! \n")
    
  }
  
  ### TODO: Daten müssen vor dieser Normalisierung nicht log-transformiert werden!
  ### TODO: groupwise normalisation?
  if (method == "lts") {
    
    DATA_norm <- vsn::vsn2(as.matrix(DATA), lts.quantile = lts.quantile)
    DATA_norm <- DATA_norm@hx
    DATA_norm <- as.data.frame(DATA_norm)
    
    DATA_norm_2 <- cbind(id_columns, DATA_norm)
    message("Normalized data successfully saved!")
    mess <- paste0(mess, "Normalized data successfully saved! \n")
  }
  
  if (method == "nonorm") { # if method == "nonorm"
    DATA_norm <- DATA
    DATA_norm_2 <- data.frame(id_columns, DATA_norm)
    
    message("Normalized data successfully saved!")
    mess <- paste0(mess, "Normalized data successfully saved! \n")
  }
  
  return(list("data" = DATA_norm, "message" = mess))
  
}