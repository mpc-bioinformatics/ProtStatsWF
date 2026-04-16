#' A method for filtering the data for PCA.
#'
#'
#' @param SE                \strong{SummarizedExperiment object} \cr
#'                         The data set containing intensities of the sample.
#' @param assay            \strong{character(1)} \cr
#'                         The name of the assay in SE containing the protein intensities.
#'                         Default is "intensity_norm", which uses the normalized
#'                         protein intensities if the SE object was generated
#'                         by the ProtStatsWF functionality.
#' @param imputeMethod    \strong{character} \cr
#'                         The imputation method. Options are "mean" or "median" or "none
#' @param propNA           \strong{numeric} \cr
#'                         The proportion of allowed missing NAs for a protein, before it is discarded. It is set to 0 if imputeMethod is "none".
#'
#'
#' @return The filtered data.
#' @export
#' @examples
#'
#'
#'

filter_PCA_data <- function(SE, assay = "intensity_norm",
                            imputeMethod = "none", propNA = 0) {

  # proportion of missing values per protein
  D_all <- SummarizedExperiment::assays(SE)[[assay]]
  mean_NA <- apply(D_all, 1, function(x) mean(is.na(x)))

  ### remove rows with too many missing values
  if (imputeMethod == "none") propNA <- 0
  index_to_keep <- mean_NA <= propNA
  SE <- SE[index_to_keep, ]

  if (nrow(SE) == 0){
    return(NULL)
  }

  ## remove or impute missing values
  #D <- SummarizedExperiment::assays(SE)$intensity_norm
  D <- SummarizedExperiment::assays(SE)[[assay]]
  if (imputeMethod != "none") {
    D <- as.data.frame(t(apply(D, 1, function(x) {
      if(anyNA(x)) {
        x[is.na(x)] <- switch(imputeMethod, mean = mean(x, na.rm = TRUE),
                              median = stats::median(x, na.rm = TRUE))
      }
      return(x)
    })))


    SummarizedExperiment::assays(SE)[[assay]] <- D
  }

  ### remove proteins/peptides with an (almost) constant value (variance near zero)
  v <- matrixStats::rowVars(as.matrix(D))

  ind_zeroVar <- (v < 1e-25)
  SE <- SE[!ind_zeroVar, ]
  return(SE)
}



