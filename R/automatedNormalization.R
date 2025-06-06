#' Automated normalization of proteomics data.
#'
#' @param DATA                 \strong{data.frame} \cr
#'                             The data set containing intensities of the sample.
#' @param method               \strong{character} \cr
#'                             The method of normalization. Options are "nonorm" (no normalization), "median", "loess", "quantile" or "lts" normalization.
#' @param is_log_transformed   \strong{logical} \cr
#'                             If \code{TRUE}, the data is log-transformed.
#' @param log_base             \strong{numeric} \cr
#'                             The base, in case the data was log-transformed.
#' @param lts.quantile         \strong{numeric} \cr
#'                             The quantile for the lts normalization.
#'
#' @return The normalized data as well as a message.
#' @export
#' @examples
#'

automatedNormalization <- function(DATA,
                                   method = "median",
                                   is_log_transformed = TRUE,
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

    mess <- paste0(mess, "Data successfully ", method ," normalized. \n")

  }

  if (method == "lts") {

    ### reverse log-transformation, if data is log-transformed
    if(is_log_transformed){
      DATA <- log_base^DATA
    }

    DATA_norm <- vsn::vsn2(as.matrix(DATA), lts.quantile = lts.quantile)
    DATA_norm <- DATA_norm@hx
    DATA_norm <- as.data.frame(DATA_norm)
    mess <- paste0(mess, "Data successfully lts normalized. \n")
  }

  if (method == "nonorm") {
    DATA_norm <- DATA
    mess <- paste0(mess, "Data successfully not normalized. \n")
  }

  return(list("data" = DATA_norm, "message" = mess))

}
