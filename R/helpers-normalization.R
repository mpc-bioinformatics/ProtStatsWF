#' Automated normalization of proteomics data.
#'
#' @param DATA **data.frame** \cr
#' Data set containing protein or peptide intensities of the samples and no
#' additional columns.
#' @param method **character(1)** \cr
#' The method of normalization. Options are "nonorm" (no normalization),
#' "median", "loess" (default), "quantile" or "lts" normalization.
#' @param is_log_transformed **logical(1)** \cr
#' Set to `TRUE`, if the data is already log-transformed. Only relevant for lts
#' normalization.
#' @param log_base **numeric(1)** \cr
#' The base, in case the data was log-transformed before normalization. Default
#' is 2.
#' @param lts.quantile **numeric(1)** \cr
#' The quantile for the lts normalization. For details check [vsn::vsn2].
#' @param verbose **logical(1)** \cr
#' If `TRUE` (default), messages are printed out.
#'
#' @return The normalized data as a data.frame.
#'
#' @importFrom checkmate assertDataFrame assertFlag assertNumber assertSubset
#' @importFrom limma normalizeBetweenArrays
#' @importFrom vsn vsn2
.normalization <- function(DATA,
                           method = "loess",
                           is_log_transformed = TRUE,  ##TODO: only used for lts, unclear if I can always get that info right.
                           log_base = 2,
                           lts.quantile = 0.8,
                           verbose = TRUE) {
  checkmate::assertDataFrame(DATA)
  checkmate::assertSubset(method, choices = c("nonorm", "median", "loess", "quantile", "lts"))
  checkmate::assertFlag(is_log_transformed)
  checkmate::assertNumber(log_base, lower = 1)
  checkmate::assertNumber(lts.quantile, lower = 0, upper = 1)
  checkmate::assertFlag(verbose)

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

    if (verbose) message("Data successfully ", method ," normalized.")
  }

  if (method == "lts") {

    ### reverse log-transformation, if data is log-transformed
    if(is_log_transformed){
      DATA <- log_base^DATA
    }

    DATA_norm <- vsn::vsn2(as.matrix(DATA), lts.quantile = lts.quantile)
    DATA_norm <- DATA_norm@hx
    DATA_norm <- as.data.frame(DATA_norm)
    if (verbose) message("Data successfully lts normalized with lts.quantile = ",
                         lts.quantile, ".")
  }

  if (method == "nonorm") {
    DATA_norm <- DATA
    if (verbose) message("No normalization.")
  }

  return(DATA_norm)
}
