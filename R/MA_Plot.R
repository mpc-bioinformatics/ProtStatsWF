#' Calculate MA plots for proteomics data
#'
#' @param sample_1 A numeric vector containing data of the first sample.
#' @param sample_2 A numeric vector containing data of the second sample.
#' @param do_log_transformation If \code{TRUE}, the data will be log-transformed.
#' @param alpha If \code{TRUE}, the data points will be transparent.
#' @param col A character containing the colors of the data points.
#' @param ... Additional arguments for affy::ma.plot.
#'
#' @return An MA plot for two samples.
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' TODO!!!
#' 
#'}
#' 
MAPlot_single <- function(sample_1, sample_2, do_log_transformation = TRUE, alpha = FALSE, col = "black", ...) {
  
  if(log) {
    sample_1 <- log2(sample_1)
    sample_2 <- log2(sample_2)
  }
  if(alpha) {
    col = alpha(col, 0.5)
  }
  
  M <- stats::na.omit(sample_1 - sample_2)
  A <- stats::na.omit((sample_1 + sample_2)/2)
  
  if (length(col) > 1) {
    na.ind <- attr(M, "na.action")
    col <- col[-na.ind]
  }
  
  
  affy::ma.plot(A = A, M = M, pch = 16, cex = 0.7, col = col, show.statistics = FALSE, ...)
}
