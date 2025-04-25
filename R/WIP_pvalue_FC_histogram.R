


#' Generates histograms for p-values, adjusted p-values and fold changes.
#'
#' @param RES               \strong{data.frame} \cr
#'                          The results from a t-test.
#' @param columnname_p      \strong{characte} \cr
#'                          The column name for p-value.
#' @param columnname_padj   \strong{characte} \cr
#'                          The columns name for adjusted p-value.
#' @param columnname_FC     \strong{characte} \cr
#'                          The column name for fold change.
#'
#' @return A list of 3 ggplot objects with a histogram for p-values, adjusted p-values and fold changes each.
#' @export
#'
#' @examples
#'

pvalue_foldchange_histogram <- function(RES,
                             columnname_p = "p",
                             columnname_padj = "padj",
                             columnname_FC = "FC",
                             base_size = 15) {

  D <- data.frame(p = RES[[columnname_p]],
                  padj = RES[[columnname_padj]],
                  FC = RES[[columnname_FC]])

<<<<<<< HEAD:R/pvalue_FC_histogram.R
  pl_hist_p <- ggplot2::ggplot(D, ggplot2::aes(x = D[["p"]])) +
=======
  p <- padj <- FC <- NULL # silence notes when checking the package

  pl_hist_p <- ggplot2::ggplot(D, ggplot2::aes(x = p)) +
>>>>>>> cdea481 (clean up to get rid of most warnings and notes when checking):R/WIP_pvalue_FC_histogram.R
    ggplot2::geom_histogram(breaks = seq(0,1, 0.05)) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::labs(title = "Histogram of p-values",
         x = "p-value",
         y = "Frequency") +
    ggplot2::xlim(0,1)

  pl_hist_padj <- ggplot2::ggplot(D, ggplot2::aes(x = D[["padj"]])) +
    ggplot2::geom_histogram(breaks = seq(0,1, 0.05)) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::labs(title = "Histogram of adjusted p-values",
         x = "Adjusted p-value",
         y = "Frequency") +
    ggplot2::xlim(0,1)

  pl_hist_FC <- ggplot2::ggplot(D, ggplot2::aes(x = log2(D[["FC"]]))) +
    ggplot2::geom_histogram() +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::labs(title = "Histogram of fold changes",
         x = "log2 Fold change",
         y = "Frequency")

  return(list("histogram_p_value" = pl_hist_p, "histogram_adjusted_p_value" = pl_hist_padj, "histogram_fold_change" = pl_hist_FC))

}


