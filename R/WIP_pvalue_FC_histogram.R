


#' Generates histograms for p-values, adjusted p-values and fold changes
#'
#' @param RES result table from t-test
#' @param columnname_p column name for p-value
#' @param columnname_padj columns name for adjusted p-value
#' @param columnname_FC column name for fold change
#'
#' @return list of 3 ggplot objects
#' @export
#'
#' @examples # TODO
pvalue_foldchange_histogram <- function(RES,
                             columnname_p = "p",
                             columnname_padj = "padj",
                             columnname_FC = "FC") {

  D <- data.frame(p = RES[[columnname_p]],
                  padj = RES[[columnname_padj]],
                  FC = RES[[columnname_FC]])

  p <- padj <- FC <- NULL # silence notes when checking the package

  pl_hist_p <- ggplot2::ggplot(D, ggplot2::aes(x = p)) +
    ggplot2::geom_histogram(breaks = seq(0,1, 0.05)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Histogram of p-values",
         x = "p-value",
         y = "Frequency") +
    ggplot2::xlim(0,1)

  pl_hist_padj <- ggplot2::ggplot(D, ggplot2::aes(x = padj)) +
    ggplot2::geom_histogram(breaks = seq(0,1, 0.05)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Histogram of p-values",
         x = "Adjusted p-value",
         y = "Frequency") +
    ggplot2::xlim(0,1)

  pl_hist_FC <- ggplot2::ggplot(D, ggplot2::aes(x = log2(FC))) +
    ggplot2::geom_histogram() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Histogram of p-values",
         x = "log2 Fold change",
         y = "Frequency")

  return(list(pl_hist_p, pl_hist_padj, pl_hist_FC))

}


