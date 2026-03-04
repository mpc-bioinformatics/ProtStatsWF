#' Barplots showing the percentage of valid values for each sample.
#'
#' @param D_long                  \strong{data.frame} \cr
#'                                The data set given in long format.
#' @param method                  \strong{character} \cr
#'                                The method used. Options are "boxplot" and "violinplot".
#' @param groupColumn              \strong{logical} \cr
#' @param groupvar_name           \strong{character} \cr
#'                                The name for the group variable.
#' @param group_colours           \strong{character vector} \cr
#'                                The hex codes for the group colors.
#' @param base_size               \strong{numeric} \cr
#'                                The base size of the font.
#' @param lwd                     \strong{numeric} \cr
#'                                The line width of the boxplot.
#' @param outlier_size            \strong{numeric} \cr
#'                                The size of the outliers.
#'
#' @return boxplots and messages
#' @export
#'
#' @examples
#' \dontrun{
#' prepared_data <- prepareData(...)
#'
#'
#' boxplot <- Boxplots(D_long = prepared_data[["D_long"]])
#' }
#'
#'

Boxplots <- function(D_long,
                     method = "boxplot",
                     groupColumn = NULL,
                     group_colours = NULL,
                     base_size = 15,
                     lwd = 0.5,
                     outlier_size = 1) {

  # select only relevant columns
  D_long <- dplyr::select(D_long, c(".feature", ".sample", "intensity_norm", group = groupColumn))

  x_axis <- sort(unique(D_long$.sample)) # save the different states for later
  D_long <- D_long[!is.na(D_long$intensity_norm),] # remove NA values


  .sample <- intensity_norm <- group <- NULL
  if (!is.null(groupColumn)) {
   pl_boxplot <- ggplot2::ggplot(data = D_long, mapping = ggplot2::aes(x = .sample, y = intensity_norm, fill = group)) +
     ggplot2::labs(fill = groupColumn)
    if (!is.null(group_colours)) pl_boxplot <- pl_boxplot + ggplot2::scale_fill_manual(values = group_colours)
  } else {
    pl_boxplot <- ggplot2::ggplot(data = D_long, mapping = ggplot2::aes(x = .sample, y = intensity_norm))
  }


  pl_boxplot <- pl_boxplot +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggplot2::ylab("Log intensity") + ggplot2::xlab("Sample") +
    ggplot2::scale_x_discrete(limits = x_axis, drop = FALSE, na.translate = TRUE)


  if (method == "violinplot") {
    pl_boxplot <- pl_boxplot + ggplot2::geom_violin()
  }
  if (method == "boxplot") {
    pl_boxplot <- pl_boxplot + ggplot2::geom_boxplot(linewidth = lwd, outlier.size = outlier_size)
  }

  return(pl_boxplot)
}

