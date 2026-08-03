#' Boxplots showing the distribution of intensities over all samples
#'
#' @param D_long **data.frame** \cr
#' The data set given in long format, output object from the [prepareDataSE]
#' function. Must contain at least the columns .feature, .sample and
#' intensity_norm.
#' @param method **character(1)** \cr
#' Type of plot, options are "boxplot" (default) and "violinplot".
#' @param groupColumn **character(1)** \cr
#' Column name of D_long containing grouping information that shoul be used for
#' colouring the boxplots.
#' @param group_colours **character** \cr
#' Vector of names or hexadecimal codes for colours per group. Colours will be
#' used in alphabetical order of the groups. By default (NULL), ggplot2 standard
#'  colours are used.
#' @param base_size **integer(1)**             \strong{numeric} \cr
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
                     groupColours = NULL,
                     baseSize = 15,
                     lwd = 0.5,
                     outlierSize = 1) {

  # select only relevant columns
  D_long <- dplyr::select(D_long, c(".feature", ".sample", "intensity_norm",  ### TODO: choosable colmn name!
                                    group = tidyselect::all_of(groupColumn)))

  x_axis <- sort(unique(D_long$.sample)) # save the different states for later
  D_long <- D_long[!is.na(D_long$intensity_norm),] # remove NA values


  .sample <- intensity_norm <- group <- NULL
  if (!is.null(groupColumn)) {
   pl_boxplot <- ggplot2::ggplot(data = D_long, mapping = ggplot2::aes(x = .sample, y = intensity_norm, fill = group)) +
     ggplot2::labs(fill = groupColumn)
    if (!is.null(groupColours)) pl_boxplot <- pl_boxplot + ggplot2::scale_fill_manual(values = groupColours)
  } else {
    pl_boxplot <- ggplot2::ggplot(data = D_long, mapping = ggplot2::aes(x = .sample, y = intensity_norm))
  }


  pl_boxplot <- pl_boxplot +
    ggplot2::theme_bw(base_size = baseSize) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggplot2::ylab("Log intensity") + ggplot2::xlab("Sample") +
    ggplot2::scale_x_discrete(limits = x_axis, drop = FALSE, na.translate = TRUE)


  if (method == "violinplot") {
    pl_boxplot <- pl_boxplot + ggplot2::geom_violin()
  }
  if (method == "boxplot") {
    pl_boxplot <- pl_boxplot + ggplot2::geom_boxplot(linewidth = lwd, outlier.size = outlierSize)
  }

  return(pl_boxplot)
}

















