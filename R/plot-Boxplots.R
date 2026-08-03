#' Boxplots showing the distribution of intensities over all samples
#'
#' @param D_long **data.frame** \cr
#'  The data set in long format, output object from [prepareDataSE].
#'  Must at least contain the columns ".sample" (sample IDs) and
#'  "intensity_norm" (protein intensities).
#'  May contain additional columns that can be used for assigning samples to
#'  groups for colouring.
#' @param method **character(1)** \cr
#'  Options are "boxplot" (default) and "violinplot".
#' @param groupColumn **character(1)** \cr
#'  Name of the column in D_long that contains the group information for the
#'  samples. If NULL (default), no grouping is applied.
#' @param groupColours **character** \cr
#'  Vector containing names of hex codes for the group colours (assigned in
#'  alphabetical order of the groups). If NULL (default), ggplot2 default
#'  colours are used.
#' @param baseSize **numeric(1)** \cr
#'  The base size of the font used for axis labels etc. Default is 15.
#' @param lwd **numeric(1)** \cr
#'  The line width for the boxplot, default is 0.5
#' @param outlierSize **numeric(1)** \cr
#'  The point size of the outliers. Default is 1. If 0, outliers are not shown.
#'
#' @return ggplot2 object containing the boxplot figure
#' @export
#'
#' @importFrom checkmate assertCharacter assertDataFrame assertNumeric
#' @importFrom checkmate assertSubset
#' @importFrom dplyr select
#' @importFrom ggplot2 aes element_text geom_boxplot geom_violin ggplot labs
#' @importFrom ggplot2 scale_fill_manual scale_x_discrete theme theme_bw xlab
#' @importFrom ggplot2 ylab
#' @importFrom tidyselect all_of
#'
#' @examples
#' file_proteins <- system.file("extdata", "proteins_HCC.csv",
#'   package = "ProtStatsWF")
#' file_clinical <- system.file("extdata", "clinical_data.csv",
#'   package = "ProtStatsWF")
#' D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
#'   proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
#'   sampleNameColumn = "Sample", verbose = FALSE)
#'
#' Boxplots(D_long = D$D_long, method = "boxplot", groupColumn = "Group")
Boxplots <- function(D_long,
                     method = "boxplot",
                     groupColumn = NULL,
                     groupColours = NULL,
                     baseSize = 15,
                     lwd = 0.5,
                     outlierSize = 1) {
  checkmate::assertDataFrame(D_long, min.cols = 2, min.rows = 1)
  checkmate::assertSubset(method, c("boxplot", "violinplot"))
  checkmate::assertSubset(groupColumn, choices = colnames(D_long))
  nr_groups <- length(levels(factor(D_long[, groupColumn])))
  checkmate::assertCharacter(groupColours, len = nr_groups, null.ok = TRUE)
  checkmate::assertNumeric(baseSize, lower = 0)
  checkmate::assertNumeric(lwd, lower = 0)
  checkmate::assertNumeric(outlierSize, lower = 0)

  # select only relevant columns
  D_long <- dplyr::select(D_long, c(".sample", "intensity_norm",
                                    group = tidyselect::all_of(groupColumn)))
  x_axis <- sort(unique(D_long$.sample)) # save the different states for later
  D_long <- D_long[!is.na(D_long$intensity_norm),] # remove NA values
  .sample <- intensity_norm <- group <- NULL
  if (!is.null(groupColumn)) {
   pl_boxplot <- ggplot2::ggplot(data = D_long,
      mapping = ggplot2::aes(x = .sample, y = intensity_norm, fill = group)) +
     ggplot2::labs(fill = groupColumn)
    if (!is.null(groupColours)) pl_boxplot <- pl_boxplot +
        ggplot2::scale_fill_manual(values = groupColours)
  } else {
    pl_boxplot <- ggplot2::ggplot(data = D_long,
      mapping = ggplot2::aes(x = .sample, y = intensity_norm))
  }
  pl_boxplot <- pl_boxplot +
    ggplot2::theme_bw(base_size = baseSize) +
    ggplot2::theme(axis.text.x =
                     ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggplot2::ylab("Log intensity") + ggplot2::xlab("Sample") +
    ggplot2::scale_x_discrete(limits = x_axis, drop = FALSE,
                              na.translate = TRUE)
  if (method == "violinplot") {
    pl_boxplot <- pl_boxplot + ggplot2::geom_violin(linewidth = lwd)
  }
  if (method == "boxplot") {
    pl_boxplot <- pl_boxplot +
      ggplot2::geom_boxplot(linewidth = lwd, outlier.size = outlierSize)
  }
  return(pl_boxplot)
}

















