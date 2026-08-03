#' Barplots showing the number of valid values for each sample.
#'
#' @param D_long          **data.frame** \cr
#'  The data set in long format, output object from [prepareDataSE].
#'  Must at least contain the columns ".sample" (sample IDs) and
#'  "intensity_norm" (protein intensities).
#'  May contain additional columns that can be used for assigning samples to
#'  groups for colouring.
#' @param groupColumn **character(1)** \cr
#'  Name of the column in D_long that contains the group information for the
#'  samples. If NULL (default), no grouping is applied.
#' @param groupColours **character** \cr
#'  Vector containing names of hex codes for the group colours (assigned in
#'  alphabetical order of the groups). If NULL (default), ggplot2 default
#'  colours are used.
#' @param baseSize **numeric(1)** \cr
#'  The base size of the font used for axis labels etc. Default is 15.
#'
#' @return A list with two elements:
#' \itemize{
#'  \item plot: ggplot2 object of the valid value plot
#'  \item table: table containing the plotted numbers
#' }
#' @export
#'
#' @importFrom checkmate assertCharacter assertDataFrame assertNumeric
#' @importFrom checkmate assertSubset
#' @importFrom dplyr group_by select summarize
#' @importFrom ggplot2 aes element_text geom_bar ggplot labs theme theme_bw xlab
#' @importFrom ggplot2 ylab scale_fill_manual
#' @importFrom magrittr %>%
#'
#' @examples
#'
#' file_proteins <- system.file("extdata", "proteins_HCC.csv",
#'   package = "ProtStatsWF")
#' file_clinical <- system.file("extdata", "clinical_data.csv",
#'   package = "ProtStatsWF")
#' D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
#'   proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
#'   sampleNameColumn = "Sample", verbose = FALSE)
#'
#' ValidValuePlot(D_long = D$D_long, groupColumn = "Group")
#'
ValidValuePlot <- function(D_long,
                           groupColumn = NULL,
                           groupColours = NULL,
                           baseSize = 15) {
  checkmate::assertDataFrame(D_long, min.cols = 2, min.rows = 1)
  checkmate::assertSubset(groupColumn, choices = colnames(D_long))
  nr_groups <- length(levels(factor(D_long[, groupColumn])))
  checkmate::assertCharacter(groupColours, len = nr_groups, null.ok = TRUE)
  checkmate::assertNumeric(baseSize, lower = 0)
  # select only relevant columns
  D_long_sel <- dplyr::select(D_long, c(".sample", "intensity_norm"))
  D_long_sel <- cbind(D_long_sel, group = D_long[, groupColumn])
  sample_levels <- levels(D_long_sel$.sample)
  #### calculate valid value table ####
  .sample <- group <- intensity_norm <- nrvalid <- NULL  # initialize variables
  valid_value_table <- D_long_sel %>%
    dplyr::group_by(.sample, group) %>%
    dplyr::summarize(nrvalid = sum(!is.na(intensity_norm)),
                     meanvalid = mean(!is.na(intensity_norm)), .groups = 'drop')
  valid_value_table$.sample <- factor(valid_value_table$.sample,
                                      levels = sample_levels)
  valid_value_plot <- ggplot2::ggplot(valid_value_table) +
    ggplot2::theme_bw(base_size = baseSize) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       vjust = 1, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ylab("Number of valid values") + ggplot2::xlab("Sample")
  #### add bars and group colours ####
  if (!is.null(groupColumn)) {
    valid_value_plot <- valid_value_plot +
      ggplot2::geom_bar(stat = "identity", ggplot2::aes(x = .sample,
                        y = nrvalid, fill = group)) +
      ggplot2::labs(fill = groupColumn)
    if (!is.null(groupColours)) valid_value_plot <- valid_value_plot +
        ggplot2::scale_fill_manual(values = groupColours)
  } else {
    valid_value_plot <- valid_value_plot +
      ggplot2::geom_bar(stat = "identity",ggplot2::aes(x = .sample,
                                                       y = nrvalid))
  }
  return(list("plot" = valid_value_plot, "table" = valid_value_table))
}



