#' Barplots showing the percentage of valid values for each sample.
#'
#' @param D_long          \strong{data.frame} \cr
#'                        The data set given in long format.
#' @param groupColumn
#' @param group_colours   \strong{character vector} \cr
#'                        The hex codes for the group colors.
#' @param base_size       \strong{numeric} \cr
#'                        The base size of the font.
#'
#' @return A tibble and a ggplot of the valid values.
#' @export
#'
#' @examples
#'
#' path <- system.file("extdata", "exampleQC.xlsx", package = "ProtStatsWF")
#' columns <- 3:11
#'
#' prepared_data <- prepareData(data_path = path, intensity_columns = columns)
#'
#' vvplot <- ValidValuePlot(D_long = prepared_data[["D_long"]])
#'
#' plot(vvplot[["plot"]])
#'
#'
#' @importFrom magrittr %>%
#'

ValidValuePlot <- function(D_long,
                           groupColumn = NULL,
                           group_colours = NULL,
                           base_size = 15) {

  # select only relevant columns
  D_long <- dplyr::select(D_long, c(".feature", ".sample", "intensity_norm", group = groupColumn))


  #### calculate valid value table and save it ####

  .sample <- group <- intensity_norm <- nrvalid <- NULL  # initialize variables
  valid_value_table <- D_long %>%
    dplyr::group_by(.sample, group) %>%
    dplyr::summarize(nrvalid = sum(!is.na(intensity_norm)), meanvalid = mean(!is.na(intensity_norm)), .groups = 'drop')

  ### add column with sample number
  #valid_value_table$sample <- limma::strsplit2(valid_value_table$name, "_")[,2]

  #### generate basic plot skeleton ####
  valid_value_plot <- ggplot2::ggplot(valid_value_table) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ylab("Number of valid values") + ggplot2::xlab("Sample")


  #### add bars and group colours ####
  if (!is.null(groupColumn)) {
    valid_value_plot <- valid_value_plot +
      ggplot2::geom_bar(stat = "identity", ggplot2::aes(x = .sample, y = nrvalid, fill = group)) +
      ggplot2::labs(fill = groupColumn)
    if (!is.null(group_colours)) valid_value_plot <- valid_value_plot + ggplot2::scale_fill_manual(values = group_colours)
  } else {
    valid_value_plot <- valid_value_plot +
      ggplot2::geom_bar(stat = "identity",ggplot2::aes(x = .sample, y = nrvalid))
  }

  ### TODO: plot sorts samples alphabetically, not by the order it appears in the data set.

  return(list("plot" = valid_value_plot, "table" = valid_value_table))
}



