#' Barplots showing the percentage of valid values for each sample.
#'
#' @param D_long          \strong{data.frame} \cr
#'                        The data set given in long format.
#' @param use_groups      \strong{logical} \cr
#'                        If \code{TRUE} data will be plotted in groups.
#' @param groupvar_name   \strong{character} \cr
#'                        The name for the group variable.
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
                           use_groups = NULL,
                           groupvar_name = "Group",
                           group_colours = NULL,
                           base_size = 15){


  mess <- ""

  if(is.null(use_groups)){
    if(is.na(D_long$group[1])){
      use_groups <- FALSE
    }
    else{use_groups <- TRUE}
  }
  else{
    if(use_groups && is.na(D_long$group[1])){
      use_groups <- FALSE
      }
  }


  #### calculate valid value table and save it ####

  name <- group <- value <- nrvalid <- NULL  # initialize variables
  valid_value_table <- D_long %>%
    dplyr::group_by(name, group) %>%
    dplyr::summarize(nrvalid = sum(!is.na(value)), meanvalid = mean(!is.na(value)), .groups = 'drop')

  ### add column with sample number
  #valid_value_table$sample <- limma::strsplit2(valid_value_table$name, "_")[,2]



  #### generate basic plot skeleton ####

  valid_value_plot <- ggplot2::ggplot(valid_value_table) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ylab("Number of valid values") + ggplot2::xlab("Sample")


  #### add bars and group colours ####

  if (use_groups) {
    valid_value_plot <- valid_value_plot +
      ggplot2::geom_bar(stat = "identity", ggplot2::aes(x = name,y = nrvalid, fill = group)) +
      ggplot2::labs(fill = groupvar_name)
    if (!is.null(group_colours)) valid_value_plot <- valid_value_plot + ggplot2::scale_fill_manual(values = group_colours)
    mess <- paste0(mess, "Valid Value Plot generated with groups. \n")
  } else {
    valid_value_plot <- valid_value_plot +
      ggplot2::geom_bar(stat = "identity",ggplot2::aes(x = name,y = nrvalid))
    mess <- paste0(mess, "Valid Value Plot generated without groups. \n")
  }


  ### TODO: plot sorts samples alphabetically, not by the order it appears in the data set.

  message(mess)

  return(list("plot" = valid_value_plot, "table" = valid_value_table, "message" = mess))
}



