
#' Barplots showing the percentage of valid values for each sample
#'
#' @param D_long the data set given in long format
#' @param output_path output path
#' @param suffix suffix for file name
#' @param groupvar_name name for the group variable
#' @param group_colours colours for the groups
#' @param plot_device device
#' @param plot_height height of the plot
#' @param plot_width width of the plot
#' @param plot_dpi resolution of the plot
#' @param base_size base size of the font size
#'
#' @return ggplot, saves a table
#' @export
#'
#' @examples
#' 
#' @importFrom magrittr %>%
#' 
ValidValuePlot <- function(D_long,
                           output_path = "", suffix = "",
                           groupvar_name = "Group", group_colours = NULL,
                           plot_device = "png",
                           plot_height = 10, plot_width = 15, plot_dpi = 300,
                           base_size = 15){


  use_groups <- !all(is.na(D_long$group))

  ## calculate valid value table and save it
  name <- group <- value <- nrvalid <- NULL
  X <- D_long %>% dplyr::group_by(name, group) %>% dplyr::summarize(nrvalid = sum(!is.na(value)), meanvalid = mean(!is.na(value)), .groups = 'drop')
  openxlsx::write.xlsx(x = X, file = paste0(output_path, "validvalues", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)


  ## generate basic plot skeleton
  pl_valid_values <- ggplot2::ggplot(X) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
          plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ylab("Number of valid values") + ggplot2::xlab("Sample")

  # add bars and group colours
  if (use_groups) {
    pl_valid_values <- pl_valid_values +
      ggplot2::geom_bar(stat = "identity", ggplot2::aes(x = name,y = nrvalid, fill = group)) +
      ggplot2::labs(fill = groupvar_name)
    if (!is.null(group_colours)) pl_valid_values <- pl_valid_values + ggplot2::scale_fill_manual(values = group_colours)
  } else {
    pl_valid_values <- pl_valid_values +
      ggplot2::geom_bar(stat = "identity",ggplot2::aes(x = name,y = nrvalid))
  }

  # save plot
  ggplot2::ggsave(paste0(output_path,"valid_value_plot", suffix,".",plot_device),
         plot = pl_valid_values, device = plot_device,
         height = plot_height, width = plot_width, dpi = plot_dpi)

  return(pl_valid_values)
}



