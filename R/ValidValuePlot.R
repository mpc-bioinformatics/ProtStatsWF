#' Barplots showing the percentage of valid values for each sample
#'
#' @param D_long A data.frame of the data set given in long format.
#' @param groupvar_name A character containing the name for the group variable.
#' @param group_colours A character vector of hex codes for the group colors.
#' @param base_size A numeric containing the base size of the font.
#'
#' @return a tibble and a ggplot of the valid values
#' @export
#'
#' @examples
#' 
#' @importFrom magrittr %>%
#' 
ValidValuePlot <- function(D_long,
                           groupvar_name = "Group", group_colours = NULL,
                           base_size = 15){


  mess <- ""
  use_groups <- !all(is.na(D_long$group))

  
  #### calculate valid value table and save it ####
  
  name <- group <- value <- nrvalid <- NULL  # initialize variables 
  valid_value_table <- D_long %>% dplyr::group_by(name, group) %>% dplyr::summarize(nrvalid = sum(!is.na(value)), meanvalid = mean(!is.na(value)), .groups = 'drop')


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
    mess <- paste0(mess, "Valid Value Plot with groups generated. \n")
  } else {
    valid_value_plot <- valid_value_plot +
      ggplot2::geom_bar(stat = "identity",ggplot2::aes(x = name,y = nrvalid))
    mess <- paste0(mess, "Valid Value Plot generated without groups. \n")
  }

  return(list("table" = valid_value_table, "plot" = valid_value_plot, "message" = mess))
}



