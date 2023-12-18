#' Barplots showing the percentage of valid values for each sample
#'
#' @param D_long A data.frame of the data set given in long format.
#' @param do_log_transformation If \code{TRUE}, the data is log-transformed
#' @param log_base A numeric containing the base used, if data is log-transformed.
#' @param method A character containing the method used. Possible are "boxplot" and "violinplot".
#' @param use_groups If \code{TRUE} data will be plotted in groups.
#' @param groupvar_name A character containing the name for the group variable.
#' @param group_colours A character vector of hex codes for the group colors.
#' @param base_size A numeric containing the base size of the font.
#'
#' @return a tibble and a ggplot of the valid values
#' @export
#'
#' @examples
#' \dontrun{
#' prepared_data <- prepareData(data_path = path, intensity_columns = intensity_cols)
#' 
#' boxplot <- Boxplots(D_long = prepared_data[["D_long"]])
#' }
#' 

Boxplots <- function(D_long,
                     do_log_transformation = FALSE, 
                     log_base = 2,
                     method = "boxplot",
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
  
  # log-transform data if necessary
  if(do_log_transformation) {
    D_long$value <- log(D_long$value, base = log_base)
  }
  
  x_axis <- sort(unique(D_long$name)) # save the different states for later
  D_long <- D_long[!is.na(D_long$value),] # remove NA values
  
  name <- value <- group <- NULL
  if (use_groups) {
   pl_boxplot <- ggplot2::ggplot(D_long, ggplot2::aes(x = name, y = value, fill = group)) +
     ggplot2::labs(fill = groupvar_name) 
    if (!is.null(group_colours)) pl_boxplot <- pl_boxplot + ggplot2::scale_fill_manual(values = group_colours)
    mess <- paste0(mess, "with groups. \n")
  } else {
    pl_boxplot <- ggplot2::ggplot(D_long, ggplot2::aes(x = name, y = value))
    mess <- paste0(mess, "without groups. \n")
  }
  
  
  pl_boxplot <- pl_boxplot +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
    ggplot2::ylab("Log intensity") + ggplot2::xlab("Sample") + 
    ggplot2::scale_x_discrete(limits = x_axis, drop = FALSE, na.translate = TRUE)
  
  
  if (method == "violinplot"){
    pl_boxplot <- pl_boxplot + ggplot2::geom_violin()
    mess <- paste0("Violin Plot generated ", mess)
  }
  
  if (method == "boxplot") {
    pl_boxplot <- pl_boxplot + ggplot2::geom_boxplot()
    mess <- paste0("Boxplot generated ", mess)
  }
  
  message(mess)
  
  return(list("plot" = pl_boxplot, "message" = mess))
}