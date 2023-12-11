#' Title
#'
#' @param D              A data.frame of the data set.
#' @param id             A data.frame of the according IDs of the data set.
#' @param groupvar1      A variable used for colors.
#' @param groupvar2      A variable used for shapes.
#' @param impute         If \code{TRUE}, missing values will be imputed.
#' @param impute_method  A character containing the imputation method ("mean" or "median")
#' @param propNA         A numeric of the proportion of allowed missing NAs for a protein, before it is discarded. 
#' @param scale.         If \code{TRUE}, the data will be scaled before computing the PCA.
#' @param PCx            The principle component for the x-axis (default: 1).
#' @param PCy            The principle component for the y-axis (default: 2).
#' @param groupvar1_name Titles of legends for colour and shape.
#' @param groupvar2_name Titles of legends for colour and shape.
#' @param group_colours  The colors for the groups.
#' @param alpha          If \code{TRUE}, the data points will be transparent.
#' @param label          If \code{TRUE}, the samples will be labeled.
#' @param label_size     A numeric containing the size of the sample labels.
#' @param xlim           Limit of the x-axis.
#' @param ylim           Limit of the y-axis.
#' @param point.size     The size of the data points.
#' @param base_size      The base size for the plot
#' @param plot_height    A numeric containing the plot height in cm.
#' @param plot_width     A numeric containing the plot width in cm.
#' @param plot_dpi       A numeric containing the "dots-per-inch" for the plot.
#' @param ...            Additional arguments for ggplot2::plot.
#'
#' @return The PCA plot for a proteomics sample 
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' TODO!!!
#' 
#'}
#' 

PCA_Plot <- function(D, id,
                     groupvar1, groupvar2 = NULL,
                     
                     impute = FALSE, impute_method = "mean", propNA = 0,
                     scale. = TRUE,
                     PCx = 1, PCy = 2,
                     
                     groupvar1_name = "group", groupvar2_name = NULL,
                     group_colours = NULL, alpha = 1,
                     label = FALSE, label_size = 4,
                     xlim = NULL, ylim = NULL,
                     
                     point.size = 4, base_size = 11,
                     plot_height = 10, plot_width = 10, plot_dpi = 300,
                     ...
                     
) {
  
  mess = ""
  
  # TODO: maybe give as an argument?
  use_groups <- !is.null(groupvar1)
  
  filtered_D <- filter_PCA_data(D)
  
  
  
  
  
  
  
  return(list("message" = mess))
}