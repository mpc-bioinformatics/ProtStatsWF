#' The main workflow for proteomics quality control
#'
#' @param data_path             A character containing the path to an .xlsx file.
#' @param intensity_columns     An integer vector containing the intensity columns of the table. 
#' @param na_strings            A character vector containing symbols to be recognized as missing values (with the exception of 0). 
#' @param zero_to_NA            If \code{TRUE}, 0 will be treated as missing value. 
#' @param do_log_transformation If \code{TRUE}, the data will be log-transformed.
#' @param log_base              A numeric containing the base used, if data is log-transformed. 
#' @param use_groups            If \code{TRUE}, the data contains groups. 
#' @param groupvar_name         A character containing the name for the group variable. 
#' @param group_colours         A character vector of hex codes for the group colors, if the data has groups. If \code{NULL}, a default color scale will be used. 
#' @param normalization_method  A character containing the method of normalization. The possible methods are no normalization "nonorm" or "median", "loess", "quantile" or "lts" normalization. 
#' @param boxplot_method        A character containing the method used for the boxplots. Possible are "boxplot" and "violinplot".
#' @param base_size             A numeric containing the base size of the font.
#'
#' @return TBD
#' @export
#'
#' @examples
#' \dontrun{
#' path <- "/Users/thisuser/Documents/dataFolder/data.xlsx" 
#' intensity_cols <- 3:17
#' 
#' result <- workflow_QC(data_path = path, intensity_columns = intensity_cols)
#'}
#' 

workflow_QC <- function(data_path, 
                        intensity_columns, 
                        na_strings = c("NA", "NaN", "Filtered","#NV"), 
                        zero_to_NA = TRUE, 
                        
                        do_log_transformation = TRUE, 
                        log_base = 2, 
                        
                        use_groups = FALSE, 
                        groupvar_name = "Group",
                        group_colours = NULL,
                        
                        normalization_method = "loess",
                        boxplot_method = "boxplot",
                        
                        base_size = 15
                        ){
  
  mess = ""
  
  prepared_data <- prepareData(data_path = data_path, intensity_columns = intensity_columns,
                               na_strings = na_strings, zero_to_NA = zero_to_NA,
                               do_log_transformation = do_log_transformation, log_base = log_base,
                               use_groups = use_groups, group_colours = group_colours,
                               normalization = normalization_method)
  mess <- paste0(mess, prepared_data[["message"]])
  
  vv_plot_data <- ValidValuePlot(D_long = prepared_data[["D_long"]],
                                          use_groups = use_groups, groupvar_name = groupvar_name, group_colours = group_colours,
                                          base_size = base_size)
  mess <- paste0(mess, vv_plot_data[["message"]])
  
  boxplot_data <- Boxplots(D_long = prepared_data[["D_long"]],
                           do_log_transformation = !do_log_transformation, log_base = log_base,
                           use_groups = use_groups, groupvar_name = groupvar_name, group_colours = group_colours,
                           base_size = base_size)
  mess <- paste0(mess, boxplot_data[["message"]])
  
  return (list("message" = mess))
}
