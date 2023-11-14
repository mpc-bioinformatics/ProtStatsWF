#' Prepare proteomics data for analysis
#' 
#' @param data_path A character containing the path to an .xlsx file.
#' @param intensity_columns An integer vector containing the intensity columns of the table.
#' @param na_strings A character vector containing symbols to be recognized as missing values (with the exception of 0).
#' @param zero_to_NA If \code{TRUE}, 0 will be treated as missing value.
#' @param log_data If \code{TRUE}, the data will be log-transformed.
#' @param log_base A numeric containing the base used, if data is log-transformed.
#' @param use_groups If \code{TRUE}, the data contains groups.
#' @param group_colours A character vector of hex codes for the group colors, if the data has groups. If \code{NULL}, a default color scale will be used.
#' @param normalization A character containing the method of normalization. The possible methods are no normalization "nonorm" or "median", "loess", "quantile" or "lts" normalization.
#' 
#' @return A list containing the prepared data and the ids of the data as data.frames as well as the groups, number of groups and group colors
#' @export
#' 
#' @examples
#'
#' 


prepareData <- function (data_path, 
                         intensity_columns, 
                         na_strings = c("NA", "NaN", "Filtered","#NV"), 
                         zero_to_NA = TRUE, 
                         log_data = TRUE, log_base = 2, 
                         use_groups = FALSE, group_colours = NULL,
                         normalization = "loess"){
  
  
  #### read and prepare data file ####
  
  D <- openxlsx::read.xlsx(data_path, na.strings = na_strings)
  
  id <- D[, -intensity_columns]
  D <- D[, intensity_columns]
  
  if(zero_to_NA) {
    D[D == 0] <- NA
  }
  
  if(log_data) {
    D <- log(D, base = log_base)
  }
  
  
  #### make data groups ####
  
  if (use_groups) {
    group <- factor(limma::strsplit2(colnames(D), "_")[,1])
  } else {
    group <- NULL
  }
  
  nr_groups <- length(levels(group))
  
  if (is.null(group_colours) & nr_groups >= 1) group_colours <- scales::hue_pal()(nr_groups)
  
  
  #### normalize the data ####
  
  D <- automatedNormalization(DATA = D, method = normalization, id_columns = id)
  
  Normalization_messages <- NULL
  Normalization_messages <- D $message
  D <- D $data
  
  
  #### calculate long form ####
  
  D_long <- tidyr::pivot_longer(data = D, cols = 1:ncol(D))
  if (use_groups) {
    D_long$group <- factor(limma::strsplit2(D_long$name, "_")[,1])
  } else {
    D_long$group <- NA
  }
  

  return (list(D, id, group, nr_groups, group_colours, D_long))
}