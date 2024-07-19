#' The main workflow for quality control of quantitative proteomics data
#'
#'
#' @param data_path         A character containing the path to an .xlsx file.
#' @param output_path       A character containing the path to an output folder.
#' 
#' 
#' @return TODO
#' @export
#'
#' @examples
#' \dontrun{
#' in_path <- "/Users/thisuser/Documents/dataFolder/data.xlsx"
#' out_path <- "/Users/thisuser/Documents/resultsFolder/"
#'
#' result <- workflow_ttest(data_path = in_path,
#'                       output_path = out_path)
#'}

workflow_ttest <- function(data_path,
                           output_path
                           ){
  
  mess = ""
  
  
  #### Prepare Data ####
  
  # just for testing purposes
  data <- prepareData(data_path = "/home/kalar_ubuntu/datasets/preprocessed_peptide_data_D1_ttest.xlsx", intensity_columns = 3:8, do_log_transformation = FALSE, use_groups = TRUE)
  
  
  #### Calculate ttest ####
  
  ttest_data <- ttest(D = data[["D"]], id = data[["ID"]], group = data[["group"]])
  
  
  #### Calculate ANOVA ####
  
  
  #### Create Volcano Plot ####
  
  
  #### Create Histogram for p-values and fold changes ####
  
  
  return(ttest_data)
  # return(list("message" = mess))
}