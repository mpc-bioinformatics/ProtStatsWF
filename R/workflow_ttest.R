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
                           output_path,
                           
                           ttest_sample = NULL, 
                           ttest_paired = FALSE,
                           ttest_log_before_test = TRUE, 
                           ttest_delog_for_FC = TRUE
                           ){
  
  mess = ""
  
  
  #### Prepare Data ####
  
  # JUST FOR TESTING PURPOSES!!! 
  output_path = "/home/kalar_ubuntu/dataresults/"
  
  data <- prepareData(data_path = "/home/kalar_ubuntu/datasets/preprocessed_peptide_data_D1_ttest.xlsx", intensity_columns = 3:8, do_log_transformation = FALSE, use_groups = TRUE)
  ttest_sample <- as.factor(c("1","2","3","1","2","3"))
  
  biomarker_data = data[["D"]][c(1,3,5,7,9,11,13),]
  biomarker_names = data[["ID"]][c(1,3,5,7,9,11,13),1]
  biomarker_group = data[["group"]][c(1,3,5,7,9,11,13)]
  
  
  
  
  
  
  #### Calculate ttest ####
  
  ttest_data <- ttest(D = data[["D"]], id = data[["ID"]], 
                      group = data[["group"]],  sample = ttest_sample, 
                      paired = ttest_paired, var.equal = FALSE,
                      log_before_test = ttest_log_before_test, delog_for_FC = ttest_delog_for_FC, log_base = 2,
                      min_obs_per_group = 3,
                      filename = paste0(output_path, "results_ttest.xlsx"),
                      min_obs_per_group_ratio = NULL)
  
  
  #### Create Boxplots of Biomarker Candidates ####
  
  Boxplots_candidates(D = biomarker_data, protein.names = biomarker_names, output_path = output_path, group = biomarker_group)
  
  
  #### Calculate ANOVA ####
  #### Create Volcano Plot ####
  #### Create Histogram for p-values and fold changes ####
  #### Create Heatmap ####
  #### Create On-Off Heatmap ####
  
  
  return(ttest_data)
  # return(list("message" = mess))
}