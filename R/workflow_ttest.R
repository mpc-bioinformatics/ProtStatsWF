


#' The workflow for t-test of quantitative proteomics data
#'
#'
#' @param data_path              A character containing the path to an .xlsx file.
#' @param output_path            A character containing the path to an output folder.
#' @param intensity_columns      An integer vector containing the intensity columns of the table.
#' @param sample                 A factor of  which column belongs to which sample.
#' @param paired                 If \code{TRUE}, a paired test will be done, otherwise an unpaired test.
#' @param var.equal              If \code{TRUE}, the variances are assumed to be equal.
#' @param log_before_test        If \code{TRUE}, the data will be log-transformed.
#' @param delog_for_FC           If \code{TRUE}, the fold change will be calculated without the log-transformation.
#' 
#' 
#' @return Message log of the workflow
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
                           intensity_columns,
                           
                           sample = NULL, 
                           paired = FALSE,
                           var.equal = FALSE,
                           log_before_test = TRUE, 
                           delog_for_FC = TRUE
                           ){
  
  mess = ""
  
  
  #### Prepare Data ####
  
  # JUST FOR TESTING PURPOSES!!! 
  output_path = "/home/kalar_ubuntu/dataresults/"
  
  data <- prepareData(data_path = "/home/kalar_ubuntu/datasets/preprocessed_peptide_data_D1_ttest.xlsx", intensity_columns = 3:8, do_log_transformation = FALSE, use_groups = TRUE)
  sample <- as.factor(c("1","2","3","1","2","3"))
  
  biomarker_data = data[["D"]][c(1,3,5,7,9,11,13),]
  biomarker_names = data[["ID"]][c(1,3,5,7,9,11,13),1]
  biomarker_group = data[["group"]][c(1,3,5,7,9,11,13)]
  
  
  #### Calculate ttest ####
  
  test_results <- ttest(D = data[["D"]], id = data[["ID"]], 
                        group = data[["group"]],  sample = sample, 
                        paired = paired, var.equal = var.equal,
                        log_before_test = log_before_test, delog_for_FC = delog_for_FC, log_base = 2,
                        min_obs_per_group = 3, min_obs_per_group_ratio = NULL,
                        filename = paste0(output_path, "results_ttest.xlsx"))
  
  
  
  
  
  
  #### Create Volcano Plot ####
  
  volcano_plot <- VolcanoPlot_ttest(RES = test_results, columnname_p = , columnname_padj = , columnname_FC = )
  
  
  
  #### Create Boxplots of Biomarker Candidates ####
  #### Create Histogram for p-values and fold changes ####
  #### Create Heatmap ####
  #### Create On-Off Heatmap ####
  
  
  return(ttest_data)
  # return(list("message" = mess))
}








#' The workflow for ANOVA of quantitative proteomics data
#'
#'
#' @param data_path              A character containing the path to an .xlsx file.
#' @param output_path            A character containing the path to an output folder.
#' @param intensity_columns      An integer vector containing the intensity columns of the table.
#' @param sample                 A factor of  which column belongs to which sample.
#' @param paired                 If \code{TRUE}, a paired test will be done, otherwise an unpaired test.
#' @param var.equal              If \code{TRUE}, the variances are assumed to be equal.
#' @param log_before_test        If \code{TRUE}, the data will be log-transformed.
#' @param delog_for_FC           If \code{TRUE}, the fold change will be calculated without the log-transformation.
#' 
#' 
#' @return Message log of the workflow
#' @export
#'
#' @examples
#' \dontrun{
#' in_path <- "/Users/thisuser/Documents/dataFolder/data.xlsx"
#' out_path <- "/Users/thisuser/Documents/resultsFolder/"
#'
#' result <- workflow_ANOVA(data_path = in_path,
#'                       output_path = out_path)
#'}

workflow_ANOVA <- function(data_path,
                           output_path,
                           intensity_columns,
                           
                           sample = NULL, 
                           paired = FALSE,
                           var.equal = FALSE,
                           log_before_test = TRUE, 
                           delog_for_FC = TRUE
){
  
  mess = ""
  
  
  #### Run ANOVA ####
  
  ANOVA_results <- ANOVA(D = data[["D"]], id = data[["ID"]], 
                         group = data[["group"]],  sample = sample, 
                         paired = paired, var.equal = var.equal,
                         log_before_test = log_before_test, delog_for_FC = delog_for_FC, log_base = 2,
                         min_obs_per_group = 3, min_perc_per_group = NULL,
                         filename = paste0(output_path, "results_ANOVA.xlsx"))
  
  
  #### Create Volcano Plot ####
  
  volcano_plot <- VolcanoPlot_ANOVA(RES = test_results, columnname_p = , columnname_padj = , columnname_FC = )
  
  
  
  #### Create Boxplots of Biomarker Candidates ####
  #### Create Histogram for p-values and fold changes ####
  #### Create Heatmap ####
  #### Create On-Off Heatmap ####
  
  
  return(ANOVA_results)
  # return(list("message" = mess))
}