#' Prepare data from an xlsx sheet
#'
#'
#' @param data_path              A character containing the path to an .xlsx file.
#' @param intensity_columns      An integer vector containing the intensity columns of the table.
#' 
#' @return A list containing an intensity data.frame, an IDs data frame and a factor of the sample groups
#'
#' @examples
#' \dontrun{
#' in_path <- "/Users/thisuser/Documents/dataFolder/data.xlsx"
#' int_cols <- 3:8
#'
#' result <- workflow_ttest(data_path = in_path, intensity_columns = int_cols)
#'}

prepareTtestData <- function(data_path,
                             intensity_columns
){
  
  D <- openxlsx::read.xlsx(data_path)
  
  id <- D[, -intensity_columns]
  D <- D[, intensity_columns]
  
  group <- factor(limma::strsplit2(colnames(D), "_")[,1])
  number_of_groups <- length(levels(group))
  
  
  return(list("D" = D, "ID" = id, "group" = group, "number_of_groups" = number_of_groups))
}




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
#' @param plot_device            A character containing the type of the output file, e.g. "pdf" or "png".
#' @param plot_height            A numeric of the plot height in cm.
#' @param plot_width             A numeric of the plot width in cm.
#' @param plot_dpi               A numeric of the "dots per inch" of the plot aka. the plot resolution.
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
                           delog_for_FC = TRUE,
                           
                           plot_device = "pdf",
                           plot_height = 15,
                           plot_width = 15,
                           plot_dpi = 300
                           ){
  
  mess = ""
  
  
  #### Prepare Data ####
  
  # data <- prepareTtestData(data_path = "C:/Users/kalar/Documents/0_Studium/WHK/Testdata/ttest/preprocessed_peptide_data_D3.xlsx" , intensity_columns = 3:8)
  # output_path = "C:/Users/kalar/Documents/0_Studium/WHK/Testdata/ttest/results/"
  data <- prepareTtestData(data_path = data_path , intensity_columns = intensity_columns)
  
  
  
  #### Calculate ttest ####
  
  test_results <- ttest(D = data[["D"]], id = data[["ID"]], 
                        group = data[["group"]],  sample = sample, 
                        paired = paired, var.equal = var.equal,
                        log_before_test = log_before_test, delog_for_FC = delog_for_FC, log_base = 2,
                        min_obs_per_group = 3, min_obs_per_group_ratio = NULL,
                        filename = paste0(output_path, "results_ttest.xlsx"))
  
  mess <- paste0(mess, 
                 ifelse(paired, "Unpaired", "Paired"), 
                 " t-test calculated with the variance assumed to be ", 
                 ifelse(var.equal, "equal", "unequal"), ". \n",
                 "Data was ", 
                 ifelse(log_before_test, "", "not"), 
                 "log-transformed before the t-test and ", 
                 ifelse(delog_for_FC, "", "not"),
                 "de-log-transformed for the fold change. \n")
  
  
  
  
  #### Create Volcano Plot ####
  
  volcano_plot <- VolcanoPlot_ttest(RES = test_results, 
                                    columnname_p = "p", columnname_padj = "p.fdr", 
                                    columnname_FC = "FC_state1_divided_by_state2")
  
  mess <- paste0(mess, "Volcano plot calculated. \n")
  
  ggplot2::ggsave(paste0(output_path, "volcano_plot", ".", plot_device), plot = volcano_plot,
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  
  
  
  #### Create Histogram for p-values and fold changes ####
  
  histograms <- pvalue_foldchange_histogram(RES = test_results, 
                                            columnname_p = "p", columnname_padj = "p.fdr", 
                                            columnname_FC = "FC_state1_divided_by_state2")
  
  mess <- paste0(mess, "p-value, adjusted p-value and fold change histograms calculated. \n")
  
  ggplot2::ggsave(paste0(output_path, "histogram_p_value", ".", plot_device), plot = histograms[["histogram_p_value"]],
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  ggplot2::ggsave(paste0(output_path, "histogram_adjusted_p_value", ".", plot_device), plot = histograms[["histogram_adjusted_p_value"]],
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  ggplot2::ggsave(paste0(output_path, "histogram_fold_change", ".", plot_device), plot = histograms[["histogram_fold_change"]],
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  
  
  
  
  
  #### Create Boxplots of Biomarker Candidates ####
  #### Create Heatmap ####
  #### Create On-Off Heatmap ####
  
  
  return(list("message" = mess))
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
  
  
  #### Prepare Data ####
  
  # data <- prepareTtestData(data_path = "C:/Users/kalar/Documents/0_Studium/WHK/Testdata/ANOVA/preprocessed_peptide_data_D1.xlsx" , intensity_columns = 3:17)
  # output_path = "C:/Users/kalar/Documents/0_Studium/WHK/Testdata/ANOVA/results/"
  data <- prepareTtestData(data_path = data_path , intensity_columns = intensity_columns)
  
  
  
  #### Calculate ANOVA ####
  
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