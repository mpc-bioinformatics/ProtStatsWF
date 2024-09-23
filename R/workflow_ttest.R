#' The workflow for t-test of quantitative proteomics data
#'
#'
#' @param data_path              A character containing the path to an .xlsx file.
#' @param output_path            A character containing the path to an output folder.
#' @param intensity_columns      An integer vector containing the intensity columns of the table.
#' @param paired                 If \code{TRUE}, a paired test will be done, otherwise an unpaired test.
#' @param var.equal              If \code{TRUE}, the variances are assumed to be equal.
#' @param log_before_test        If \code{TRUE}, the data will be log-transformed.
#' @param delog_for_FC           If \code{TRUE}, the fold change will be calculated without the log-transformation.
#' 
#' @param suffix                 A character if the filenames should contain a suffix.
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
                           
                           paired = FALSE,
                           var.equal = FALSE,
                           log_before_test = TRUE, 
                           delog_for_FC = TRUE,
                           
                           max_valid_values_off = NULL,
                           min_valid_values_on = NULL,
                           
                           suffix = "",
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
                        group = data[["group"]],  sample = data[["sample"]], 
                        paired = paired, var.equal = var.equal,
                        log_before_test = log_before_test, delog_for_FC = delog_for_FC, log_base = 2,
                        min_obs_per_group = 3, min_obs_per_group_ratio = NULL,
                        filename = paste0(output_path, "results_ttest", suffix, ".xlsx"))
  
  mess <- paste0(mess, 
                 ifelse(paired, "Paired", "Unaired"), 
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
  
  ggplot2::ggsave(paste0(output_path, "volcano_plot", suffix, ".", plot_device), plot = volcano_plot,
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  
  mess <- paste0(mess, "Volcano plot calculated. \n")
  
  
  
  #### Create Histogram for p-values and fold changes ####
  
  histograms <- pvalue_foldchange_histogram(RES = test_results, 
                                            columnname_p = "p", columnname_padj = "p.fdr", 
                                            columnname_FC = "FC_state1_divided_by_state2")
  
  ggplot2::ggsave(paste0(output_path, "histogram_p_value", suffix, ".", plot_device), plot = histograms[["histogram_p_value"]],
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  ggplot2::ggsave(paste0(output_path, "histogram_adjusted_p_value", suffix, ".", plot_device), plot = histograms[["histogram_adjusted_p_value"]],
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  ggplot2::ggsave(paste0(output_path, "histogram_fold_change", suffix, ".", plot_device), plot = histograms[["histogram_fold_change"]],
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  
  mess <- paste0(mess, "p-value, adjusted p-value and fold change histograms calculated. \n")
  
  
  #### Get significant candidates ####
  
  candidates <- as.character(calculate_significance_categories_ttest(p = test_results[["p"]], 
                                                        p_adj = test_results[["p.fdr"]],
                                                        fc = test_results[["FC_state1_divided_by_state2"]]))

  candidates <- which(candidates == "significant after FDR correction")
  
  mess <- paste0(mess, "There are ", length(candidates), " candidates, which were significant after FDR correction. \n")
  
  
  #### Create Boxplots of Biomarker Candidates ####

  Boxplots_candidates(D = data[["D"]][candidates, ], 
                      protein.names = data[["ID"]][candidates, "protein"],
                      group = data[["group"]],
                      suffix = suffix,
                      output_path = paste0(output_path))
  
  mess <- paste0(mess, "Boxplots made for the candidates. \n")
  
  
  
  #### Create Heatmap ####
  
  t_heatmap <- Heatmap_with_groups(D = data[["D"]][candidates, ], 
                                   id = data[["ID"]][candidates, ],
                                   groups = data[["group"]])
  
  grDevices::pdf(paste0(output_path, "/heatmap", suffix, ".pdf"), height = plot_height, width = plot_width)
  graphics::plot(t_heatmap)
  grDevices::dev.off()
  
  mess <- paste0(mess, "Heatmap made for the candidates. \n")
  
  
  
  #### Create On-Off Heatmap ####
  
  
  return(list("message" = mess))
}








#' The workflow for ANOVA of quantitative proteomics data
#'
#'
#' @param data_path              A character containing the path to an .xlsx file.
#' @param output_path            A character containing the path to an output folder.
#' @param intensity_columns      An integer vector containing the intensity columns of the table.
#' @param paired                 If \code{TRUE}, a paired test will be done, otherwise an unpaired test.
#' @param var.equal              If \code{TRUE}, the variances are assumed to be equal.
#' @param log_before_test        If \code{TRUE}, the data will be log-transformed.
#' @param delog_for_FC           If \code{TRUE}, the fold change will be calculated without the log-transformation.
#' 
#' @param suffix                 A character if the filenames should contain a suffix.
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
                           
                           paired = FALSE,
                           var.equal = FALSE,
                           log_before_test = TRUE, 
                           delog_for_FC = TRUE,
                           
                           suffix = ""
){
  
  mess = ""
  
  
  #### Prepare Data ####
  
  # data <- prepareTtestData(data_path = "C:/Users/kalar/Documents/0_Studium/WHK/Testdata/ANOVA/preprocessed_peptide_data_D1.xlsx" , intensity_columns = 3:17)
  # output_path = "C:/Users/kalar/Documents/0_Studium/WHK/Testdata/ANOVA/results/"
  data <- prepareTtestData(data_path = data_path , intensity_columns = intensity_columns)
  
  
  
  #### Calculate ANOVA ####
  
  ANOVA_results <- ANOVA(D = data[["D"]], id = data[["ID"]], 
                         group = data[["group"]],  sample = data[["sample"]], 
                         paired = paired, var.equal = var.equal,
                         log_before_test = log_before_test, delog_for_FC = delog_for_FC, log_base = 2,
                         min_obs_per_group = 3, min_perc_per_group = NULL,
                         filename = paste0("results_ANOVA.xlsx", suffix))
  
  
  #### Create Volcano Plot ####
  
  # volcano_plot <- VolcanoPlot_ANOVA(RES = test_results, columnname_p = , columnname_padj = , columnname_FC = )
  
  # mess <- paste0(mess, "Volcano plot calculated. \n")
  
  # ggplot2::ggsave(paste0(output_path, "volcano_plot", ".", plot_device), plot = volcano_plot,
  #                        device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  
  
  
  
  
  
  #### Create Boxplots of Biomarker Candidates ####
  #### Create Heatmap ####
  #### Create On-Off Heatmap ####
  
  
  return(ANOVA_results)
  # return(list("message" = mess))
}