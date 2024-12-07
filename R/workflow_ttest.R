#' The workflow for a t-test of quantitative proteomics data.
#'
#' @param data_path              \strong{character} \cr
#'                               The path to an .xlsx file containing the input data.
#' @param output_path            \strong{character} \cr
#'                               The path to the output folder.
#' @param intensity_columns      \strong{integer} \cr
#'                               A vector containing the numbers of the intensity columns in the table.
#' 
#' @param paired                 \strong{logical} \cr
#'                               If \code{TRUE}, a paired test will be done, otherwise an unpaired test.
#' @param var.equal              \strong{logical} \cr
#'                               If \code{TRUE}, the variances are assumed to be equal.
#' @param log_before_test        \strong{logical} \cr
#'                               If \code{TRUE}, the data will be log-transformed.
#' @param delog_for_FC           \strong{logical} \cr
#'                               If \code{TRUE}, the fold change will be calculated without the log-transformation.
#' 
#' @param significant_after_FDR  \strong{logical} \cr
#'                               If \code{TRUE}, candidates for the boxplots and heatmap need to be significant after FDR correction, otherwise all significant candidates will be used.
#' @param max_valid_values_off   \strong{integer} \cr
#'                               The maximum number of valid values to be an off protein.
#' @param min_valid_values_on    \strong{integer} \cr
#'                               The minimum number of valid values to be an on protein.
#' 
#' @param suffix                 \strong{character} \cr
#'                               The suffix of the file names should have one.
#' @param plot_device            \strong{character} \cr
#'                               The type of the output file, e.g. "pdf" or "png".
#' @param plot_height            \strong{numeric} \cr
#'                               The plot height in cm.
#' @param plot_width             \strong{numeric} \cr
#'                               The plot width in cm.
#' @param plot_dpi               \strong{integer} \cr
#'                               The "dots per inch" of the plot aka. the plot resolution.
#' 
#' 
#' @return Returns a message log of the workflow. The log contains an overview of the settings and gives some information e.g number of significant candidates or on-off-proteins.
#' @export
#' 
#' @seealso [workflow_ANOVA()] in case of more than two groups in the sample.\cr
#'          Functions used in this workflow: 
#'          [prepareTtestData()], [ttest()], [VolcanoPlot_ttest()], [pvalue_foldchange_histogram()], 
#'          [calculate_significance_categories_ttest()], [Boxplots_candidates()], 
#'          [Heatmap_with_groups()], [calculate_onoff()]
#'
#' @examples
#' 
#' # 1. Set the character of your data path, leading to an .xlsx file.
#' in_path <- "C:/Users/thisuser/Documents/dataFolder/data.xlsx"
#' 
#' # 2. Set the integer vector of the columns, which contain the intensities.
#' int_col <- 3:8
#' 
#' # 3. Set the character of the output path, leading to a folder for the results.
#' out_path <- "C:/Users/thisuser/Documents/resultsFolder/"
#' 
#' # 4. Run the ttest with the parameters you set.
#' \dontrun{
#' result <- workflow_ttest(data_path = in_path,
#'                          output_path = out_path,
#'                          intensity_columns = int_col) }
#'                          

workflow_ttest <- function(data_path,
                           output_path,
                           intensity_columns,
                           
                           paired = FALSE,
                           var.equal = FALSE,
                           log_before_test = TRUE,
                           delog_for_FC = TRUE,
                           
                           significant_after_FDR = TRUE,
                           max_valid_values_off = 0,
                           min_valid_values_on = NULL,
                           
                           suffix = "",
                           plot_device = "pdf",
                           plot_height = 15,
                           plot_width = 15,
                           plot_dpi = 300
                           ){
  
  mess = ""
  
  
  #### Prepare Data ####
  
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
  
  fc_col_name <- paste0("FC_", levels(data[["group"]])[[1]], "_divided_by_", levels(data[["group"]])[[2]])
  
  
  
  #### Create Volcano Plot ####
  
  volcano_plot <- VolcanoPlot_ttest(RES = test_results, 
                                    columnname_p = "p", columnname_padj = "p.fdr", 
                                    columnname_FC = fc_col_name)
  
  ggplot2::ggsave(paste0(output_path, "volcano_plot", suffix, ".", plot_device), plot = volcano_plot,
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  
  mess <- paste0(mess, "Volcano plot calculated. \n")
  
  
  
  #### Create Histogram for p-values and fold changes ####
  
  histograms <- pvalue_foldchange_histogram(RES = test_results, 
                                            columnname_p = "p", columnname_padj = "p.fdr", 
                                            columnname_FC = fc_col_name)
  
  ggplot2::ggsave(paste0(output_path, "histogram_p_value", suffix, ".", plot_device), plot = histograms[["histogram_p_value"]],
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  ggplot2::ggsave(paste0(output_path, "histogram_adjusted_p_value", suffix, ".", plot_device), plot = histograms[["histogram_adjusted_p_value"]],
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  ggplot2::ggsave(paste0(output_path, "histogram_fold_change", suffix, ".", plot_device), plot = histograms[["histogram_fold_change"]],
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  
  mess <- paste0(mess, "p-value, adjusted p-value and fold change histograms calculated. \n")
  
  
  #### Get significant candidates ####
  
  significance <- calculate_significance_categories_ttest(p = test_results[["p"]], 
                                                          p_adj = test_results[["p.fdr"]],
                                                          fc = test_results[[fc_col_name]])
  
  candidates <- as.character(significance)
  
  if(significant_after_FDR){
    candidates <- which(candidates == "significant after FDR correction")
  }else{
    candidates <- which(candidates == "significant" | candidates == "significant after FDR correction")
  }
  
  mess <- paste0(mess, "There are ", length(candidates), " candidates, which were significant", 
                 ifelse(significant_after_FDR, " after FDR correction. \n", ". \n"))
  
  
  
  #### Create Boxplots of Biomarker Candidates ####

  Boxplots_candidates(D = data[["D"]][candidates, ], 
                      protein.names = data[["ID"]][candidates, "protein"],
                      group = data[["group"]],
                      suffix = suffix,
                      output_path = paste0(output_path))
  
  mess <- paste0(mess, "Boxplots made for the candidates. \n")
  
  
  
  #### Create Heatmap ####
  
  set.seed(14) 
  
  t_heatmap <- Heatmap_with_groups(D = data[["D"]][candidates, ], 
                                   id = data[["ID"]][candidates, ],
                                   groups = data[["group"]])
  
  grDevices::pdf(paste0(output_path, "heatmap", suffix, ".pdf"), height = plot_height, width = plot_width)
  graphics::plot(t_heatmap[["heatmap"]])
  grDevices::dev.off()
  
  openxlsx::write.xlsx(cbind(data[["ID"]][candidates, ], zscore = t_heatmap[["data_as_matrix"]]), paste0(output_path, "heatmap_data", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)
  
  mess <- paste0(mess, "Heatmap made for the candidates. \n")
  
  
  
  #### Create On-Off Heatmap ####
  
  if(is.null(min_valid_values_on)){
    min_valid_values_on <- length(intensity_columns)
  }
  
  t_on_off_heatmap <- calculate_onoff(D = data[["D"]], 
                                      id = data[["ID"]],
                                      group = data[["group"]],
                                      max_vv_off = max_valid_values_off,
                                      min_vv_on = min_valid_values_on,
                                      protein_id_col = 1)
  
  grDevices::pdf(paste0(output_path, "on_off_heatmap", suffix, ".pdf"), height = plot_height, width = plot_width)
  graphics::plot(t_on_off_heatmap)
  grDevices::dev.off()
  
  mess <- paste0(mess, "On-Off-Heatmap made. \n", "There were ", sum(t_on_off_heatmap[["isonoff"]]), " on/off proteins.")
  
  
  
  #### Save message log ####
  
  cat(mess, file = paste0(output_path, "message_log_ttest", suffix, ".txt"))

  return(list("message" = mess))
}








#' The workflow for ANOVA of quantitative proteomics data
#'
#' @param data_path              \strong{character} \cr
#'                               The path to an .xlsx file containing the input data.
#' @param output_path            \strong{character} \cr
#'                               The path to the output folder.
#' @param intensity_columns      \strong{integer} \cr
#'                               A vector containing the numbers of the intensity columns in the table.
#' 
#' @param paired                 \strong{logical} \cr
#'                               If \code{TRUE}, a paired test will be done, otherwise an unpaired test.
#' @param var.equal              \strong{logical} \cr
#'                               If \code{TRUE}, the variances are assumed to be equal.
#' @param log_before_test        \strong{logical} \cr
#'                               If \code{TRUE}, the data will be log-transformed.
#' @param delog_for_FC           \strong{logical} \cr
#'                               If \code{TRUE}, the fold change will be calculated without the log-transformation.
#' 
#' @param significant_after_FDR  \strong{logical} \cr
#'                               If \code{TRUE}, candidates for the boxplots and heatmap need to be significant after FDR correction, otherwise all significant candidates will be used.
#' @param max_valid_values_off   \strong{integer} \cr
#'                               The maximum number of valid values to be an off protein.
#' @param min_valid_values_on    \strong{integer} \cr
#'                               The minimum number of valid values to be an on protein.
#' 
#' @param suffix                 \strong{character} \cr
#'                               The suffix of the file names should have one.
#' @param plot_device            \strong{character} \cr
#'                               The type of the output file, e.g. "pdf" or "png".
#' @param plot_height            \strong{numeric} \cr
#'                               The plot height in cm.
#' @param plot_width             \strong{numeric} \cr
#'                               The plot width in cm.
#' @param plot_dpi               \strong{integer} \cr
#'                               The "dots per inch" of the plot aka. the plot resolution.
#' 
#' 
#' @return Message log of the workflow
#' @export
#'
#' @seealso [workflow_ttest()] in case of more than two groups in the sample.\cr
#'          Functions used in this workflow: 
#'          [prepareTtestData()], [ttest()], [VolcanoPlot_ttest()], [pvalue_foldchange_histogram()], 
#'          [calculate_significance_categories_ttest()], [Boxplots_candidates()], 
#'          [Heatmap_with_groups()], [calculate_onoff()]
#'          
#'          
#' @examples
#' 
#' 
#' # 1. Set the character of your data path, leading to an .xlsx file.
#' in_path <- "C:/Users/thisuser/Documents/dataFolder/data.xlsx"
#' 
#' # 2. Set the integer vector of the columns, which contain the intensities.
#' int_col <- 3:17
#' 
#' # 3. Set the character of the output path, leading to a folder for the results.
#' out_path <- "C:/Users/thisuser/Documents/resultsFolder/"
#' 
#' # 4. Run the ANOVA with the parameters you set.
#' \dontrun{
#' result <- workflow_ANOVA(data_path = in_path,
#'                          output_path = out_path,
#'                          intensity_columns = int_col) }
#'

workflow_ANOVA <- function(data_path,
                           output_path,
                           intensity_columns,
                           
                           paired = FALSE,
                           var.equal = TRUE,
                           log_before_test = TRUE,
                           delog_for_FC = TRUE,
                           
                           significant_after_FDR = TRUE,
                           max_valid_values_off = 0,
                           min_valid_values_on = NULL,
                           
                           suffix = "",
                           plot_device = "pdf",
                           plot_height = 15,
                           plot_width = 15,
                           plot_dpi = 300
                           ){
  
  mess = ""
  
  
  #### Prepare Data ####
  
  data <- prepareTtestData(data_path = data_path , intensity_columns = intensity_columns)
  
  mess <- paste0(mess, "Data file contained ", length(levels(data[["group"]])), " groups with ", length(levels(data[["sample"]])), " samples. \n")
  
  
  
  #### Calculate ANOVA ####
  
  ANOVA_results <- ANOVA(D = data[["D"]], id = data[["ID"]], 
                         group = data[["group"]],  sample = data[["sample"]], 
                         paired = paired, var.equal = var.equal,
                         log_before_test = log_before_test, delog_for_FC = delog_for_FC, log_base = 2,
                         min_obs_per_group = 3, min_perc_per_group = NULL,
                         filename = paste0(output_path, "results_ANOVA", suffix, ".xlsx"))
  
  mess <- paste0(mess, "ANOVA calculated. \n")
  
  
  
  #### Create Volcano Plot ####
  
  p_posthoc_columns <- grep("p.posthoc.", colnames(ANOVA_results))
  fc_columns <- grep("FC_", colnames(ANOVA_results))
  
  # filter FC columns because FC is calculated "twice" e.g. for state1_divided_state2 and state1_divided_state2
  if(fc_columns[[1]]%%2 == 0){
    fc_columns <- fc_columns[fc_columns %% 2 == 0]
  }else{
    fc_columns <- fc_columns[fc_columns %% 2 != 0]
  }
  
  volcano_plots <- VolcanoPlot_ANOVA(RES = ANOVA_results, 
                                    columnname_p_ANOVA = "p.anova",
                                    columnname_p_ANOVA_adj = "p.anova.fdr",
                                    columns_FC = fc_columns,
                                    columns_p_posthoc = p_posthoc_columns )
  
  mess <- paste0(mess, length(volcano_plots), " volcano plots calculated. \n")
  
  grDevices::pdf(paste0(output_path, "volcano_plots", suffix ,".pdf"), height = plot_height, width = plot_width)
  for (v_plot in volcano_plots) {
    graphics::plot(x = v_plot)
  }
  grDevices::dev.off()
  
  rm(volcano_plots, v_plot)
  
  
  
  #### Create Histogram for p-values and fold changes ####
  
  histograms <- pvalue_foldchange_histogram(RES = ANOVA_results, 
                                            columnname_p = "p.anova", columnname_padj = "p.anova.fdr", 
                                            columnname_FC = colnames(ANOVA_results)[[fc_columns[[1]]]])
  
  ggplot2::ggsave(paste0(output_path, "histogram_p_value", suffix, ".", plot_device), plot = histograms[["histogram_p_value"]],
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  ggplot2::ggsave(paste0(output_path, "histogram_adjusted_p_value", suffix, ".", plot_device), plot = histograms[["histogram_adjusted_p_value"]],
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)
  
  mess <- paste0(mess, "p-value and adjusted p-value histograms calculated. \n")
  
  rm(histograms)
  
  
  
  #### Get significant candidates ####
  
  candidates <- list()
  
  for (i in 1:length(p_posthoc_columns)) {
    significance <- factor(calculate_significance_categories_ANOVA(p_anova = ANOVA_results[["p.anova"]],
                                                                                        p_anova_adj = ANOVA_results[["p.anova.fdr"]],
                                                                                        p_posthoc = ANOVA_results[[p_posthoc_columns[[i]]]],
                                                                                        fc = ANOVA_results[[fc_columns[[i]]]]))
    
    candidates[[i]] <- as.character(significance)
    
    if(significant_after_FDR){
      candidates[[i]] <- which(candidates[[i]] == "significant after FDR correction")
    }else{
      candidates[[i]] <- which(candidates[[i]] == "significant" | candidates[[i]] == "significant after FDR correction")
    }
    
    
  }
  
  union_candidates <- unique(unlist(candidates))
  
  mess <- paste0(mess, "There are ", length(union_candidates), " in the union of candidates, which were significant", 
                 ifelse(significant_after_FDR, " after FDR correction. \n", ". \n"))
  
  
  
  #### Create Boxplots of Biomarker Candidates ####
  
  Boxplots_candidates(D = data[["D"]][union_candidates, ], 
                      protein.names = data[["ID"]][union_candidates, "protein"],
                      group = data[["group"]],
                      suffix = suffix,
                      output_path = paste0(output_path))
  
  mess <- paste0(mess, "Boxplots made for the ", length(union_candidates) ," candidates. \n")
  
  
  
  #### Create Heatmap ####
  
  set.seed(27) #set seed so colors can be differentiated well: 14, 27, 101
  
  t_heatmap <- Heatmap_with_groups(D = data[["D"]][union_candidates, ], 
                                   id = data[["ID"]][union_candidates, ],
                                   groups = data[["group"]])
  
  grDevices::pdf(paste0(output_path, "heatmap", suffix, ".pdf"), height = plot_height, width = plot_width)
  graphics::plot(t_heatmap[["heatmap"]])
  grDevices::dev.off()
  
  remaining_candidates <- as.integer(row.names(t_heatmap[["data_as_matrix"]]))
  openxlsx::write.xlsx(cbind(data[["ID"]][remaining_candidates, ], zscore = t_heatmap[["data_as_matrix"]]), paste0(output_path, "heatmap_data", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)
  
  mess <- paste0(mess, "Heatmap made for the union of all candidates. \n")
  
  rm(t_heatmap, remaining_candidates)
  
  
  
  #### Create On-Off Heatmap ####
  
  if(is.null(min_valid_values_on)){
    min_valid_values_on <- length(intensity_columns)
  }
  
  t_on_off_heatmap <- calculate_onoff(D = data[["D"]], 
                                      id = data[["ID"]],
                                      group = data[["group"]],
                                      max_vv_off = max_valid_values_off,
                                      min_vv_on = min_valid_values_on,
                                      protein_id_col = 1)
  
  grDevices::pdf(paste0(output_path, "on_off_heatmap", suffix, ".pdf"), height = plot_height, width = plot_width)
  graphics::plot(t_on_off_heatmap)
  grDevices::dev.off()
  
  mess <- paste0(mess, "On-Off-Heatmap made. \n", "There were ", sum(t_on_off_heatmap[["isonoff"]]), " on/off proteins.")
  
  rm(t_on_off_heatmap)
  
  
  
  #### Save message log ####
  
  cat(mess, file = paste0(output_path, "message_log_anova", suffix, ".txt"))
  
  
  
  return(list("message" = mess))
}