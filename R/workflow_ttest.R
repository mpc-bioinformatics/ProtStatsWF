#' t-test workflow
#'
#' @description
#' Workflow for t-test analysis of quantitative proteomics data.
#'
#' @details
#' This function performs a t-test to compare two experimental groups in a quantitative proteomics dataset.
#' Ideally, the data should be already normalized before performing the t-test (e.g. by using the [workflow_QC()] function).
#' It is recommended to log-transform the data before the t-test (either use already log-transformed data and \code{log_before_test = FALSE}
#' or log-transform the data on the fly with \code{log_before_test = TRUE}).
#' The function generates a volcano plot, a histogram of p-values and fold changes,
#' boxplots of the significant candidates, and a heatmap of the significant candidates.
#'
#'
#'
#' @param data_path              \strong{character} \cr
#'                               The path to an .xlsx file containing the input data.
#' @param output_path            \strong{character} \cr
#'                               The path to the output folder.
#' @param intensity_columns      \strong{integer vector} \cr
#'                               The numbers of the intensity columns in the table.
#'
#' @param paired                 \strong{logical} \cr
#'                               If \code{TRUE}, a paired test will be done, otherwise an unpaired test.
#' @param var.equal              \strong{logical} \cr
#'                               If \code{TRUE}, the variances are assumed to be equal.
#' @param log_before_test        \strong{logical} \cr
#'                               If \code{TRUE}, the data will be log-transformed.
#' @param delog_for_FC           \strong{logical} \cr
#'                               If \code{TRUE}, the fold change will be calculated without the log-transformation.
#' @param p_value_zeros_to_min   \strong{logical} \cr
#'                               If \code{TRUE}, then \code{p_values == 0} will be set to the next smallest value of the p-values.
#'
#' @param volcano_base_size      \strong{numeric} \cr
#'                               The base size of the volcano plot.
#'
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
#' @param column_name_protein    \strong{character} \cr
#'                               The column name containing the proteins.
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
                           p_value_zeros_to_min = TRUE,

                           volcano_base_size = 25,

                           significant_after_FDR = TRUE,
                           max_valid_values_off = 0,
                           min_valid_values_on = NULL,

                           suffix = "",
                           plot_device = "pdf",
                           plot_height = 15,
                           plot_width = 15,
                           plot_dpi = 300,

                           column_name_protein = "Protein"
                           ) {

  mess = ""


  #### Prepare Data ####

  data <- prepareTtestData(data_path = data_path , intensity_columns = intensity_columns)



  #### Calculate ttest ####

  test_results <<- ttest(D = data[["D"]], id = data[["ID"]],
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


  if(p_value_zeros_to_min){

    p_value_zero <- which(test_results$p == 0)

    if(length(p_value_zero) > 0){
      next_smallest_value <- sort(unique(test_results$p))[2]

      test_results$p[test_results$p == 0] <- next_smallest_value

      mess <- paste0(mess, "There were ", length(p_value_zero), " p_values, which were 0. ",
                     "They were set to the next smallest occuring value ", next_smallest_value, ". \n")
    }
  }


  fc_col_name <- paste0("FC_", levels(data[["group"]])[[1]], "_divided_by_", levels(data[["group"]])[[2]])



  #### Create Volcano Plot ####

  volcano_plot <- VolcanoPlot_ttest(RES = test_results,
                                    columnname_p = "p", columnname_padj = "p.fdr",
                                    columnname_FC = fc_col_name, base_size = volcano_base_size)

  ggplot2::ggsave(paste0(output_path, "volcano_plot", suffix, ".", plot_device), plot = volcano_plot,
                  device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi)

  mess <- paste0(mess, "Volcano plot calculated. \n")



  #### Create Histogram for p-values and fold changes ####

  histograms <- ProtStatsWF::pvalue_foldchange_histogram(RES = test_results,
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

  significance <- ProtStatsWF::calculate_significance_categories_ttest(p = test_results[["p"]],
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

  ProtStatsWF::Boxplots_candidates(D = data[["D"]][candidates, ],
                      protein.names = data[["ID"]][candidates, column_name_protein],
                      group = data[["group"]],
                      suffix = suffix,
                      output_path = paste0(output_path),
                      log_data = log_before_test)

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



