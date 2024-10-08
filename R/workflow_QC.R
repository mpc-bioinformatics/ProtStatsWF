#' The main workflow for quality control of quantitative proteomics data
#'
# path parameters
#' @param data_path         A character containing the path to an .xlsx file.
#' @param output_path       A character containing the path to an output folder.
#'
# mandatory parameters
#' @param intensity_columns An integer vector containing the intensity columns of the table.
#' @param normalization_method A character containing the method of normalization. The possible methods are no normalization "nonorm" or "median", "loess", "quantile" or "lts" normalization.
#' @param lts_quantile          A numeric containing the quantile for the lts normalization if \code{normalization = "lts"}, default is 0.8.
#' @param use_groups    If \code{TRUE}, group information encoded in the column names are used. Default is \code{TRUE}.
#'
# additional parameters
#' @param na_strings A character vector containing symbols to be recognized as missing values (with the exception of 0).
#' @param zero_to_NA If \code{TRUE}, 0 will be treated as missing value.
#' @param do_log_transformation If \code{TRUE}, the data will be log-transformed.
#' @param log_base              A numeric containing the base used, if data is log-transformed.
#' @param groupvar_name A character containing the name for the group variable.
#' @param group_colours A character vector of hex codes for the group colors, if the data has groups. If \code{NULL}, a default color scale will be used.
### TODO: ideally, group_colours would be a named vector!
#'
#' @param suffix      A character containing the suffix for the output files. Needs to start with an underscore.
#'
# general plot parameters
#' @param base_size   A numeric containing the base size of the font.
#' @param plot_device A character containing the type of the output file, e.g. "pdf" or "png".
#' @param plot_height A numeric of the plot height in cm.
#' @param plot_width  A numeric of the plot width in cm.
#' @param plot_dpi    A numeric of the "dots per inch" of the plot aka. the plot resolution.
#'
#'
# Boxplot parameters
#' @param boxplot_method A character containing the method used. Possible are "boxplot" and "violinplot". Default is "boxplot".
#'
# MA-Plot parameters
#' @param generate_MAplots If \code{TRUE}, MA plots will be generated, if \code{FALSE} they will not be generated (mostly for debugging purposes).
#' @param MA_maxPlots A numeric containing the maximum number of MA plots that should be generated. Defaults is 5000.
#' @param MA_alpha    If \code{TRUE}, the data points of the MA plots will be transparent.
#' @param MA_sampling A numeric containing the sampling rate for MA-Plots. Useful to sample part of the data set for data sets on peptide/feature level with many data points.
#'
# PCA parameters:
#' @param PCA_impute         If \code{TRUE}, missing values will be imputed.
#' @param PCA_impute_method  A character containing the imputation method ("mean" or "median")
#' @param PCA_propNA         A numeric of the proportion of allowed missing NAs for a protein, before it is discarded.
#' @param PCA_scale.         If \code{TRUE}, the data will be scaled before computing the PCA.
#' @param PCA_PCx            The principle component for the x-axis (default: 1).
#' @param PCA_PCy            The principle component for the y-axis (default: 2).
#' @param PCA_groupvar1_name Titles of legends for colour.
#' @param PCA_alpha          If \code{TRUE}, the data points of the PCA plot will be transparent.
#' @param PCA_label          If \code{TRUE}, the samples will be labeled.
#' @param PCA_label_seed     A numeric, which sets the seed for the label.
#' @param PCA_label_size     A numeric containing the size of the sample labels.
#' @param PCA_xlim           Limit of the x-axis.
#' @param PCA_ylim           Limit of the y-axis.
#' @param PCA_point.size     The size of the data points.
#'
#'
#'
#'
#'
#'
#' @return saves several plots and excel files
#' @export
#'
#' @examples
#' \dontrun{
#' in_path <- "/Users/thisuser/Documents/dataFolder/data.xlsx"
#' intensity_cols <- 3:17
#' out_path <- "/Users/thisuser/Documents/resultsFolder/"
#'
#' result <- workflow_QC(data_path = in_path,
#'                       intensity_columns = intensity_cols,
#'                       output_path = out_path)
#'}
#'



workflow_QC <- function(data_path,
                        output_path,

                        intensity_columns,
                        normalization_method = "loess",
                        lts_quantile = 0.8,
                        use_groups = TRUE,

                        na_strings = c("NA", "NaN", "Filtered","#NV"),
                        zero_to_NA = TRUE,
                        do_log_transformation = TRUE,
                        log_base = 2,

                        groupvar_name = "Group",
                        group_colours = NULL,

                        base_size = 15,
                        plot_device = "pdf",
                        plot_height = 10,
                        plot_width = 15,
                        plot_dpi = 300,

                        suffix = "_",


                        boxplot_method = "boxplot",

                        generate_MAplots = TRUE,
                        MA_maxPlots = 5000,
                        MA_alpha = FALSE,
                        MA_sampling = 1,

                        #PCA_groupvar1 = "group",
                        #PCA_groupvar2 = NULL,

                        PCA_impute = FALSE, PCA_impute_method = "mean", PCA_propNA = 0,
                        PCA_scale. = TRUE,
                        PCA_PCx = 1, PCA_PCy = 2,
                        PCA_groupvar1_name = "group",
                        #PCA_groupvar2_name = NULL,
                        PCA_alpha = 1, PCA_label = FALSE, PCA_label_seed = NA, PCA_label_size = 4,
                        PCA_xlim = NULL, PCA_ylim = NULL, PCA_point.size = 4

                        ){

  mess = ""


  #### Prepare Data ####

  prepared_data <- prepareData(data_path = data_path, intensity_columns = intensity_columns,
                               na_strings = na_strings, zero_to_NA = zero_to_NA,
                               do_log_transformation = do_log_transformation, log_base = log_base,
                               use_groups = use_groups, group_colours = group_colours,
                               normalization = normalization_method, lts_quantile = lts_quantile)


  mess <- paste0(mess, prepared_data[["message"]])#

  group <- prepared_data$group

  write.csv(x = prepared_data$ID, file = paste0(output_path, "/ID", suffix, ".csv"), row.names = FALSE)
  write.csv(x = prepared_data$D, file = paste0(output_path, "/D_norm_wide", suffix, ".csv"), row.names = FALSE)
  write.csv(x = prepared_data$D_long, file = paste0(output_path, "/D_norm_long", suffix, ".csv"), row.names = FALSE)

  openxlsx::write.xlsx(x = cbind(prepared_data$ID, prepared_data$D), file = paste0(output_path, "/D_norm_ID", suffix, ".xlsx"),
                       rowNames = FALSE, overwrite = TRUE, keepNA = TRUE)


  #### Calculate Valid Value Plot ####

  vv_plot_data <- ValidValuePlot(D_long = prepared_data[["D_long"]],
                                          use_groups = use_groups, groupvar_name = groupvar_name, group_colours = group_colours,
                                          base_size = base_size)
  #### reorder valid values table to stay in the same order as the original data ####
  cnames <- colnames(prepared_data$D)
  vv_plot_data$table$name <- factor(vv_plot_data$table$name, levels = cnames)
  vv_plot_data$table <- vv_plot_data$table[order(vv_plot_data$table$name),]

  mess <- paste0(mess, vv_plot_data[["message"]])


  ggplot2::ggsave(paste0(output_path, "/valid_value_plot", suffix, ".", plot_device), plot = vv_plot_data[["plot"]],
         device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi, units = "cm")
  write.csv(x = vv_plot_data$table, file = paste0(output_path, "/D_validvalues", suffix, ".csv"), row.names = FALSE)



  #### Calculate Valid Value Plot ####

  boxplot_data <- Boxplots(D_long = prepared_data[["D_long"]],
                           do_log_transformation = FALSE, log_base = log_base,
                           use_groups = use_groups, groupvar_name = groupvar_name, group_colours = group_colours,

                           base_size = base_size, method = boxplot_method)

  mess <- paste0(mess, boxplot_data[["message"]])

  ggplot2::ggsave(paste0(output_path, "/boxplot", suffix, ".", plot_device), plot = boxplot_data[["plot"]],

         device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi, units = "cm")


  #### Calculate MA Plot ####

  if (generate_MAplots) {
  ma_data <- MA_Plots(D = prepared_data[["D"]],
                      do_log_transformation = FALSE,
                      output_path = output_path, suffix = suffix,
                      labels = 1:ncol(prepared_data[["D"]]), labels2 = colnames(prepared_data[["D"]]),
                      maxPlots = MA_maxPlots, alpha = MA_alpha,
                      plot_height = plot_height, plot_width = plot_width, sampling = MA_sampling)

  mess <- paste0(mess, ma_data)
  }


  #### Calculate PCA Plot ####


  ### depending of groups should be used or not, the group variable is used in the PCA function
  if (use_groups) {
    PCA_groupvar1 <- group
  } else {
    PCA_groupvar1 <- NULL
  }


  pca_data <- PCA_Plot(D = prepared_data[["D"]],
                       groupvar1 = group,
                       groupvar2 = NULL,
                       impute = PCA_impute, impute_method = PCA_impute_method, propNA = PCA_propNA,
                       scale. = PCA_scale.,
                       PCx = PCA_PCx, PCy = PCA_PCy,
                       groupvar1_name = PCA_groupvar1_name,
                       groupvar2_name = NULL,
                       group_colours = group_colours, PCA_alpha = 1,
                       label = PCA_label, PCA_label_seed = NA, PCA_label_size = 4,
                       xlim = PCA_xlim, ylim = PCA_ylim,
                       point.size = PCA_point.size, base_size = base_size)

  mess <- paste0(mess, pca_data[["message"]])


  ggplot2::ggsave(paste0(output_path, "/PCA_plot", suffix, ".", plot_device), plot = pca_data[["plot"]],
         device = plot_device, height = plot_height, width = plot_width, dpi = plot_dpi, units = "cm")
  write.csv(x = pca_data$D_PCA_plot, file = paste0(output_path, "/D_PCA", suffix, ".csv"), row.names = FALSE)


  return (list("message" = mess))
}
