#' QC Workflow (QCQuant)
#'
#' @description
#' Workflow for quality control of quantitative proteomics data
#'
#' @details
#' This function performs quality control of quantitative proteomics data.
#' As a first step, the data can be normalized.
#' Then the following plots are generated: a valid value plot, boxplots, MA-plots and a PCA plot.
#'
#'
#'
#' @param data_path              \strong{character} \cr
#'                               The path to an .xlsx file containing the input data.
#'
#'
#' @param filetype **character(1)** \cr Type of input file: "csv" or "tsv" or "txt" or "xlsx".
#' @param sep **character(1)** \cr The field separator, e.g. " " for blanks, "," for comma or "\\t" for tab. Default is ",".
#' @param dec **character(1)** \cr Decimal separator, e.g. "," for comma or "." for dot. Default is ".".
#' @param header **logical(1)** \cr If TRUE, first line is counted as column names.
#' @param sheet **integer(1)** \cr Sheet number (only needed for xlsx files, default is to use the first sheet).
#' @param output_path            \strong{character} \cr
#'                               The path to the output folder.
#' @param output_type **character(1)** \cr Type of input file: "csv" or "tsv" or "xlsx".
#'
# mandatory parameters
#' @param intensity_columns      \strong{integer vector} \cr
#'                               The numbers of the intensity columns in the table.
#' @param normalization_method   \strong{character} \cr
#'                               The method of normalization. Options are "nonorm" (no normalization), "median", "loess", "quantile" or "lts" normalization.
#' @param lts_quantile           \strong{numeric} \cr
#'                               The quantile for the lts normalization if \code{normalization = "lts"}.
#' @param use_groups             \strong{logical} \cr
#'                               If \code{TRUE}, group information encoded in the column names is used.
#'
# additional parameters
#' @param na_strings             \strong{character} \cr
#'                               A vector containing the symbols to be recognized as missing values (with the exception of 0).
#' @param na_out                 \strong{character} \cr
#'                               NA values will be converted to this character before writing results.
#' @param zero_to_NA             \strong{logical} \cr
#'                               If \code{TRUE}, 0 will be treated as missing value.
#' @param do_log_transformation  \strong{logical} \cr
#'                               If \code{TRUE}, the data will be log-transformed.
#' @param log_base               \strong{numeric} \cr
#'                               The base used, if \code{do_log_transformation = TRUE}.
#' @param groupvar_name          \strong{character} \cr
#'                               The name for the group variable.
#' @param group_colours          \strong{character vector} \cr
#'                               The hex codes for the group colors, if the data has groups. If \code{NULL}, a default color scale will be used.
#'
#' @param suffix                 \strong{character} \cr
#'                               The suffix for the output files. It needs to start with an underscore.
#'
# general plot parameters
#' @param plot                   \strong{logical} \cr
#'                               If TRUE (default), plots are generated. If FALSE, plots are skipped (may be useful for debugging purposes and reduces runtime for large data sets).
#' @param base_size              \strong{numeric} \cr
#'                               The base size of the font.
#' @param plot_device            \strong{character} \cr
#'                               The type of the output file. Options are "pdf" or "png".
#' @param plot_height            \strong{numeric} \cr
#'                               The plot height in cm.
#' @param plot_width             \strong{numeric} \cr
#'                               The plot width in cm.
#' @param plot_dpi               \strong{numeric} \cr
#'                               The "dots per inch" of the plot aka. the plot resolution.
#'
#'
# Boxplot parameters
#' @param boxplot_method          \strong{character} \cr
#'                                The method used. Options are "boxplot" and "violinplot".
#'
# MA-Plot parameters
#' @param generate_MAplots       \strong{logical} \cr
#'                               If \code{TRUE}, MA plots will be generated, if \code{FALSE} they will not be generated (mostly for debugging purposes).
#' @param MA_maxPlots            \strong{integer} \cr
#'                               The maximum number of MA plots that should be generated.
#' @param MA_alpha               \strong{logical} \cr
#'                               If \code{TRUE}, the data points of the MA plots will be transparent.
#' @param MA_sampling            \strong{numeric} \cr
#'                               The sampling rate for MA-Plots. Useful to sample part of the data set for data sets on peptide/feature level with many data points.
#'
# PCA parameters:
#' @param PCA_impute             \strong{logical} \cr
#'                               If \code{TRUE}, missing values will be imputed for the PCA.
#' @param PCA_impute_method      \strong{character} \cr
#'                               The imputation method. Options are "mean" or "median".
#' @param PCA_propNA             \strong{numeric} \cr
#'                               The proportion of allowed missing NAs for a protein, before it is discarded.
#' @param PCA_scale.             \strong{logical} \cr
#'                               If \code{TRUE}, the data will be scaled before computing the PCA.
#' @param PCA_PCx                \strong{numeric} \cr
#'                               The principle component for the x-axis.
#' @param PCA_PCy                \strong{numeric} \cr
#'                               The principle component for the y-axis.
#' @param PCA_groupvar1_name     \strong{character} \cr
#'                               The titles of legends for colour.
#' @param PCA_alpha              \strong{logical} \cr
#'                               If \code{TRUE}, the data points of the PCA plot will be transparent.
#' @param PCA_label              \strong{logical} \cr
#'                               If \code{TRUE}, the samples will be labeled.
#' @param PCA_label_seed         \strong{numeric} \cr
#'                               The seed for the label.
#' @param PCA_label_size         \strong{numeric} \cr
#'                               The size of the sample labels.
#' @param PCA_xlim               \strong{numeric} \cr
#'                               The limit of the x-axis.
#' @param PCA_ylim               \strong{numeric} \cr
#'                               The limit of the y-axis.
#' @param PCA_point.size         \strong{numeric} \cr
#'                               The size of the data points.
#'
#'
#'
#'
#'
#'
#' @return The workflow saves several plots and excel files and returns a message log of the workflow.
#' @export
#'
#' @seealso Functions used in this workflow:
#'          [prepareData()], [ValidValuePlot()], [Boxplots()], [MA_Plots()], [PCA_Plot()]
#'
#' @examples
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
#' # 4. Run the QC workflow with the parameters you set.
#' \dontrun{
#' result <- workflow_QC(data_path = in_path,
#'                          output_path = out_path,
#'                          intensity_columns = int_col) }
#'



workflow_QC <- function(data_path,
                        filetype = "xlsx",
                        sep = ",",
                        dec = ".",
                        header = TRUE,
                        sheet = 1,
                        output_path,
                        output_type = "xlsx",

                        intensity_columns,
                        normalization_method = "loess",
                        lts_quantile = 0.8,
                        use_groups = TRUE,

                        na_strings = c("NA", "NaN", "Filtered","#NV"),
                        na_out = "NA",
                        zero_to_NA = TRUE,
                        do_log_transformation = TRUE,
                        log_base = 2,

                        groupvar_name = "Group",
                        group_colours = NULL,

                        suffix = "_",

                        plot = TRUE,
                        base_size = 15,
                        plot_device = "pdf",
                        plot_height = 10,
                        plot_width = 15,
                        plot_dpi = 300,

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

  prepared_data <- prepareData(data_path = data_path,
                               filetype = filetype,
                               sep = sep,
                               dec = dec,
                               header = header,
                               sheet = sheet,
                               intensity_columns = intensity_columns,
                               na_strings = na_strings, zero_to_NA = zero_to_NA,
                               do_log_transformation = do_log_transformation, log_base = log_base,
                               use_groups = use_groups, group_colours = group_colours,
                               normalization = normalization_method, lts_quantile = lts_quantile)


  mess <- paste0(mess, prepared_data[["message"]])

  group <- prepared_data$group

  utils::write.csv(x = prepared_data$ID, file = paste0(output_path, "/ID", suffix, ".csv"), row.names = FALSE)
  utils::write.csv(x = prepared_data$D, file = paste0(output_path, "/D_norm_wide", suffix, ".csv"), row.names = FALSE)
  utils::write.csv(x = prepared_data$D_long, file = paste0(output_path, "/D_norm_long", suffix, ".csv"), row.names = FALSE)


  if (output_type == "xlsx") {
    openxlsx::write.xlsx(x = cbind(prepared_data$ID, prepared_data$D), file = paste0(output_path, "/D_norm_ID", suffix, ".xlsx"),
                         rowNames = FALSE, overwrite = TRUE, keepNA = TRUE, na.string = na_out)
  }
  if (output_type == "csv") {
    utils::write.csv(x = cbind(prepared_data$ID, prepared_data$D), file = paste0(output_path, "/D_norm_ID", suffix, ".csv"),
                     row.names = FALSE, na = na_out)
  }
  if (output_type == "tsv") {
    utils::write.table(x = cbind(prepared_data$ID, prepared_data$D), file = paste0(output_path, "/D_norm_ID", suffix, ".tsv"),
                       row.names = FALSE, sep = "\t", na = na_out)
  }


  if (plot) {
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
    utils::write.csv(x = vv_plot_data$table, file = paste0(output_path, "/D_validvalues", suffix, ".csv"), row.names = FALSE)



    #### Calculate Boxlots ####

    boxplot_data <- Boxplots(D_long = prepared_data[["D_long"]],
                             do_log_transformation = FALSE, log_base = log_base,
                             use_groups = use_groups, groupvar_name = groupvar_name, group_colours = group_colours,
                             base_size = base_size, method = boxplot_method, lwd = 0.5)

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
    utils::write.csv(x = pca_data$D_PCA_plot, file = paste0(output_path, "/D_PCA", suffix, ".csv"), row.names = FALSE)
    utils::write.csv(x = pca_data$filtered_D, file = paste0(output_path, "/PCA_data_after_imputation", suffix, ".csv"), row.names = FALSE)


  }

  return(list("message" = mess))
}
