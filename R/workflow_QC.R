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
#' @param dataPath              \strong{character} \cr
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
#' @param base_size              \strong{numeric} \cr
#'                               The base size of the font.
#' @param plot_device            \strong{character} \cr
#'                               The type of the output file. Options are "pdf" or "png".
#' @param plot_height_BP_VV      \strong{numeric} \cr
#'                               The plot height for boxplots and valid-value plots in cm.
#' @param plot_width_BP_VV       \strong{numeric} \cr
#'                               The plot width for boxplots and valid-value plots in cm. The optimal width is highly affected by the number of samples.
#' @param plot_height_PCA_MA     \strong{numeric} \cr
#'                               The plot height for PCA and MA-plots in cm.
#' @param plot_width_PCA_MA      \strong{numeric} \cr
#'                               The plot width for PCA and MA-plots in cm.
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
#'                               If \code{TRUE}, MA plots will be generated, if \code{FALSE} they will not be generated (mostly for debugging purposes and time reduction in large datasets).
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



workflow_QC <- function(dataPath,
                        intensityColumns,
                        proteinNameColumn = "Protein",
                        sampleInfoPath = NULL,
                        sampleNameColumn = "sampleName",
                        groupColumn = NULL,
                        group2Column = NULL,

                        fileType = "xlsx",
                        sep = ",",
                        dec = ".",
                        header = TRUE,
                        sheet = 1,
                        output_path,
                        output_type = "xlsx",


                        normMethod = "loess",
                        ltsQuantile = 0.8,

                        NAStrings = c("NA", "NaN", "Filtered","#NV"),
                        na_out = "NA",
                        zeroToNA = TRUE,
                        doLogTrans = TRUE,
                        logBase = 2,

                        groupvar_name = "Group",
                        group_colours = NULL,

                        suffix = "",
                        verbose = TRUE,

                        base_size = 15,
                        plot_device = "pdf",
                        plot_height_BP_VV = 10,
                        plot_width_BP_VV = 15,
                        plot_height_PCA_MA = 15,
                        plot_width_PCA_MA = 15,
                        plot_dpi = 300,

                        boxplot_method = "boxplot",

                        generate_MAplots = TRUE,
                        MA_maxPlots = 5000,
                        MA_alpha = FALSE,
                        MA_sampling = 1,

                        PCA_impute = FALSE,
                        PCA_impute_method = "mean",
                        PCA_propNA = 0,
                        PCA_scale. = TRUE,
                        PCA_PCx = 1,
                        PCA_PCy = 2,
                        PCA_groupvar1_name = "group",
                        PCA_alpha = 1,
                        PCA_label = FALSE,
                        PCA_label_seed = NA,
                        PCA_label_size = 4,
                        PCA_xlim = NULL,
                        PCA_ylim = NULL,
                        PCA_point.size = 4

){



  #### Prepare Data ####

  prepared_data <- prepareDataSE(dataPath,
                                 intensityColumns,
                                 proteinNameColumn = proteinNameColumn,
                                 sampleInfoPath = sampleInfoPath,
                                 sampleNameColumn = sampleNameColumn,
                                 doLogTrans = doLogTrans,
                                 logBase = logBase,
                                 normMethod = normMethod,
                                 ltsQuantile = ltsQuantile,

                                 fileType = fileType,
                                 sep = sep,
                                 dec = dec,
                                 header = header,
                                 sheet = sheet,
                                 zeroToNA = zeroToNA,
                                 NAStrings = NAStrings,
                                 verbose = verbose)

  if (output_type == "xlsx") {
    exportSE(prepared_data$SE, file = file.path(output_path, paste0("D_norm", suffix, ".xlsx")))
  }
           
  # prepare group colours
  group <- summarizedExperiment::colData(prepared_data$SE)[, groupColumn]
  nr_groups <- length(levels(group))
  if (is.null(group_colours) & nr_groups >= 1) group_colours <- scales::hue_pal()(nr_groups)

  #### Calculate Valid Value Plot ####
  vv_plot <- ValidValuePlot(D_long = prepared_data$D_long,
                                 groupColumn = groupColumn,
                                 group_colours = group_colours,
                                 base_size = base_size)

  #### TODO: reorder valid values table to stay in the same order as the original data ####
  #cnames <- colnames(prepared_data$D)
  #vv_plot_data$table$name <- factor(vv_plot_data$table$name, levels = cnames)
  #vv_plot_data$table <- vv_plot_data$table[order(vv_plot_data$table$name),]

  ggplot2::ggsave(file.path(output_path, paste0("valid_value_plot", suffix, ".", plot_device)), plot = vv_plot$plot,
                  device = plot_device, height = plot_height_BP_VV, width = plot_width_BP_VV, dpi = plot_dpi, units = "cm")
  utils::write.csv(x = vv_plot$table, file = file.path(output_path, paste0("D_validvalues", suffix, ".csv")), row.names = FALSE)



  #### Calculate Boxlots ####

  boxplots <- Boxplots(D_long = prepared_data[["D_long"]],
                           groupColumn = groupColumn,
                           group_colours = group_colours,
                           base_size = base_size, method = boxplot_method, lwd = 0.5)

  ggplot2::ggsave(file.path(output_path, paste0("boxplot", suffix, ".", plot_device)), plot = boxplots,
                  device = plot_device, height = plot_height_BP_VV, width = plot_width_BP_VV, dpi = plot_dpi, units = "cm")





  #### Calculate MA Plot ####



  if (generate_MAplots) {
    ma_data <- MA_Plots(D = prepared_data$SE,
                        output_path = output_path, suffix = suffix,
                        labels = 1:ncol(prepared_data[["D"]]), labels2 = colnames(prepared_data[["D"]]),
                        maxPlots = MA_maxPlots, alpha = MA_alpha,
                        plot_height = plot_height_PCA_MA, plot_width = plot_width_PCA_MA, sampling = MA_sampling)
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


  ggplot2::ggsave(file.path(output_path, paste0("PCA_plot", suffix, ".", plot_device)), plot = pca_data[["plot"]],
                  device = plot_device, height = plot_height_PCA_MA, width = plot_width_PCA_MA, dpi = plot_dpi, units = "cm")
  utils::write.csv(x = pca_data$D_PCA_plot, file = file.path(output_path, paste0("D_PCA", suffix, ".csv")), row.names = FALSE)
  utils::write.csv(x = pca_data$filtered_D, file = file.path(output_path, paste0("PCA_data_after_imputation", suffix, ".csv")), row.names = FALSE)


  return(invisible(NULL))
}
