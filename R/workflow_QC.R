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
#' @param D              \strong{list} \cr
#'                               Output from [prepareDataSE()], containing the data
#'                               as a SummarizedExperiment object and in long format.
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
# additional parameters
#' @param na_out                 \strong{character} \cr
#'                               NA values will be converted to this character before writing results.
#' @param groupName              \strong{character} \cr
#'                               The name for the first group variable (used for colour in plots).
#' @param group2Name              \strong{character} \cr
#'                               The name for the second group variable (used for shape in the PCA-plot).
#' @param group_colours          \strong{character vector} \cr
#'                               The hex codes for the group colors, if the data has groups. If \code{NULL}, the default color scale from ggplot2 will be used.
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
#'                               The plot resolution in "dots per inch".
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
#' @param PCA_scale              \strong{logical} \cr
#'                               If \code{TRUE}, the data will be scaled before computing the PCA.
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



workflow_QC <- function(D,
                        #intensityColumns,
                        #proteinNameColumn = "Protein",
                        #sampleInfoPath = NULL,
                        #sampleNameColumn = "sampleName",
                        groupColumn = NULL,
                        group2Column = NULL,
                        groupColours = NULL,

                        #fileType = "xlsx",
                        #sep = ",",
                        #dec = ".",
                        #header = TRUE,
                        #sheet = 1,

                        #NAStrings = c("NA", "NaN", "Filtered","#NV"),
                        #zeroToNA = TRUE,
                        #doLogTrans = TRUE,
                        #logBase = 2,
                        #normMethod = "loess",
                        #ltsQuantile = 0.8,

                        outPath,
                        outType = "xlsx",
                        suffix = "",
                        NAOut = "NA",
                        verbose = TRUE,

                        baseSize = 15,
                        plotDevice = "pdf",
                        plotHeight_BP_VV = 10,
                        plotWidth_BP_VV = 15,
                        plotHeight_PCA_MA = 15,
                        plotWidth_PCA_MA = 15,
                        plotDPI = 300,

                        boxplotMethod = "boxplot",

                        MAMaxPlots = 5000, # TODO: if = 0, MA plost are skipped
                        MAAlpha = 1, # TODO: give degree of alpha, not only TRUE/FALSE

                        PCAImputeMethod = "mean", ## TODO: none means no imputation
                        PCAPropNA = 0,
                        PCAScale = TRUE,
                        PCAAlpha = 1,
                        PCALabel = FALSE,
                        PCALabelSeed = NA,
                        PCALabelSize = 4,
                        PCAXlim = NULL,
                        PCAYlim = NULL,
                        PCAPointSize = 4
){

# TODO: check all variables that are not used in another main function in this workflow

  #### Prepare Data ####

  # D <- prepareDataSE(dataPath,
  #                                intensityColumns,
  #                                proteinNameColumn = proteinNameColumn,
  #                                sampleInfoPath = sampleInfoPath,
  #                                sampleNameColumn = sampleNameColumn,
  #                                doLogTrans = doLogTrans,
  #                                logBase = logBase,
  #                                normMethod = normMethod,
  #                                ltsQuantile = ltsQuantile,
  #
  #                                fileType = fileType,
  #                                sep = sep,
  #                                dec = dec,
  #                                header = header,
  #                                sheet = sheet,
  #                                zeroToNA = zeroToNA,
  #                                NAStrings = NAStrings,
  #                                verbose = verbose)


  ### TODO: na out and other output types

  # prepare group colours
  group <- SummarizedExperiment::colData(D$SE)[, groupColumn]
  nr_groups <- length(levels(group))
  if (is.null(groupColours) & nr_groups >= 1) groupColours <- scales::hue_pal()(nr_groups)

  #### Calculate Valid Value Plot ####
  vv_plot <- ValidValuePlot(D_long = D$D_long,
                                 groupColumn = groupColumn,
                                 groupColours = groupColours,
                                 baseSize = baseSize)

  ggplot2::ggsave(file.path(outPath, paste0("valid_value_plot", suffix, ".", plotDevice)),
                  plot = vv_plot$plot, device = plotDevice, height = plotHeight_BP_VV,
                  width = plotWidth_BP_VV, dpi = plotDPI, units = "cm")

  ### TODO: different output file types
  utils::write.csv(x = vv_plot$table, file = file.path(outPath,
                    paste0("D_validvalues", suffix, ".csv")), row.names = FALSE)



  #### Calculate Boxlots ####

  boxplots <- Boxplots(D_long = D$D_long,
                           groupColumn = groupColumn,
                           groupColours = groupColours,
                           baseSize = baseSize,
                           method = boxplotMethod, lwd = 0.5)

  ggplot2::ggsave(file.path(outPath, paste0("boxplot", suffix, ".", plotDevice)),
                  plot = boxplots, device = plotDevice, height = plotHeight_BP_VV,
                  width = plotWidth_BP_VV, dpi = plotDPI, units = "cm")





  #### Calculate MA Plot ####



  if (MAMaxPlots > 0) {
    ma_data <- MA_Plots(D = D$SE,
                        outPath = outPath, suffix = suffix,
                        labels = 1:ncol(D[["D"]]),  # TODO
                        labels2 = colnames(D[["D"]]),  # TODO
                        maxPlots = MAMaxPlots, alpha = MAAlpha,
                        plot_height = plotHeight_PCA_MA,
                        plot_width = plotWidth_PCA_MA,
                        sampling = 1, verbose = verbose)
  }


  #### Calculate PCA Plot ####


  pca_data <- PCA_Plot(D = preparedData$SE,
                       groupForColour = groupColumn,
                       groupForShape = group2Column,
                       #impute = PCA_impute,
                       imputeMethod = PCAImputeMethod,
                       propNA = PCAPropNA,
                       scale. = PCAScale,
                       PCx = 1, PCy = 2,
                       groupvar1_name = groupColumn,
                       groupvar2_name = group2Column,
                       groupColours = group_colours,
                       alpha = PCAAlpha,
                       label = PCALabel, PCALabelSeed = NA, PCALabelSize = 4,
                       xlim = PCAXlim, ylim = PCAYlim,
                       pointSize = PCAPointSize, baseSize = baseSize)

### TODO:: extract and save loadings
  # Loadings <- as.data.frame(pca$rotation)
  # if (!is.null(id)) {
  #   Loadings <- cbind(id, Loadings)


  ggplot2::ggsave(file.path(outPath, paste0("PCA_plot", suffix, ".", plotDevice)), plot = pca_data[["plot"]],
                  device = plotDevice, height = plot_height_PCA_MA, width = plot_width_PCA_MA, dpi = plot_dpi, units = "cm")
  utils::write.csv(x = pca_data$D_PCA_plot, file = file.path(outPath, paste0("D_PCA", suffix, ".csv")), row.names = FALSE)
  utils::write.csv(x = pca_data$filtered_D, file = file.path(outPath, paste0("PCA_data_after_imputation", suffix, ".csv")), row.names = FALSE)


  return(invisible(NULL))
}
