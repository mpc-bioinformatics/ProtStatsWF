#' QC Workflow (QCQuant)
#'
#'
#' @description
#' Workflow for quality control of quantitative proteomics data
#'
#' @details
#' This function performs quality control of quantitative proteomics data.
#' The following plots are generated: a valid value plot, boxplots, MA-plots and a PCA plot.
#'
#' @param D **list** Result from [prepareDatSE()] containing the imported data.
#' @param groupColumn **character(1)** Column name in the colData(D$SE) that contains the group
#' information for colouring the plots. Default is NULL (no group information for colouring).
#' @param group2Column **character(1)** Column name in the colData(D$SE) that
#' contains a secondary group information (used in PCA_plot as shape). Default is NULL.
#' @param groupColours **character** Vector of colours for the groups. Default is NULL
#' (default ggplot2 colour palette will be used.)
#' @param outPath **character(1)** Path to save the output files. The corresponding folder must exist.
#' @param outType **character(1)** Type of output file. Default is xlsx.
#' @param suffix **character(1)** Suffix to add to the output file names. Default is "".
#' @param NAOut **character(1)** String to represent missing values in the output files. Default is "NA".
#' @param verbose **logical(1)** Whether to print messages during the workflow. Default is TRUE.
#' @param baseSize **numeric(1)** Base size for the plots. Default is 15.
#' @param plotDevice **character(1)** Device to use for saving plots. Default is "pdf".
#' @param plotHeight_BP_VV **numeric(1)** Height of the boxplot and valid value plots in cm. Default is 10.
#' @param plotWidth_BP_VV **numeric(1)** Width of the boxplot and valid value plots in cm. Default is 15.
#' @param plotHeight_PCA_MA **numeric(1)** Height of the PCA and MA plots in cm. Default is 15.
#' @param plotWidth_PCA_MA **numeric(1)** Width of the PCA and MA plots in cm. Default is 15.
#' @param plotDPI **numeric(1)** DPI for the plots. Default is 300.
#' @param boxplotMethod **character(1)** Method for boxplot. Default is "boxplot".
#' @param MAMaxPlots **numeric(1)** Maximum number of MA plots to be drawn. Default is 5000.
#' @param MAAlpha **numeric(1)** Alpha value for MA plots. Default is 1.
#' @param PCAImputeMethod **character(1)** Method for missing value imputation for the PCA. Default is "mean".
#' @param PCAPropNA **numeric(1)** Allowed proportion of missing values for PCA. Default is 0.
#' Proteins with more missing values will be removed from the PCA. The others will be imputed.
#' @param PCAScale **logical(1)** Whether to scale data for the PCA. Default is TRUE.
  #' @param PCAAlpha **numeric(1)** Alpha value (transparency) for points in the PCA. Default is 1.
#' @param PCALabel **logical(1)** Whether to label PCA points. Default is FALSE.
#' @param PCALabelSeed **numeric(1)** Seed for data point labels in the PCA. Default is NA.
#' @param PCALabelSize **numeric(1)** Size of data point labels in the PCA. Default is 4.
#' @param PCAXlim **numeric(2)** Limits for the x-axis in the PCA plot. Default is NULL.
#' @param PCAYlim **numeric(2)** Limits for the y-axis in the PCA plot. Default is NULL.
#' @param PCAPointSize **numeric(1)** Size of data points in the PCA. Default is 4.
#'
#' @return The workflow saves several plots and files and returns a message log of the workflow.
#' @export
#'
#' @seealso Functions used in this workflow:
#'          [prepareDataSE()], [ValidValuePlot()], [Boxplots()], [MA_Plots()], [PCA_Plot()].
#'
#' @examples

workflow_QC <- function(D,
                        groupColumn = NULL,
                        group2Column = NULL,
                        groupColours = NULL,

                        outPath,
                        outType = "xlsx",
                        suffix = "",  ## TODO: if it doesnt start with underscore, add it automatically"
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

  ### TODO: na out and other output types

  # prepare group colours
  ### TODO: group colours should be a named vector if possible
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
    ma_data <- MA_Plots(D = SummarizedExperiment::assay(D$SE),  # TODO: give SE directly, not assay
                        outPath = outPath, suffix = suffix,
                        #labels = 1:ncol(SummarizedExperiment::assay(D$SE)),
                        #labels2 = colnames(SummarizedExperiment::assay(D$SE)),
                        maxPlots = MAMaxPlots, alpha = MAAlpha,
                        plotHeight = plotHeight_PCA_MA,
                        plotWidth = plotWidth_PCA_MA,
                        sampling = 1, verbose = verbose)
  }


  #### Calculate PCA Plot ####


  pca_data <- PCA_Plot(SE = D$SE,
                       groupForColour = groupColumn,
                       groupForShape = group2Column,
                       #impute = PCA_impute,
                       imputeMethod = PCAImputeMethod,
                       propNA = PCAPropNA,
                       scale. = PCAScale,
                       PCx = 1, PCy = 2,
                       #groupvar1_name = groupColumn,
                       #groupvar2_name = group2Column,
                       groupColours = groupColours,
                       alpha = PCAAlpha,
                       label = PCALabel,
                       labelSeed = PCALabelSeed,
                       labelSize = PCALabelSize,
                       xlim = PCAXlim, ylim = PCAYlim,
                       pointSize = PCAPointSize, baseSize = baseSize)

### TODO:: extract and save loadings
  # Loadings <- as.data.frame(pca$rotation)
  # if (!is.null(id)) {
  #   Loadings <- cbind(id, Loadings)


  ggplot2::ggsave(file.path(outPath, paste0("PCA_plot", suffix, ".", plotDevice)), plot = pca_data[["plot"]],
                  device = plotDevice, height = plotHeight_PCA_MA, width = plotWidth_PCA_MA, dpi = plotDPI, units = "cm")
  utils::write.csv(x = pca_data$D_PCA_plot, file = file.path(outPath, paste0("D_PCA", suffix, ".csv")), row.names = FALSE)
  utils::write.csv(x = pca_data$filtered_D, file = file.path(outPath, paste0("PCA_data_after_imputation", suffix, ".csv")), row.names = FALSE)

  return(invisible(NULL))
}
