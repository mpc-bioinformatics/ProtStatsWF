#' QC Workflow (QCQuant)
#'
#'
#' @description
#' Workflow for quality control of quantitative proteomics data
#'
#' @details
#' This function performs quality control of quantitative proteomics data.
#' The following plots are generated: a valid value plot, boxplots, MA-plots and
#' a PCA plot.
#'
#' @param D **list** \cr
#' Result from [prepareDataSE()] containing the imported data.
#' @param groupColumn **character(1)** \cr
#' Name of the column that contains the group information for colouring the
#' plots. Default is NULL (no group information for colouring).
#' @param group2Column **character(1)** \cr
#' Name of the column that contains a secondary group information (used only in
#' PCA_plot for shape of the data points). Default is NULL.
#' @param outPath **character(1)** \cr
#' Path to save the output files. The corresponding folder must exist.
#' @param outType **character(1)** \cr
#' Type of output tables. Default is "xlsx".
#' @param suffix **character(1)** \cr
#' Suffix to add to the output file names. Default is "". Should ideally start
#' with an underscore "_".
#' @param NAOut **character(1)** \cr
#' String to represent missing values in the output files. Default is "NA".
#' @param groupColours **character** \cr
#' Vector of colours for the groups. Default is NULL (default ggplot2 colour
#' palette will be used).
#' @param baseSize **numeric(1)** \cr
#' Base size for the plots. Default is 15.
#' @param plotDevice **character(1)** \cr
#' Device to use for saving plots. Default is "pdf".
#' @param plotHeight_BP **numeric(1)** \cr
#' Height of the boxplots in cm. Default is 10.
#' @param plotWidth_BP **numeric(1)** \cr
#' Width of the boxplots in cm. Default is 15.
#' @param plotHeight_VV **numeric(1)** \cr
#' Height of the valid value plot in cm. Default is 10.
#' @param plotWidth_VV **numeric(1)** \cr
#' Width of the valid value plot in cm. Default is 15.
#' @param plotHeight_PCA **numeric(1)** \cr
#' Height of the PCA plot in cm. Default is 15.
#' @param plotWidth_PCA **numeric(1)** \cr
#' Width of the PCA plot in cm. Default is 15.
#' @param plotHeight_MA **numeric(1)** \cr
#' Height of the MA plots in cm. Default is 15.
#' @param plotWidth_MA **numeric(1)** \cr
#' Width of the MA plots in cm. Default is 15.
#' @param plotDPI **numeric(1)** \cr
#' DPI for the plots. Default is 300.
#' @param boxplotMethod **character(1)** \cr
#' Method for boxplot. Options are "boxplot" (default) and "violinplot".
#' @param MAMaxPlots **numeric(1)** \cr
#' The maximum number of MA plots that should be generated. Default is 5000.
#' This setting reduces the time and amount of space needed for MA-Plots if
#' the number of samples is high.
#' @param MAAlpha **numeric(1)** \cr
#' Alpha value for MA plots. Default is 1 (no transparency).
#' @param PCAImputeMethod **character(1)** \cr
#' The missing value imputation method used for PCA. Options are "mean" or
#' "median" or "none" (default).
#' @param PCAPropNA **numeric(1)** \cr
#' The proportion of allowed missing values for PCA for a protein, before it is
#' discarded. It is automatically set to 0 if imputeMethod is "none".
#' @param PCAAlpha **numeric(1)** \cr
#' Alpha value (transparency) for points in the PCA. Default is 1.
#' @param PCALabel **logical(1)** \cr
#' Whether to label PCA points. Default is FALSE.
#' @param PCALabelSeed **numeric(1)** \cr
#' Seed for random number generator used for label placement in PCA. Default is
#' NA (no seed), which may lead to different label placement for multiple
#' executions of this function.
#' @param PCALabelSize **numeric(1)** \cr
#' Size of data point labels in the PCA. Default is 4.
#' @param PCAPointSize **numeric(1)** \cr
#'  Size of data points in the PCA. Default is 4.
#' @param verbose **logical(1)** \cr
#' Whether to print messages during the workflow. Default is TRUE.
#'
#' @return The workflow saves several plots and tables.
#' @export
#'
#' @importFrom checkmate assertCharacter assertDirectoryExists assertList
#' @importFrom checkmate assertNumber assertSubset
#' @importFrom ggplot2 ggsave
#' @importFrom scales hue_pal
#' @importFrom SummarizedExperiment assay colData
#' @importFrom utils write.csv
#'
#' @seealso Functions used in this workflow:
#'          [prepareDataSE()], [ValidValuePlot()], [Boxplots()], [MAPlots()],
#'          [PCA_Plot()].
#'
#' @examples
workflow_QC <- function(D,
                        groupColumn = NULL,
                        group2Column = NULL,

                        outPath,
                        outType = "xlsx",  ## TODO: does this work?
                        suffix = "",
                        NAOut = "NA", ## TODO: does this work?

                        groupColours = NULL,
                        baseSize = 15,
                        plotDevice = "pdf",
                        plotHeight_BP = 10,
                        plotWidth_BP = 15,
                        plotHeight_VV = 10,
                        plotWidth_VV = 15,
                        plotHeight_PCA = 15,
                        plotWidth_PCA = 15,
                        plotHeight_MA = 15,
                        plotWidth_MA = 15,
                        plotDPI = 300,

                        boxplotMethod = "boxplot",
                        MAMaxPlots = 5000,
                        MAAlpha = 1,
                        PCAImputeMethod = "none",
                        PCAPropNA = 0,
                        PCAAlpha = 1,
                        PCALabel = FALSE,
                        PCALabelSeed = NA,
                        PCALabelSize = 4,
                        PCAPointSize = 4,
                        verbose = TRUE
){
  checkmate::assertList(D)
  checkmate::assertSubset(c("SE", "D_long"), names(D))
  checkmate::assertDirectoryExists(outPath, access = "w")
  checkmate::assertSubset(outType, c("xlsx", "csv"))
  checkmate::assertCharacter(suffix, len = 1)
  checkmate::assertCharacter(NAOut, len = 1)
  checkmate::assertSubset(plotDevice, c("pdf", "jpeg", "tiff", "png", "svg"))
  checkmate::assertNumber(plotHeight_BP, lower = 0)
  checkmate::assertNumber(plotWidth_BP, lower = 0)
  checkmate::assertNumber(plotHeight_VV, lower = 0)
  checkmate::assertNumber(plotWidth_VV, lower = 0)
  checkmate::assertNumber(plotHeight_PCA, lower = 0)
  checkmate::assertNumber(plotWidth_PCA, lower = 0)
  checkmate::assertNumber(plotHeight_MA, lower = 0)
  checkmate::assertNumber(plotWidth_MA, lower = 0)
  checkmate::assertNumber(plotDPI, lower = 0)

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
                  plot = vv_plot$plot, device = plotDevice, height = plotHeight_VV,
                  width = plotWidth_VV, dpi = plotDPI, units = "cm")

  utils::write.csv(x = vv_plot$table, file = file.path(outPath,
                    paste0("D_validvalues", suffix, ".csv")), row.names = FALSE)



  #### Calculate Boxlots ####

  boxplots <- Boxplots(D_long = D$D_long,
                           groupColumn = groupColumn,
                           groupColours = groupColours,
                           baseSize = baseSize,
                           method = boxplotMethod, lwd = 0.5)

  ggplot2::ggsave(file.path(outPath, paste0("boxplot", suffix, ".", plotDevice)),
                  plot = boxplots, device = plotDevice, height = plotHeight_BP,
                  width = plotWidth_BP, dpi = plotDPI, units = "cm")





  #### Calculate MA Plot ####



  if (MAMaxPlots > 0) {
    ma_data <- MAPlots(D = SummarizedExperiment::assay(D$SE),  ### TODO: define which assay
                        outPath = outPath, suffix = suffix,
                        maxPlots = MAMaxPlots, alpha = MAAlpha,
                        plotHeight = plotHeight_MA,
                        plotWidth = plotWidth_MA,
                        verbose = verbose)
  }


  #### Calculate PCA Plot ####


  pca_data <- PCA_Plot(SE = D$SE,
                       groupForColour = groupColumn,
                       groupForShape = group2Column,
                       #impute = PCA_impute,
                       imputeMethod = PCAImputeMethod,
                       propNA = PCAPropNA,
                       scale = TRUE,
                       PCx = 1, PCy = 2,
                       #groupvar1_name = groupColumn,
                       #groupvar2_name = group2Column,
                       groupColours = groupColours,
                       alpha = PCAAlpha,
                       label = PCALabel,
                       labelSeed = PCALabelSeed,
                       labelSize = PCALabelSize,
                       xlim = NULL, ylim = NULL,
                       pointSize = PCAPointSize, baseSize = baseSize,
                       verbose = verbose)

  # Loadings <- as.data.frame(pca$rotation)
  # if (!is.null(id)) {
  #   Loadings <- cbind(id, Loadings)


  ggplot2::ggsave(file.path(outPath, paste0("PCA_plot", suffix, ".", plotDevice)), plot = pca_data[["plot"]],
                  device = plotDevice, height = plotHeight_PCA, width = plotWidth_PCA, dpi = plotDPI, units = "cm")
  utils::write.csv(x = pca_data$D_PCA_plot, file = file.path(outPath, paste0("D_PCA", suffix, ".csv")), row.names = FALSE)
  utils::write.csv(x = pca_data$filtered_D, file = file.path(outPath, paste0("PCA_data_after_imputation", suffix, ".csv")), row.names = FALSE)

  return(invisible(NULL))
}
