#' t-test workflow
#'
#' @description
#' Workflow for t-test analysis of quantitative proteomics data.
#'
#' @details
#' This function performs a t-test to compare two experimental groups in a
#' quantitative proteomics dataset. The input \code{D} should be the list
#' returned by [prepareData()], which handles data import, log-transformation,
#' and normalisation. Because the assay data is already log-transformed by
#' [prepareData()], \code{logBeforeTest} defaults to \code{FALSE}.
#' The function generates a volcano plot, histograms of p-values and fold
#' changes as well as boxplots and a heatmap of the significant biomarker
#' candidates.
#'
#' @param D **list** \cr
#' Result from [prepareData()] containing the prepared data.
#' @param groupColumn **character(1)** \cr
#' Name of the column that contains the groups that are compared with a t-test.
#' Exactly two groups must be present.
#' @param sampleColumn **character(1)** \cr
#' Name of the column that identifies matched samples for a paired test.
#' Required only when \code{paired = TRUE}. Default is \code{NULL}.
#' @param proteinNameColumn **character(1)** \cr
#' Name of the column that contains protein identifiers (e.g. protein
#' accessions, gene names). Default is \code{"Protein"}. This column will be
#' used to label candidate proteins in boxplots and heatmap.
#' @param outputPath **character(1)** \cr
#' Path to the output folder. The folder must already exist.
#' @param suffix **character(1)** \cr
#' Suffix to add to the output file names. Default is "". Should ideally start
#' with an underscore "_".
#' @param paired **logical(1)** \cr
#' If \code{TRUE}, a paired test is performed. Default is \code{FALSE}
#' (unpaired t-test).
#' @param varEqual **logical(1)** \cr
#' If \code{TRUE}, variances in both groups are assumed equal. Default is
#' \code{FALSE} (assuming unequal variances).
#' @param logBeforeTest **logical(1)** \cr
#' If \code{TRUE}, data will be log-transformed before the test. Defaults to
#' \code{FALSE} because usually data prepared with [prepareData()] is already
#' log-transformed.
#' @param delogForFC **logical(1)** \cr
#' If \code{TRUE}, fold changes are computed on the original (de-log) scale.
#' Default is \code{TRUE}.
#' @param groupColours **character** \cr
#' Vector of colours for the two groups. Default is \code{NULL}
#'   (default ggplot2 colour palette).
#' @param significantAfterFDR **logical(1)** \cr
#' If \code{TRUE}, only proteins significant after FDR correction are shown in
#' boxplots and heatmap. Default is \code{TRUE}.
#' @param thresFC **numeric(1)** \cr
#' Fold change threshold used to classify a candidate as significant. Default is 2.
#' @param thresP **numeric(1)** \cr
#' P-value threshold used to classify a candidate as significant. Default is 0.05.
#' @param pValueZerosToMin **logical(1)** If \code{TRUE}, p-values equal to 0 are replaced
#'   by the next smallest observed p-value. Default is \code{TRUE}.
#' @param volcanoBaseSize **numeric(1)** \cr
#' Base size for the volcano plot. Default is 25.
#' @param histogramBaseSize **numeric(1)** \cr
#' Base size for the p-value and fold change histograms. Default is 15.
#' @param heatmapTextSize **numeric(1)** \cr
#' Text size for the heatmap. Default is 15.
#' @param maxValidValuesOff **integer(1)** Maximum number of valid values for a protein to be
#'   classified as "off". Default is \code{0}.
#' @param minValidValuesOn **integer(1)** Minimum number of valid values for a protein to be
#'   classified as "on". Default is \code{NULL} (set automatically to the smallest group size).
#' @param plotDevice **character(1)** \cr
#' Device to use for saving plots. Default is "pdf".
#' @param plotHeight_Volcano **numeric(1)** \cr
#' Height of the volcano plot in cm. Default is \code{15}.
#' @param plotWidth_Volcano **numeric(1)** \cr
#' Width of the volcano plot in cm. Default is \code{15}.
#' @param plotHeight_Histogram **numeric(1)** \cr
#' Height of the p-value and fold change histograms in cm. Default is \code{15}.
#' @param plotWidth_Histogram **numeric(1)** \cr
#' Width of the p-value and fold change histograms in cm. Default is \code{15}.
#' @param plotHeight_Boxplot **numeric(1)** \cr
#' Height of the boxplots of biomarker candidates in cm. Default is \code{15}.
#' @param plotWidth_Boxplot **numeric(1)** \cr
#' Width of the boxplots of biomarker candidates in cm. Default is \code{15}.
#' @param plotHeight_Heatmap **numeric(1)** \cr
#' Height of the heatmap in cm. Default is \code{15}.
#' @param plotWidth_Heatmap **numeric(1)** \cr
#' Width of the heatmap in cm. Default is \code{15}.
#' @param plotDPI **integer(1)** \cr
#' Plot resolution in DPI. Default is \code{300}.
#' @param assayName **character(1)** \cr
#' Name of the assay in \code{D$SE} to use as input data.
#' Default is \code{"intensity_norm"}, which corresponds to the output of
#' [prepareData()]. Only change if you do not directly use the output from
#' [prepareData()].
#' @param verbose **logical(1)** \cr
#' Whether to print messages and progress bars during the workflow. Default is TRUE.
#'
#' @return A list with one element \code{"message"}: a character string log of the workflow
#'   summarising settings and results. All output files are written to \code{outputPath}.
#' @export
#'
#' @importFrom checkmate assertCharacter assertClass assertDirectoryExists assertFlag
#' @importFrom checkmate assertInt assertList assertNumber assertSubset
#' @importFrom ComplexHeatmap draw
#' @importFrom ggplot2 ggsave
#' @importFrom grDevices dev.off pdf
#' @importFrom openxlsx write.xlsx
#' @importFrom pbapply pboptions
#' @importFrom scales hue_pal
#' @importFrom SummarizedExperiment assay colData rowData
#'
#' @seealso [workflow_ANOVA()] for more than two groups.\cr
#'          Functions used in this workflow:
#'          [prepareData()], [ttest()], [VolcanoPlot_ttest()], [pvalueFCHistogram()],
#'          [.calcSignCat_ttest()], [BoxplotsCandidates()],
#'          [heatmap()], [calculate_onoff()]
#'
#' @examples
#' \dontrun{
#' file_proteins <- system.file("extdata", "proteins_HCC.csv", package = "ProtStatsWF")
#' file_clinical <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")
#'
#' D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
#'                    proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
#'                    sampleNameColumn = "Sample", fileType = "csv")
#'
#' result <- workflow_ttest(D = D, groupColumn = "Group",
#'                          outputPath = tempdir())
#' }

workflow_ttest <- function(D,
                           groupColumn,
                           outputPath,

                           sampleColumn = NULL,
                           assayName = "intensity_norm",
                           groupColours = NULL,
                           proteinNameColumn = "Protein",

                           paired = FALSE,
                           varEqual = FALSE,
                           logBeforeTest = FALSE,
                           delogForFC = TRUE,
                           pValueZerosToMin = TRUE,

                           volcanoBaseSize = 25,
                           histogramBaseSize = 15,
                           heatmapTextSize = 15,

                           significantAfterFDR = TRUE,
                           thresFC = 2,
                           thresP = 0.05,
                           maxValidValuesOff = 0,
                           minValidValuesOn = NULL,

                           suffix = "",
                           plotDevice = "pdf",
                           plotHeight_Volcano = 15,
                           plotWidth_Volcano = 15,
                           plotHeight_Histogram = 15,
                           plotWidth_Histogram = 15,
                           plotHeight_Boxplot = 15,
                           plotWidth_Boxplot = 15,
                           plotHeight_Heatmap = 15,
                           plotWidth_Heatmap = 15,
                           plotDPI = 300,
                           verbose = TRUE
                           ) {

  if (!requireNamespace("amap", quietly = TRUE)) {
    stop("Package \"amap\" must be installed to cluster the heatmap.",
      call. = FALSE)
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package \"circlize\" must be installed for the heatmap legend.",
      call. = FALSE)
  }

  checkmate::assertList(D)
  checkmate::assertSubset("SE", names(D))
  checkmate::assertClass(D$SE, "SummarizedExperiment")
  checkmate::assertCharacter(assayName, len = 1)
  checkmate::assertSubset(assayName, names(SummarizedExperiment::assays(D$SE)))
  checkmate::assertCharacter(groupColumn, len = 1)
  checkmate::assertSubset(groupColumn, colnames(SummarizedExperiment::colData(D$SE)))
  checkmate::assertCharacter(sampleColumn, len = 1, null.ok = TRUE)
  if (!is.null(sampleColumn)) {
    checkmate::assertSubset(sampleColumn, colnames(SummarizedExperiment::colData(D$SE)))
  }
  checkmate::assertCharacter(proteinNameColumn, len = 1)
  checkmate::assertSubset(proteinNameColumn, colnames(SummarizedExperiment::rowData(D$SE)))
  checkmate::assertDirectoryExists(outputPath, access = "w")
  checkmate::assertCharacter(suffix, len = 1)
  checkmate::assertFlag(varEqual)
  checkmate::assertFlag(logBeforeTest)
  checkmate::assertFlag(delogForFC)
  checkmate::assertFlag(pValueZerosToMin)
  checkmate::assertFlag(significantAfterFDR)
  checkmate::assertInt(maxValidValuesOff, lower = 0)
  checkmate::assertInt(minValidValuesOn, lower = 0, null.ok = TRUE)
  checkmate::assertSubset(plotDevice, c("pdf", "jpeg", "tiff", "png", "svg"))
  checkmate::assertNumber(plotHeight_Volcano, lower = 0)
  checkmate::assertNumber(plotWidth_Volcano, lower = 0)
  checkmate::assertNumber(plotHeight_Histogram, lower = 0)
  checkmate::assertNumber(plotWidth_Histogram, lower = 0)
  checkmate::assertNumber(plotHeight_Boxplot, lower = 0)
  checkmate::assertNumber(plotWidth_Boxplot, lower = 0)
  checkmate::assertNumber(plotHeight_Heatmap, lower = 0)
  checkmate::assertNumber(plotWidth_Heatmap, lower = 0)
  checkmate::assertNumber(plotDPI, lower = 0)
  checkmate::assertFlag(verbose)

  #### Extract data from SummarizedExperiment ####

  DATA   <- as.data.frame(SummarizedExperiment::assay(D$SE, assayName))
  ID     <- as.data.frame(SummarizedExperiment::rowData(D$SE))
  group  <- droplevels(factor(SummarizedExperiment::colData(D$SE)[, groupColumn]))
  sample <- if (!is.null(sampleColumn)) factor(SummarizedExperiment::colData(D$SE)[, sampleColumn]) else NULL

  if (length(levels(group)) != 2) {
    stop("workflow_ttest requires exactly 2 groups in '", groupColumn, "', but found: ",
         paste(levels(group), collapse = ", "))
  }

  if (is.null(groupColours)) {
    groupColours <- scales::hue_pal()(length(levels(group)))
  } else {
    checkmate::assertCharacter(groupColours, len = length(levels(group)))
  }

  if (!verbose) {
    old_pbo <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(old_pbo), add = TRUE)
  }


  #### Calculate ttest ####

  test_results <- ttest(SE = D$SE,
                        assay = assayName,
                        groupColumn = groupColumn,
                        sampleColumn = sampleColumn,
                        paired = paired, varEqual = varEqual,
                        logBeforeTest = logBeforeTest, delogForFC = delogForFC, logBase = 2,
                        minObs = 3, minObsRatio = NULL, verbose = verbose)

  openxlsx::write.xlsx(test_results,
                       file = file.path(outputPath, paste0("results_ttest", suffix, ".xlsx")),
                       overwrite = TRUE, keepNA = TRUE)
  if (verbose) message(ifelse(paired, "Paired", "Unpaired"),
                       " t-test complete. Results saved.")

  if (pValueZerosToMin) {
    p_value_zero <- which(test_results$p == 0)
    if (length(p_value_zero) > 0) {
      next_smallest_value <- sort(unique(test_results$p))[2]
      test_results$p[test_results$p == 0] <- next_smallest_value
    }
  }


  fc_col_name <- paste0("FC_", levels(group)[[1]], "_divided_by_", levels(group)[[2]])


  #### Create Volcano Plot ####

  volcano_plot <- VolcanoPlot_ttest(RES = test_results,
                                    columnNameP = "p", columnNamePadj = "p.fdr",
                                    columnNameFC = fc_col_name, baseSize = volcanoBaseSize,
                                    thresFC = thresFC, thresP = thresP)

  ggplot2::ggsave(file.path(outputPath, paste0("volcano_plot", suffix, ".", plotDevice)),
                  plot = volcano_plot, device = plotDevice,
                  height = plotHeight_Volcano, width = plotWidth_Volcano, dpi = plotDPI, units = "cm")
  if (verbose) message("Volcano plot saved.")

  #### Create Histograms for p-values and fold changes ####

  histograms <- pvalueFCHistogram(RES = test_results,
                                            columnP = "p", columnPadj = "p.fdr",
                                            columnFC = fc_col_name, baseSize = histogramBaseSize)

  ggplot2::ggsave(file.path(outputPath, paste0("histogram_p_value", suffix, ".", plotDevice)),
                  plot = histograms[["histogram_p_value"]],
                  device = plotDevice, height = plotHeight_Histogram, width = plotWidth_Histogram,
                  dpi = plotDPI, units = "cm")
  ggplot2::ggsave(file.path(outputPath, paste0("histogram_adjusted_p_value", suffix, ".", plotDevice)),
                  plot = histograms[["histogram_adjusted_p_value"]],
                  device = plotDevice, height = plotHeight_Histogram, width = plotWidth_Histogram,
                  dpi = plotDPI, units = "cm")
  ggplot2::ggsave(file.path(outputPath, paste0("histogram_fold_change", suffix, ".", plotDevice)),
                  plot = histograms[["histogram_fold_change"]],
                  device = plotDevice, height = plotHeight_Histogram, width = plotWidth_Histogram,
                  dpi = plotDPI, units = "cm")
  if (verbose) message("p-value, adjusted p-value and fold change histograms saved.")

  #### Get significant candidates ####

  significance <- .calcSignCat_ttest(p = test_results$p,
                                     pAdj = test_results$p.fdr,
                                     fc = test_results[[fc_col_name]],
                                     thresFC = thresFC, thresP = thresP)

  candidates <- as.character(significance)

  if (significantAfterFDR) {
    candidates <- which(candidates == "significant after FDR correction")
  } else {
    candidates <- which(candidates == "significant" | candidates == "significant after FDR correction")
  }
  if (verbose) message("Found ", length(candidates), " significant candidate",
                       ifelse(length(candidates) == 1, "", "s"),
                       ifelse(significantAfterFDR, " after FDR correction.", "."))


  #### Create Boxplots of Biomarker Candidates ####

  if (length(candidates) > 0) {
    BoxplotsCandidates(SE = D$SE[candidates, ],
                        assay = assayName,
                        groupColumn = groupColumn,
                        proteinNameColumn = proteinNameColumn,
                        groupColours = groupColours,
                        suffix = suffix,
                        outputPath = outputPath,
                        plotDevice = plotDevice,
                        plotHeight = plotHeight_Boxplot,
                        plotWidth = plotWidth_Boxplot,
                        plotDPI = plotDPI,
                        verbose = verbose)
    if (verbose) message("Boxplots saved.")
  }

  #### Create Heatmap ####

  if (length(candidates) > 1) {
    # set.seed(14)
    t_heatmap <- heatmap(D = DATA[candidates, ],
                                     id = ID[candidates, ],
                                     groups = data.frame(Group = group),
                                     textSize = heatmapTextSize,
                                     verbose = verbose)

    if (!is.null(t_heatmap)) {
      grDevices::pdf(file.path(outputPath, paste0("heatmap", suffix, ".pdf")),
                     height = plotHeight_Heatmap/2.54, width = plotWidth_Heatmap/2.54)
      ComplexHeatmap::draw(t_heatmap[["heatmap"]])
      grDevices::dev.off()
      if (verbose) message("Heatmap saved.")
    }

  }


  #### Calculate on/off proteins ####

  if (is.null(minValidValuesOn)) {
    minValidValuesOn <- min(table(group))
  }

  on_off <- calculate_onoff(D = DATA,
                            id = ID,
                            group = group,
                            maxValidValuesOff = maxValidValuesOff,
                            minValidValuesOn = minValidValuesOn,
                            proteinNamesColumn = which(colnames(ID) == proteinNameColumn))

  openxlsx::write.xlsx(on_off, file = file.path(outputPath, paste0("results_onoff", suffix, ".xlsx")),
                       overwrite = TRUE, keepNA = TRUE)
  if (verbose) message("On/off analysis complete. Results saved.")


  return(list(test_results = test_results, significance = significance))
}



