#' t-test workflow
#'
#' @description
#' Workflow for t-test analysis of quantitative proteomics data.
#'
#' @details
#' This function performs a t-test to compare two experimental groups in a
#' quantitative proteomics dataset. The input \code{D} should be the list
#' returned by [prepareDataSE()], which handles data import, log-transformation,
#' and normalisation. Because the assay data is already log-transformed by
#' [prepareDataSE()], \code{logBeforeTest} defaults to \code{FALSE}.
#' The function generates a volcano plot, histograms of p-values and fold
#' changes as well as boxplots and a heatmap of the significant biomarker
#' candidates.
#'
#' @param D **list** \cr
#' Result from [prepareDataSE()] containing the prepared data.
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
#' \code{FALSE} because usually data prepared with [prepareDataSE()] is already
#' log-transformed.
#' @param delogForFC **logical(1)** \cr
#' If \code{TRUE}, fold changes are computed on the original (de-log) scale.
#' Default is \code{TRUE}.
#' @param groupColours **character** \cr
#' Vector of colours for the two groups. Default is \code{NULL}
#'   (default ggplot2 colour palette).

### TODO: hier fehlt p-value and fold change cutoff!
#' @param significantAfterFDR **logical(1)** \cr
#' If \code{TRUE}, only proteins significant after FDR correction are shown in
#' boxplots and heatmap. Default is \code{TRUE}.
#' @param pValueZerosToMin **logical(1)** If \code{TRUE}, p-values equal to 0 are replaced
#'   by the next smallest observed p-value. Default is \code{TRUE}.
#' @param baseSize **numeric(1)** \cr
#' Base size for the plots. Default is 15.
# TODO: use in all plots, not only volcano!
#'
#'
#' @param maxValidValuesOff **integer(1)** Maximum number of valid values for a protein to be
#'   classified as "off". Default is \code{0}.
#' @param minValidValuesOn **integer(1)** Minimum number of valid values for a protein to be
#'   classified as "on". Default is \code{NULL} (set automatically to the smallest group size).

#' @param plotDevice **character(1)** \cr
#' Device to use for saving plots. Default is "pdf".
#' @param plotHeight **numeric(1)** \cr
#' Plot height in cm. Default is \code{15}.
#' @param plotWidth **numeric(1)** \cr
#' Plot width in cm. Default is \code{15}.
# TODO: separate height/width for different plots.
#' @param plotDPI **integer(1)** \cr
#' Plot resolution in DPI. Default is \code{300}.
#' @param assayName **character(1)** \cr
#' Name of the assay in \code{D$SE} to use as input data.
#' Default is \code{"intensity_norm"}, which corresponds to the output of
#' [prepareDataSE()]. Only change if you do not directly use the output from
#' [prepareDataSE()].
#' @param verbose **logical(1)** \cr
#' Whether to print messages and progress bars during the workflow. Default is TRUE.
#'
#' @return A list with one element \code{"message"}: a character string log of the workflow
#'   summarising settings and results. All output files are written to \code{outputPath}.
#' @export
#'
#' @seealso [workflow_ANOVA()] for more than two groups.\cr
#'          Functions used in this workflow:
#'          [prepareDataSE()], [ttest()], [VolcanoPlot_ttest()], [pvalue_foldchange_histogram()],
#'          [calculate_significance_categories_ttest()], [Boxplots_candidates()],
#'          [Heatmap_with_groups()], [calculate_onoff()]
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

                           significantAfterFDR = TRUE,
                           thresFC = 2,
                           thresP = 0.05,
                           maxValidValuesOff = 0,
                           minValidValuesOn = NULL,

                           suffix = "",
                           plotDevice = "pdf",
                           plotHeight = 15,
                           plotWidth = 15,
                           plotDPI = 300,
                           verbose = TRUE
                           ) {

  #### Extract data from SummarizedExperiment ####

  DATA   <- as.data.frame(SummarizedExperiment::assay(D$SE, assayName))
  ID     <- as.data.frame(SummarizedExperiment::rowData(D$SE))
  group  <- droplevels(factor(SummarizedExperiment::colData(D$SE)[, groupColumn]))
  sample <- if (!is.null(sampleColumn)) factor(SummarizedExperiment::colData(D$SE)[, sampleColumn]) else NULL

  if (length(levels(group)) != 2) {
    stop("workflow_ttest requires exactly 2 groups in '", groupColumn, "', but found: ",
         paste(levels(group), collapse = ", "))
  }

  if (is.null(groupColours)) groupColours <- scales::hue_pal()(length(levels(group)))

  if (!verbose) {
    old_pbo <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(old_pbo), add = TRUE)
  }


  #### Calculate ttest ####

  test_results <- ttest(D = DATA, id = ID,
                        group = group, sample = sample,
                        paired = paired, varEqual = varEqual,
                        logBeforeTest = logBeforeTest, delogForFC = delogForFC, logBase = 2,
                        minObsPerGroup = 3, minObsPerGroupRatio = NULL)

  openxlsx::write.xlsx(test_results,
                       file = file.path(outputPath, paste0("results_ttest", suffix, ".xlsx")),
                       overwrite = TRUE, keepNA = TRUE)
  if (verbose) message(ifelse(paired, "Paired", "Unpaired"), " t-test complete. Results saved.")

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
                                    columnNameFC = fc_col_name, base_size = volcanoBaseSize,
                                    thresFC = thresFC, thresP = thresP)

  ggplot2::ggsave(file.path(outputPath, paste0("volcano_plot", suffix, ".", plotDevice)),
                  plot = volcano_plot, device = plotDevice,
                  height = plotHeight, width = plotWidth, dpi = plotDPI)
  if (verbose) message("Volcano plot saved.")

  #### Create Histograms for p-values and fold changes ####

  histograms <- pvalue_foldchange_histogram(RES = test_results,
                                            columnNameP = "p", columnNamePadj = "p.fdr",
                                            columnNameFC = fc_col_name)

  ggplot2::ggsave(file.path(outputPath, paste0("histogram_p_value", suffix, ".", plotDevice)),
                  plot = histograms[["histogram_p_value"]],
                  device = plotDevice, height = plotHeight, width = plotWidth, dpi = plotDPI)
  ggplot2::ggsave(file.path(outputPath, paste0("histogram_adjusted_p_value", suffix, ".", plotDevice)),
                  plot = histograms[["histogram_adjusted_p_value"]],
                  device = plotDevice, height = plotHeight, width = plotWidth, dpi = plotDPI)
  ggplot2::ggsave(file.path(outputPath, paste0("histogram_fold_change", suffix, ".", plotDevice)),
                  plot = histograms[["histogram_fold_change"]],
                  device = plotDevice, height = plotHeight, width = plotWidth, dpi = plotDPI)
  if (verbose) message("p-value, adjusted p-value and fold change histograms saved.")

  #### Get significant candidates ####

  significance <- calculate_significance_categories_ttest(p = test_results[["p"]],
                                                          pAdj = test_results[["p.fdr"]],
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
    Boxplots_candidates(D = DATA[candidates, ],
                        proteinNames = ID[candidates, proteinNameColumn],
                        group = group,
                        groupColours = groupColours,
                        suffix = suffix,
                        outputPath = outputPath,
                        logData = logBeforeTest)
    if (verbose) message("Boxplots saved.")
  }

  #### Create Heatmap ####

  if (length(candidates) > 1) {
    # set.seed(14)
    t_heatmap <- Heatmap_with_groups(D = DATA[candidates, ],
                                     id = ID[candidates, ],
                                     groups = group,
                                     verbose = verbose)

    if (!is.null(t_heatmap)) {
      grDevices::pdf(file.path(outputPath, paste0("heatmap", suffix, ".pdf")),
                     height = plotHeight, width = plotWidth)
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



