#' Calculate one single MA plot for comparison of two samples
#'
#' @param sample1 **numeric** \cr
#' Vector of intensities for the first sample.
#' @param sample2 **numeric** \cr
#' Vector of intensities for the second sample..
#' @param alpha **numeric(1)** \cr
#' Transparency factor for data points. Default is 1, no transparency.
#' @param point_color **character** \cr
#' The color of the data points. Either a single value (colour for all points)
#' or vector containing a colour for each data point.
#' @param ... \cr
#' Additional arguments for [affy::ma.plot].
#'
#' @return Generates the MA plot for two samples.
#'
#' @seealso [MAPlots()] for calculation of MA-Plots for a whole dataset.
#'
#' @importFrom affy ma.plot
#' @importFrom scales alpha
#' @importFrom stats na.omit
#'
#'
#' @examples
#' file_proteins <- system.file("extdata", "proteins_HCC.csv",
#'   package = "ProtStatsWF")
#' file_clinical <- system.file("extdata", "clinical_data.csv",
#'   package = "ProtStatsWF")
#' D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
#'   proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
#'   sampleNameColumn = "Sample", verbose = FALSE)
#'
#' Intensities <- SummarizedExperiment::assay(D$SE)
#' s1 <- Intensities[,1]
#' s2 <- Intensities[,2]
#'
#' MA_Plot_single(sample_1 = s1, sample_2 = s2)
MAPlotSingle <- function(sample1, sample2,
                          alpha = 1, pointColour = "black", ...) {

  if (alpha) pointColour = scales::alpha(pointColour, alpha)

  M <- stats::na.omit(sample1 - sample2)
  A <- stats::na.omit((sample1 + sample2)/2)

  if (length(pointColour) > 1) {
    na.ind <- attr(M, "na.action")
    pointColour <- pointColour[-na.ind]
  }

  affy::ma.plot(A = A, M = M, pch = 16, cex = 0.7, col = pointColour,
                show.statistics = FALSE, ...)
}



#' Calculate MA plots for proteomics data set.
#'
#' @param D **matrix** \cr
#' Data set containing peptide or protein intensities in wide format.
#' @param outputPath **character(1)** \cr
#' Path to a folder for saving a pdf file with the MA-Plots. Folder must exist.
#' @param suffix **character(1)** \cr
#' Suffix for the output file. Should start with an underscore. Default is "".
#' @param labels **character** \cr
#' Vector containing sample labels for each sample. Default is
#' 'as.character(1:ncol(D))', i.e. the respective column number.
#' @param labels2 **character** \cr
#' Vector containing a second set of sample labels which are printed as a
#' second line. Default is 'colnames(D)', i.e. the respective column names.
#' @param maxPlots **integer(1)** \cr
#' The maximum number of MA plots that should be generated. Default is 5000.
#' This setting reduces the time and amount of space needed for MA-Plots if
#' the number of samples is high.
#' @param alpha **numeric(1)** \cr
#' Transparency factor for data points. Default is 1, no transparency.
#' @param plot_height **numeric(1)** \cr
#' The height of the resulting MA plots in cm. Default is 15.
#' @param plot_width **numeric(1)** \cr
#' The width of the resulting MA plots in cm. Default is 15.
#' @param ... \cr
#' Additional arguments for affy::ma.plot.
#' @param verbose **logical(1)** \cr
#' If TRUE (default), messages are printed out and a progress bar is shown.
#'
#' @return A pdf file containing the MA plots for all sample combinations.
#' @export
#'
#' @seealso [MAPlotSingle()] for internal function that generates a single plot.
#'
#' @importFrom checkmate assertCharacter assertDirectoryExists assertFlag
#' @importFrom checkmate assertMatrix assertNumeric
#' @importFrom grDevices dev.off pdf
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @examples
#' file_proteins <- system.file("extdata", "proteins_HCC.csv",
#'   package = "ProtStatsWF")
#' file_clinical <- system.file("extdata", "clinical_data.csv",
#'   package = "ProtStatsWF")
#' D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
#'   proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
#'   sampleNameColumn = "Sample", verbose = FALSE)
#'
#' Intensities <- SummarizedExperiment::assay(D$SE)
#'
#' MAPlots(Intensities)
MAPlots <- function(D,
                    outPath, suffix = "",
                    labels = as.character(1:ncol(D)), labels2 = colnames(D),
                    maxPlots = 5000,
                    alpha = 1,
                    plotHeight = 15, plotWidth = 15,
                    verbose = TRUE,
                    ...) {
  checkmate::assertMatrix(D, min.cols = 2, min.rows = 1)
  checkmate::assertDirectoryExists(outPath)
  checkmate::assertCharacter(suffix, len = 1)
  checkmate::assertCharacter(labels, len = ncol(D))
  checkmate::assertCharacter(labels2, len = ncol(D))
  checkmate::assertNumeric(maxPlots, lower = 0)
  checkmate::assertNumeric(alpha, lower = 0, upper = 1)
  checkmate::assertNumeric(plotHeight, lower = 0)
  checkmate::assertNumeric(plotWidth, lower = 0)
  checkmate::assertFlag(verbose)

  number_states <- max(as.integer(as.factor(colnames(D))))
  number_plots <- choose(number_states,2)

  if (number_plots > maxPlots & verbose) {
    message("Number of MA-Plots (", number_plots, ") is higher than maxPlots (",
            maxPlots, ").\nIncrease maxPlots to plot all MA-plots.")
  }

  if (verbose) {
    message("Generating MA plots...")
    pb <- utils::txtProgressBar(min = 0, max = min(maxPlots, number_plots),
                                char = "#", style = 3)
  }

  filename <- paste0("MA_Plots", suffix, ".pdf")
  grDevices::pdf(file.path(outPath, filename), height = plotHeight/2.54,
                 width = plotWidth/2.54)

  num <- 0
  br <- 0 # indicator for breaking outer loop
  for (i in 1:(ncol(D) - 1)) {
    if (br == 1) {
      break
    }

    for (j in (i + 1):ncol(D)) {
      # if maximum number of plots is reached, stop.
      if (num > maxPlots) {
        br <- 1
        break
      }

      if (is.null(labels2)) {
        main = paste(labels[i], labels[j])
      } else  {
        main = paste(labels[i], labels[j], "\n", labels2[i], labels2[j])
      }

      num <- num + 1
      if (verbose) utils::setTxtProgressBar(pb, num)

      MAPlotSingle(D[,i], D[, j], main = main, ...)
    }
  }
  grDevices::dev.off()
  if (verbose) {
    close(pb)
    message(number_plots, " MA plots generated.")
  }
  return(invisible(NULL))
}



