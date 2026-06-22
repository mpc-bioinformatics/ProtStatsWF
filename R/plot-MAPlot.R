#' Calculate one MA plot for proteomics data.
#'
#' @param sample_1                \strong{numeric vector} \cr
#'                                The data of the first sample.
#' @param sample_2                \strong{numeric vector} \cr
#'                                The data of the second sample.
#' @param alpha                   \strong{logical} \cr
#'                                If \code{TRUE}, the data points will be transparent.
#' @param point_color             \strong{character} \cr
#'                                The color of the data points.
#' @param sampling                \strong{numeric} \cr
#'                                The sampling rate. Useful to sample part of the data set for data sets on peptide/feature level with many data points.
#' @param ...                     Additional arguments for affy::ma.plot.
#'
#' @return Generates the MA plot for two samples.
#'
#' @seealso [MA_Plots()]
#'
#' @examples
#' \dontrun{
#' prepared_data <- prepareData(...)
#' data <- prepared_data[["D"]]
#'
#' s1 <- D[,1]
#' s1 <- D[,2]
#'
#' MA_Plot_single(sample_1 = s1, sample_2 = s2)
#'}
#'

MA_Plot_single <- function(sample_1, sample_2,
                          alpha = FALSE,
                          pointColour = "black",
                          sampling = 1,
                          ...) {

  if(alpha) {
    pointColour = alpha(pointColour, 0.5)
  }

  M <- stats::na.omit(sample_1 - sample_2)
  A <- stats::na.omit((sample_1 + sample_2)/2)

  ## sample only parts of the data points for data sets with many data points
  if (sampling < 1) {
    ind_sample <- sample(1:length(M), size = ceiling(length(M) * sampling))
    M <- M[ind_sample]
    A <- A[ind_sample]
  }


  if (length(pointColour) > 1) {
    na.ind <- attr(M, "na.action")
    pointColour <- pointColour[-na.ind]
  }


  affy::ma.plot(A = A, M = M, pch = 16, cex = 0.7, col = pointColour, show.statistics = FALSE, ...)
}



#' Calculate MA plots for proteomics data set.
#'
#' @param D                       \strong{data.frame} \cr
#'                                The data set containing intensities of the sample.
#' @param output_path             \strong{character} \cr
#'                                The path to a folder for the output.
#' @param suffix                  \strong{character} \cr
#'                                The suffix, with which the output file will be named. Should start with an underscore.
#' @param labels                  \strong{character} \cr
#'                                The sample labels for the title of the MA-Plot.
#' @param labels2                 \strong{character} \cr
#'                                The second line in sample title (e.g. group membership).
#' @param maxPlots                \strong{integer} \cr
#'                                The maximum number of MA plots that should be generated.
#' @param alpha                   \strong{logical} \cr
#'                                If \code{TRUE}, the data points will be transparent.
#' @param plot_height             \strong{numeric} \cr
#'                                The height of the resulting MA plots.
#' @param plot_width              \strong{numeric} \cr
#'                                The width of the resulting MA plots.
#' @param sampling                \strong{numeric} \cr
#'                                The sampling rate. Useful to sample part of the data set for data sets on peptide/feature level with many data points.
#' @param ...                     Additional arguments for affy::ma.plot.
#' @param verbose If TRUE, messages are printed out.
#'
#' @return A pdf file containing the MA plots for all sample combinations.
#' @export
#'
#' @seealso [MA_Plot_single()]
#'
#' @examples
#' \dontrun{
#' prepared_data <- prepareData(...)
#' out_path <- "/Users/thisuser/Documents/resultsFolder/"
#'
#' return_message <- MA_Plots(D = prepared_data[["D"]], output_path = out_path)
#'}
#'

MA_Plots <- function(D,
                    outPath = NULL, suffix = "",
                    labels = 1:ncol(D), labels2 = colnames(D),
                    maxPlots = 5000,
                    alpha = FALSE,
                    plotHeight = 15, plotWidth = 15,
                    sampling = 1, verbose = TRUE,
                    ...) {


  number_states <- max(as.integer(as.factor(colnames(D))))
  number_plots <- choose(number_states,2)

  if (number_plots > maxPlots) {
    message("Number of MA-Plots (", number_plots, ") is higher than maxPlots (", maxPlots, ").\nPlease increase maxPlots to plot all MA-plots.")
  }

  ## TODO: disable progress bar is verbose = FALSE
  num <- 0
  pb <- utils::txtProgressBar(min = 0,max = number_plots,char = "#",style = 3)

  if (!is.null(outPath)) {
    filename <- paste0("MA_Plots", suffix, ".pdf")
    grDevices::pdf(file.path(outPath, filename), height = plotHeight/2.54, width = plotWidth/2.54)
  }

  for(i in 1:(ncol(D)-1)) {
    for (j in (i + 1):ncol(D)) {

      # if maximum number of plots is reached, stop.
      if (num > maxPlots) {
        grDevices::dev.off()
        message(maxPlots, " MA plots generated.")
      }

      if (is.null(labels2)) {
        main = paste(labels[i], labels[j])
      } else  {
        main = paste(labels[i], labels[j], "\n", labels2[i], labels2[j])
      }

      num <- num + 1
      utils::setTxtProgressBar(pb, num)

      MA_Plot_single(D[,i], D[, j], main = main, sampling = sampling, ...)
    }
  }

  if (!is.null(outPath)) {
    grDevices::dev.off()
  }
  close(pb)

  if (verbose) message(number_plots, " MA plots generated.")

  return(invisible(NULL))
}



