#' Calculate one MA plot for proteomics data.
#'
#' @param sample_1                \strong{numeric vector} \cr
#'                                The data of the first sample.
#' @param sample_2                \strong{numeric vector} \cr
#'                                The data of the second sample.
#' @param do_log_transformation   \strong{logical} \cr
#'                                If \code{TRUE}, the data will be log-transformed.
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
                          do_log_transformation = FALSE,
                          alpha = FALSE,
                          point_color = "black",
                          sampling = 1,
                          ...) {


  if(do_log_transformation) {
    sample_1 <- log2(sample_1)
    sample_2 <- log2(sample_2)
  }
  if(alpha) {
    point_color = alpha(point_color, 0.5)
  }

  M <- stats::na.omit(sample_1 - sample_2)
  A <- stats::na.omit((sample_1 + sample_2)/2)

  ## sample only parts of the data points for data sets with many data points
  if (sampling < 1) {
    ind_sample <- sample(1:length(M), size = ceiling(length(M) * sampling))
    M <- M[ind_sample]
    A <- A[ind_sample]
  }


  if (length(point_color) > 1) {
    na.ind <- attr(M, "na.action")
    point_color <- point_color[-na.ind]
  }


  affy::ma.plot(A = A, M = M, pch = 16, cex = 0.7, col = point_color, show.statistics = FALSE, ...)
}



#' Calculate MA plots for proteomics data set.
#'
#' @param D                       \strong{data.frame} \cr
#'                                The data set containing intensities of the sample.
#' @param do_log_transformation   \strong{logical} \cr
#'                                If \code{TRUE}, the data will be log-transformed.
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
                    do_log_transformation = FALSE,
                    output_path = "", suffix = "",
                    labels = 1:ncol(D), labels2 = colnames(D),
                    maxPlots = 5000,
                    alpha = FALSE,
                    plot_height = 15, plot_width = 15,
                    sampling = 1,
                    ...) {


  number_states <- max(as.integer(as.factor(colnames(D))))
  number_plots <- choose(number_states,2)

  mess <- ""
  if (number_plots > maxPlots) {
    mess <- paste0(mess, "Number of MA-Plots (", number_plots, ") is higher than maxPlots (", maxPlots, ").\nPlease increase maxPlots to plot all MA-plots. \n")
  }

  num <- 0
  pb <- utils::txtProgressBar(min = 0,max = number_plots,char = "#",style = 3)

  filename <- paste0("MA_Plots", suffix, ".pdf")
  grDevices::pdf(file.path(output_path, filename), height = plot_height/2.54, width = plot_width/2.54)

  for(i in 1:(ncol(D)-1)) {
    for (j in (i + 1):ncol(D)) {

      # if maximum number of plots is reached, stop.
      if (num > maxPlots) {
        grDevices::dev.off()
        mess <- paste0(mess, maxPlots, " MA plots generated. \n")
        return("message" = mess)
      }

      if (is.null(labels2)) {
        main = paste(labels[i], labels[j])
      } else  {
        main = paste(labels[i], labels[j], "\n", labels2[i], labels2[j])
      }

      num <- num + 1
      utils::setTxtProgressBar(pb, num)

      MA_Plot_single(D[,i], D[, j], do_log_transformation = do_log_transformation, main = main, sampling = sampling, ...)
    }
  }

  grDevices::dev.off()
  close(pb)
  mess <- paste0(mess, number_plots, " MA plots generated. \n")

  message(mess)

  return("message" = mess)
}



