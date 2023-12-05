#' Calculate one MA plot for proteomics data
#'
#' @param sample_1              A numeric vector containing data of the first sample.
#' @param sample_2              A numeric vector containing data of the second sample.
#' @param do_log_transformation If \code{TRUE}, the data will be log-transformed.
#' @param alpha                 If \code{TRUE}, the data points will be transparent.
#' @param point_color           A character containing the colors of the data points.
#' @param ...                   Additional arguments for affy::ma.plot.
#'
#' @return An MA plot for two samples.
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' TODO!!!
#' 
#'}
#' 

MA_Plot_single <- function(sample_1, sample_2, 
                          do_log_transformation = FALSE, 
                          alpha = FALSE, 
                          point_color = "black", ...) {
  
  
  if(do_log_transformation) {
    sample_1 <- log2(sample_1)
    sample_2 <- log2(sample_2)
  }
  if(alpha) {
    point_color = alpha(point_color, 0.5)
  }
  
  M <- stats::na.omit(sample_1 - sample_2)
  A <- stats::na.omit((sample_1 + sample_2)/2)
  
  if (length(point_color) > 1) {
    na.ind <- attr(M, "na.action")
    point_color <- point_color[-na.ind]
  }
  
  
  affy::ma.plot(A = A, M = M, pch = 16, cex = 0.7, col = point_color, show.statistics = FALSE, ...)
}



#' Calculate MA plots for proteomics data set
#'
#' @param D           A data.frame of the data set.
#' @param do_log_transformation If \code{TRUE}, the data will be log-transformed.
#' @param output_path A character containing the path to a folder.
#' @param suffix      A character containing the suffix, with which the output file will be named.
#' @param labels      The sample labels for the title of the MA-Plot.
#' @param labels2     The second line in sample title (e.g. group membership).
#' @param maxPlots    A numeric containing the maximum number of MA plots that should be generated.
#' @param alpha       If \code{TRUE}, the data points will be transparent.
#' @param plot_height The height of the resulting MA plots.
#' @param plot_width  The width of the resulting MA plots.
#' @param ...         Additional arguments for affy::ma.plot.
#'
#' @return MA plots for all sample combinations.
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' TODO!!!
#' 
#'}
#' 

MA_Plots <- function(D,
                    do_log_transformation = FALSE,
                    output_path = "", suffix = "",
                    labels = 1:ncol(D), labels2 = colnames(D),
                    maxPlots = 5000,
                    alpha = FALSE,
                    plot_height = 15, plot_width = 15, ...) {
  
  
  number_states <- max(as.integer(as.factor(colnames(D))))
  number_plots <- choose(number_states,2)
  
  mess <- ""
  if (number_plots > maxPlots) {
    mess <- paste0(mess, "Number of MA-Plots (", number_plots, ") is higher than maxPlots (", maxPlots, ").\nPlease increase maxPlots to plot all MA-plots. \n")
  }
  
  num <- 0
  pb <- utils::txtProgressBar(min = 0,max = number_plots,char = "#",style = 3)
  
  grDevices::pdf(paste0(output_path, "MA_Plots_", suffix, ".pdf"), height = plot_height/2.54, width = plot_width/2.54)
  
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

      MA_Plot_single(D[,i], D[, j], log = do_log_transformation, main = main, ...)
    }
  }
  
  grDevices::dev.off()
  close(pb)
  mess <- paste0(mess, number_plots, " MA plots generated. \n")
  
  
  return("message" = mess)
}



