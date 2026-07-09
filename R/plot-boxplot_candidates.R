
#' Boxplots for biomarker candidates including the data point (jitter).
#'
#' @param D               \strong{data.frame} \cr
#'                        The data set containing only protein intensities,
#'                        already filtered for interesting candidates. Usually
#'                        the intensities should already be log-transformed.
#' @param proteinNames    \strong{character vector} \cr
#'                        The protein names.
#' @param group           \strong{character factor} \cr
#'                        The group membership of the data.
#' @param groupColours    \strong{character vector} \cr
#'                        The hex codes or colours names for the group colors.
#' @param plotDevice      \strong{character} \cr
#'                        The plot device. Options are "pdf" or "png.
#' @param plotHeight      \strong{numeric} \cr
#'                        The height if plot in cm.
#' @param plotWidth       \strong{numeric} \cr
#'                        The width of plot in cm.
#' @param plotDPI         \strong{numeric} \cr
#'                        The plot resolution (only for png).
#' @param groupVarName    \strong{character} \cr
#'                        The name of the group variable (displayed in legend).
#' @param outputPath      \strong{character} \cr
#'                        The path to a folder for the output.
#' @param suffix          \strong{character} \cr
#'                        The suffix for the output file.
#' @param plotNALevel     \strong{logical} \cr
#'                        If \code{TRUE}, data points will be plotted if the group variable is NA.
#'
#' @return Nothing, saves pdf or png files with boxplots to the output folder.
#' @export
#'
#' @examples
#'

Boxplots_candidates <- function(D,
                                proteinNames,
                                group = NULL,
                                groupColours = NULL,
                                logData = TRUE,
                                logBase = 2,
                                plotDevice = "pdf",
                                plotHeight = 15,
                                plotWidth = 15,
                                plotDPI = 200,
                                groupVarName = "Group",
                                outputPath,
                                suffix = NULL,
                                plotNALevel = FALSE) {


  if (plotDevice == "pdf") {
    grDevices::pdf(paste0(outputPath, "/boxplots_candidates", suffix, ".pdf"),
                   height = plotHeight/2.54,
                   width = plotWidth/2.54)
  }

  pb <- pbapply::startpb(min = 0, max = nrow(D))
  for (i in c(1:nrow(D))) {

    # prepare data for single protein
    data <- data.frame(value = unlist(D[i,]), group = group)
    # if (logData) {
    #   data$value <- log(data$value, logBase)
    # }

    if (!plotNALevel) {
      data <- stats::na.omit(data)
      if (nrow(data) == 0) {
        warning(paste0("Only missing values for row ", i, " in data set."))
        next
      }
    }
    data$group <- factor(data$group, levels = levels(group))

    value <- NULL # to silence notes while checking package

    plot <- ggplot2::ggplot(data = data, ggplot2::aes(x = group, y = value, fill = group)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::labs(title = proteinNames[i]) +
      ggplot2::geom_jitter() +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(fill = groupVarName, y = "log2(intensity)", x = groupVarName) +
      ggplot2::ggtitle(proteinNames[i])

    if (is.null(groupColours)) {
      groupColours <- scales::hue_pal()(length(levels(group)))
    }
    names(groupColours) <- levels(data$group)
    plot <- plot + ggplot2::scale_fill_manual(values = groupColours, breaks = names(groupColours), drop = FALSE)


    if (plotDevice == "png") {
      protein_names <- gsub("/", "_", proteinNames[i])
      grDevices::png(paste0(outputPath, "/boxplots_candidates", "_", protein_names, suffix, ".png"),
                     height = plotHeight, width = plotWidth, units = "cm", res = 300)
    }
    plot(plot)
    if (plotDevice == "png") grDevices::dev.off()


    pbapply::setpb(pb, i)
  }

  if (plotDevice == "pdf") grDevices::dev.off()

  return(invisible(NULL))
}



