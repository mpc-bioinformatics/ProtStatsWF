
#' Boxplots for biomarker candidates including the data point (jitter).
#'
#' @param D               \strong{data.frame} \cr
#'                        The data set containing only protein intensities, already filtered for interesting candidates.
#' @param protein.names   \strong{character vector} \cr
#'                        The protein names.
#' @param group           \strong{character factor} \cr
#'                        The group membership of the data.
#' @param group_colours   \strong{character vector} \cr
#'                        The hex codes or colours names for the group colors.
#' @param log_data        \strong{logical} \cr
#'                        If \code{TRUE}, the data will be log-transformed.
#' @param log_base        \strong{numeric} \cr
#'                        The base used, if \code{log_data = TRUE}.
#' @param plot_device     \strong{character} \cr
#'                        The plot device. Options are "pdf" or "png.
#' @param plot_height     \strong{numeric} \cr
#'                        The height if plot in cm.
#' @param plot_width      \strong{numeric} \cr
#'                        The width of plot in cm.
#' @param plot_dpi        \strong{numeric} \cr
#'                        The plot resolution (only for png).
#' @param groupvar_name   \strong{character} \cr
#'                        The name of the group variable (displayed in legend).
#' @param output_path     \strong{character} \cr
#'                        The path to a folder for the output.
#' @param suffix          \strong{character} \cr
#'                        The suffix for the output file.
#' @param plot_NA_level   \strong{logical} \cr
#'                        If \code{TRUE}, data points will be plotted if they are NA.
#'
#' @return Nothing, saves pdf or png files with boxplots to the output folder.
#' @export
#'
#' @examples 
#' 

Boxplots_candidates <- function(D,
                                protein.names,
                                group = NULL,
                                group_colours = NULL,
                                log_data = TRUE,
                                log_base = 2,
                                plot_device = "pdf",
                                plot_height = 15,
                                plot_width = 15,
                                plot_dpi = 200,
                                groupvar_name = "Group",
                                output_path,
                                suffix = NULL,
                                plot_NA_level = FALSE) {


  if (plot_device == "pdf") {
    grDevices::pdf(paste0(output_path, "/boxplots_candidates", suffix, ".pdf"),
        height = plot_height/2.54,
        width = plot_width/2.54)
  }

  pb <- pbapply::startpb(min = 0, max = nrow(D))
  for (i in c(1:nrow(D))) {

    # prepare data for single protein
    data <- data.frame(value = unlist(D[i,]), group = group)
    if (log_data) {
      data$value <- log(data$value, log_base)
    }

    if (!plot_NA_level) {
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
      ggplot2::labs(title = protein.names[i]) +
      ggplot2::geom_jitter() +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(fill = groupvar_name, y = "log2(intensity)", x = groupvar_name) +
      ggplot2::ggtitle(protein.names[i])


    if (is.null(group_colours)) {
      group_colours <- scales::hue_pal()(length(levels(group)))
    }
    names(group_colours) <- levels(data$group)
    plot <- plot + ggplot2::scale_fill_manual(values = group_colours, breaks = names(group_colours), drop = FALSE)


    if (plot_device == "png") {
      protein_names <- gsub("/", "_", protein.names[i])
      grDevices::png(paste0(output_path, "/boxplots_candidates", "_", protein_names, suffix, ".png"),
                     height = plot_height, width = plot_width, units = "cm", res = 300)
    }
    plot(plot)
    if (plot_device == "png") grDevices::dev.off()


    pbapply::setpb(pb, i)
  }
  invisible(NULL)

  if (plot_device == "pdf") grDevices::dev.off()

  return(invisible(NULL))
}


