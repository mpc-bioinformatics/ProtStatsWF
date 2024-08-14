
#' Boxplots for biomarker candidates including the data point (jitter)
#'
#' @param D data set containing only protein intensities, already filtered for interesting candidates
#' @param protein.names character vector containing protein names
#' @param group factor containing the groups
#' @param group_colours character vector of hex codes or colours names for the group colors
#' @param log_data TRUE, if the data should be log-transformed
#' @param log_base base for log-transformation
#' @param plot_device plot device, "pdf" or "png
#' @param plot_height height if plot in cm
#' @param plot_width width of plot in cm
#' @param plot_dpi plot resolution (only for png)
#' @param groupvar_name name of the group variable (displayed in legend)
#' @param output_path output path
#' @param suffix suffix for the output file
#' @param plot_NA_level TRUE, if data points should be plotted if they are NA
#'
#' @return nothing, saves pdf or png files with boxplots
#' @export
#'
#' @examples # TODO
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
      ggplot2::ggtitle(protein.names[i]) #+
      #ggplot2::scale_x_discrete(drop = FALSE)


    if (is.null(group_colours)) {
      group_colours <- scales::hue_pal()(length(levels(group)))
    }
    names(group_colours) <- levels(data$group)
    plot <- plot + ggplot2::scale_fill_manual(values = group_colours, breaks = names(group_colours), drop = FALSE)


    if (plot_device == "png") {
      grDevices::png(paste0(output_path, "/boxplots_candidates", "_", protein.names[i], suffix, ".png"),
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




### TODO: Verbindungslinien zwischen den Punkten, die zum gleichen Sample gehören?
