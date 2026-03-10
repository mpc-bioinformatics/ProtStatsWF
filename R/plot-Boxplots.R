#' Boxplots showing the distribution of intensities over all samples
#'
#' @param D_long                  \strong{data.frame} \cr
#'                                The data set given in long format.
#' @param method                  \strong{character} \cr
#'                                The method used. Options are "boxplot" and "violinplot".
#' @param groupColumn              \strong{logical} \cr
#' @param groupvar_name           \strong{character} \cr
#'                                The name for the group variable.
#' @param group_colours           \strong{character vector} \cr
#'                                The hex codes for the group colors.
#' @param base_size               \strong{numeric} \cr
#'                                The base size of the font.
#' @param lwd                     \strong{numeric} \cr
#'                                The line width of the boxplot.
#' @param outlier_size            \strong{numeric} \cr
#'                                The size of the outliers.
#'
#' @return boxplots and messages
#' @export
#'
#' @examples
#' \dontrun{
#' prepared_data <- prepareData(...)
#'
#'
#' boxplot <- Boxplots(D_long = prepared_data[["D_long"]])
#' }
#'
#'

Boxplots <- function(D_long,
                     method = "boxplot",
                     groupColumn = NULL,
                     group_colours = NULL,
                     base_size = 15,
                     lwd = 0.5,
                     outlier_size = 1) {

  # select only relevant columns
  D_long <- dplyr::select(D_long, c(".feature", ".sample", "intensity_norm", group = groupColumn))

  x_axis <- sort(unique(D_long$.sample)) # save the different states for later
  D_long <- D_long[!is.na(D_long$intensity_norm),] # remove NA values


  .sample <- intensity_norm <- group <- NULL
  if (!is.null(groupColumn)) {
   pl_boxplot <- ggplot2::ggplot(data = D_long, mapping = ggplot2::aes(x = .sample, y = intensity_norm, fill = group)) +
     ggplot2::labs(fill = groupColumn)
    if (!is.null(group_colours)) pl_boxplot <- pl_boxplot + ggplot2::scale_fill_manual(values = group_colours)
  } else {
    pl_boxplot <- ggplot2::ggplot(data = D_long, mapping = ggplot2::aes(x = .sample, y = intensity_norm))
  }


  pl_boxplot <- pl_boxplot +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggplot2::ylab("Log intensity") + ggplot2::xlab("Sample") +
    ggplot2::scale_x_discrete(limits = x_axis, drop = FALSE, na.translate = TRUE)


  if (method == "violinplot") {
    pl_boxplot <- pl_boxplot + ggplot2::geom_violin()
  }
  if (method == "boxplot") {
    pl_boxplot <- pl_boxplot + ggplot2::geom_boxplot(linewidth = lwd, outlier.size = outlier_size)
  }

  return(pl_boxplot)
}


















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



