#' Simple volcano plot from p-value, fold change and significance category
#'
#' @param p                       \strong{numeric factor} \cr
#'                                The values for the p-values.
#' @param FC                      \strong{numeric factor} \cr
#'                                The values for the fold changes.
#' @param significance_category   \strong{character factor} \cr
#'                                The significance categories for the volcano plot.
#'
#' @param log_base_fc             \strong{numeric} \cr
#'                                The base for the fold changes log-transformation.
#' @param log_base_p              \strong{numeric} \cr
#'                                The base for the p-values log-transformation.
#' @param thres_p                 \strong{numeric} \cr
#'                                The threshold for the p-values.
#' @param thres_fc                \strong{numeric} \cr
#'                                The threshold for the fold changes.
#' @param colour1                 \strong{character} \cr
#'                                The color for not significant proteins.
#' @param colour2                 \strong{character} \cr
#'                                The color for significant proteins.
#' @param colour3                 \strong{character} \cr
#'                                The color for significant proteins after FDR correction.
#' @param symmetric_x             \strong{logical} \cr
#'                                If \code{TRUE}, x-axis limits will be made symmetric (not used if xlim is defined).
#' @param legend_position         \strong{character} \cr
#'                                The positioning of the legend.
#'                                Options are "none", "left", "right", "bottom", "top" and "inside".
#' @param base_size               \strong{numeric} \cr
#'                                The base size of the font.
#' @param xlim                    \strong{numeric} \cr
#'                                The limits for x-axis.
#' @param ylim                    \strong{numeric} \cr
#'                                The limits for y-axis.
#' @param alpha                   \strong{numeric} \cr
#'                                The transparency of the data points.
#' @param point_size              \strong{numeric} \cr
#'                                The size of the data points.
#'
#' @return A ggplot showing the volcano plot.
#' @export
#'
#' @seealso [VolcanoPlot_ttest()], [VolcanoPlot_ANOVA()]
#'
#' @examples
#'

VolcanoPlot <- function(p,
                        FC,
                        significance_category,
                        log_base_fc = 2,
                        log_base_p = 10,
                        thres_p = 0.05,
                        thres_fc = 2,
                        colour1 = "grey",
                        colour2 = "black",
                        colour3 = "orange",
                        symmetric_x = FALSE,
                        legend_position = "bottom",
                        base_size = NULL,
                        xlim = NULL,
                        ylim = NULL, 
                        alpha = 0.5,
                        point_size = 3) {


  ### transform p-values and fold changes and thresholds
  ### TODO: What if p-value is exactly zero?
  transformed_FC <- log(FC, base = log_base_fc) # default: log2(FC)
  transformed_p <- -log(p, base = log_base_p) # default: -log10(p)

  log_thres_p <- -log(thres_p, base = log_base_p)
  log_thres_fc <- log(thres_fc, base = log_base_fc)

  RES <<- data.frame(transformed_FC = transformed_FC,
                    transformed_p = transformed_p,
                    significance = significance_category)


  significance <- RES$significance

  plot <- ggplot2::ggplot(data = RES, ggplot2::aes(x = transformed_FC, y = transformed_p, colour = significance)) +
    ggplot2::geom_point(alpha = alpha, show.legend = TRUE, size = point_size) +
    ggplot2::scale_colour_manual(values = c("not significant" = colour1,
                                            "significant" = colour2,
                                            "significant after FDR correction" = colour3),
                                 drop = FALSE,
                                 na.translate = FALSE) +

    ggplot2::xlab(paste0("log",log_base_fc,"(FC)")) +
    ggplot2::ylab(paste0("-log",log_base_p,"(p)")) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = point_size*1.5)))


  if (!is.null(base_size)) {
    plot <- plot + ggplot2::theme_bw(base_size = base_size)
  } else {
    plot <- plot + ggplot2::theme_bw()
  }

  plot <- plot + ggplot2::theme(legend.position = legend_position)

  ### if xlim is not set, use max/min FC and make it symmetric
  if (symmetric_x) {
    xlim_tmp <- max(abs(RES$transformed_FC), na.rm = TRUE)
    xlim <- c(-xlim_tmp, xlim_tmp)
  }
  if (!is.null(xlim)) {
    plot <- plot + ggplot2::xlim(xlim)
  }
  if (!is.null(ylim)) {
    plot <- plot + ggplot2::ylim(ylim)
  }

  ## draw threshold lines
  plot <- plot + ggplot2::geom_hline(yintercept = log_thres_p, linetype = "dotted") +
    ggplot2::geom_vline(xintercept =  log_thres_fc, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = -log_thres_fc, linetype = "dotted")

  return(plot)

}



################################################################################
################################################################################
################################################################################



#' Volcano plot for a t-test result.
#'
#' @param RES               \strong{data.frame} \cr
#'                          The results from a t-test.
#' @param columnname_p      \strong{character} \cr
#'                          The column name for p-value.
#' @param columnname_padj   \strong{character} \cr
#'                          The columns name for adjusted p-value.
#' @param columnname_FC     \strong{character} \cr
#'                          The column name for fold change.
#'                          
#' @param thres_fc          \strong{numeric} \cr
#'                          The threshold for fold change.
#' @param thres_p           \strong{numeric} \cr
#'                          The threshold for p-value.
#' @param log_base_fc       \strong{numeric} \cr
#'                          The base for the fold changes log-transformation.
#' @param log_base_p        \strong{numeric} \cr
#'                          The base for the p-values log-transformation.
#' @param is_FC_log         \strong{logical} \cr
#'                          If \code{TRUE}, fold change is already log-transformed.
#' @param is_p_log          \strong{logical} \cr
#'                          If \code{TRUE}, p-value is already log-transformed.
#'                          
#' @param show_thres_line   \strong{logical} \cr
#'                          If \code{TRUE}, threshold lines will be shown.
#' @param groupname1        \strong{character} \cr
#'                          The name of first group.
#' @param groupname2        \strong{character} \cr
#'                          The name of second group.
#'                          
#' @param plot_height       \strong{numeric} \cr
#'                          The height of plot.
#' @param plot_width        \strong{numeric} \cr
#'                          The width of plot.
#' @param plot_dpi          \strong{integer} \cr
#'                          The resolution of plot.
#' @param plot_device       \strong{character} \cr
#'                          The plot device that is used for the resulting plot.
#'                          Options are "pdf" and "png".
#'                          
#' @param output_path       \strong{character} \cr
#'                          The path for output file.
#' @param suffix            \strong{character} \cr
#'                          The suffix for output file.
#' @param add_annotation    \strong{logical} \cr
#'                          If \code{TRUE}, annotation will be added.
#'                          
#' @param ...               Additional arguments for the plot.
#'                         
#'
#' @return A ggplot object of the volcano plot from a ttest result.
#' @export
#'
#' @seealso [VolcanoPlot()], [VolcanoPlot_ANOVA()], [add_labels()]
#'
#' @examples
#'

VolcanoPlot_ttest <- function(RES,
                        columnname_p = "p",
                        columnname_padj = "padj",
                        columnname_FC = "FC",
                        
                        thres_fc = 2,
                        thres_p = 0.05,
                        log_base_fc = 2,
                        log_base_p = 10,
                        is_FC_log = FALSE,
                        is_p_log = FALSE,
                        
                        show_thres_line = TRUE,
                        groupname1 = "group1",
                        groupname2 = "group2",
                        
                        plot_height = 15,
                        plot_width = 15,
                        plot_dpi = 300,
                        plot_device="pdf",
                        
                        output_path = NULL,
                        suffix = NULL,
                        add_annotation = TRUE,
                        ...) {


  # make check work
  #if(getRversion() >= "2.15.1")  utils::globalVariables(c("transformed_FC", "transformed_p", "significance"), add = FALSE)


  p <- RES[,columnname_p]
  padj <- RES[,columnname_padj]
  FC <- RES[,columnname_FC]

  if (is_FC_log) {
    FC <- log_base_fc^FC
  }
  if (is_p_log) {
    p <- log_base_p^(-p)
    padj <- log_base_p^(-padj)
  }



  RES$significance <- calculate_significance_categories_ttest(p = p,
                                                              p_adj = padj,
                                                              fc = FC,
                                                              thres_fc = thres_fc,
                                                              thres_p = thres_p)


  plot <- VolcanoPlot(p = p,
                      FC = FC,
                      significance_category = RES$significance,
                      ...)


  ## TODO: annotation of number of significant proteins
  # # here count the number of significant after FDR correction (>0 and <=0)
  # newX <- RES %>% filter(group %in% "significant after FDR correction")
  # lessthan0 <- 0
  # greaterthan0 <- 0
  #
  # newX <- t(newX$transformed_FC)
  #
  # for( x in newX){
  #   ifelse(x > 0, greaterthan0<-greaterthan0+1, lessthan0<-lessthan0+1)
  # }
  #
  # ### TODO: make annotation optional
  # ### TODO: add arrows (in current solution, arrows are too low)
  # if (add_annotation) {
  #   plot <- plot +
  #     annotate("text", x= xlim[1]+ (xlim[2]/10), y=ylim[2], label= paste0(groupname1, ": " , lessthan0))+#, " \U2191")) +
  #     annotate("text", x = xlim[2]-(xlim[2]/10),  y=ylim[2], label = paste0(groupname2 ,": " , greaterthan0))#+#, " \U2191"))
  # }


  #openxlsx::write.xlsx(x = RES, file=paste0(output_path,"transformedData", suffix, ".xlsx"), overwrite = TRUE,keepNA = TRUE)
  #ggplot2::ggsave(paste0(output_path,"Volcano_Plot", suffix, ".",plot_device),
  #       plot = plot, device = plot_device,
  #       height = plot_height, width = plot_width, dpi = plot_dpi, units = "cm")
  return(plot = plot) # its not necessary to return the data also, they are already present in the plot object
}



################################################################################
################################################################################
################################################################################



#' Volcano Plots for an ANOVA result.
#'
#' @param RES                    \strong{data.frame} \cr
#'                               The results from an ANOVA.
#' @param columnname_p_ANOVA     \strong{character} \cr
#'                               The column name for the p-values.
#' @param columnname_p_ANOVA_adj \strong{character} \cr
#'                               The column name for the adjusted p-values.
#' @param columns_FC             \strong{integer vector} \cr
#'                               The column indices of the fold changes.
#' @param columns_p_posthoc      \strong{integer vector} \cr
#'                               The column indices for the posthoc p-values.
#' @param log_base_fc            \strong{numeric} \cr
#'                               The base for the fold changes log-transformation.
#' @param log_base_p             \strong{numeric} \cr
#'                               The base for the p-values log-transformation.
#' @param thres_p                \strong{numeric} \cr
#'                               The threshold for p-value.
#' @param thres_fc               \strong{numeric} \cr
#'                               The threshold for fold change.
#' @param colour1                \strong{character} \cr
#'                               The color for not significant proteins.
#' @param colour2                \strong{character} \cr
#'                               The color for significant proteins.
#' @param colour3                \strong{character} \cr
#'                               The color for significant proteins after FDR correction.
#' @param symmetric_x            \strong{logical} \cr
#'                               If \code{TRUE}, x-axis limits will be made symmetric (not used if xlim is defined).
#' @param legend_position        \strong{character} \cr
#'                               The positioning of the legend.
#'                               Options are "none", "left", "right", "bottom", "top" and "inside".
#' @param base_size              \strong{numeric} \cr
#'                               The base size for theme.
#' @param xlim                   \strong{numeric} \cr
#'                               The limits for x-axis.
#' @param ylim                   \strong{numeric} \cr
#'                               The limits for y-axis.
#' @param add_labels             \strong{logical} \cr
#'                               If \code{TRUE}, labels will be added.
#'
#' @return A list of ggplots of the volcano plot from an ANOVA result.
#' @export
#'
#' @seealso [VolcanoPlot()], [VolcanoPlot_ttest()], [add_labels()]
#'
#' @examples
#'

VolcanoPlot_ANOVA <- function(RES,
                              columnname_p_ANOVA = "p.anova",
                              columnname_p_ANOVA_adj = "p.anova.fdr",
                              columns_FC, # spaltennnummer
                              columns_p_posthoc, # spaltennnummer
                              log_base_fc = 2,
                              log_base_p = 10,
                              thres_p = 0.05,
                              thres_fc = 2,
                              colour1 = "grey",
                              colour2 = "black",
                              colour3 = "orange",
                              symmetric_x = FALSE,
                              legend_position = "bottom",
                              base_size = NULL,
                              xlim = NULL,
                              ylim = NULL,
                              add_labels = FALSE) {

  nr_comparisons <- length(columns_p_posthoc)
  if (length(columns_FC) != nr_comparisons) stop("columns_FC and columns_p_posthoc must have the same length!")


  #nr_groups <- length(columns_p_posthoc)
  #if (length(columns_FC) != nr_groups) stop("columns_FC and columns_p_posthoc must have the same length!")

  p_anova <- RES[, columnname_p_ANOVA]
  p_anova_adj <- RES[, columnname_p_ANOVA_adj]
  #nr_comparisons <- choose(n = nr_groups, k = 2) # pairwise comparisons between two groups

  # names of the comparisons
  comp_names <- colnames(RES)[columns_FC]
  comp_names <- substring(comp_names, 4) # remove "FC_" at beginning
  comp_names <- stringr::str_replace_all(comp_names, "devided_by_", "vs")

  #comp_names <- stringr::str_replace_all(comp_names, "FC_", "")
  #comp_names <- stringr::str_replace_all(comp_names, "_", " ")

  Volcano_plots <- list()
  for (i in 1:nr_comparisons) {

    p_posthoc <- RES[, columns_p_posthoc[i]]
    fc <- RES[, columns_FC[i]]

    X <- cbind(p_anova, p_anova_adj, p_posthoc, fc)

    significance <- calculate_significance_categories_ANOVA(p_posthoc = p_posthoc, p_anova_adj = p_anova_adj, p_anova = p_anova, fc = fc, thres_fc=2, thres_p=0.05)

    X <- cbind(X, significance)


    plot <- VolcanoPlot(p = p_posthoc,
                        FC = fc,
                        significance_category = significance,
                        log_base_fc = log_base_fc,
                        log_base_p = log_base_p,
                        thres_p = thres_p,
                        thres_fc = thres_fc,
                        colour1 = colour1,
                        colour2 = colour2,
                        colour3 = colour3,
                        symmetric_x = symmetric_x,
                        legend_position = legend_position,
                        base_size = base_size,
                        xlim = xlim,
                        ylim = ylim)

    plot <- plot + ggplot2::ggtitle(comp_names[i])

    print(plot$data)


    if (add_labels) {
      plot <- add_labels(plot,
                 label_type = "FDR",
                 ind = NULL,
                 protein_name_column = NULL,
                 protein_names = RES$Gene.names) ### TODO: verallgemeinern
    }

    Volcano_plots[[i]] <- plot

  }
  return(Volcano_plots)
}




################################################################################
################################################################################
################################################################################


### which to label:
## significant after FDR correction
## significant
## top 10 on each side?




#' Add labels to a volcano plot.
#'
#' @param RES_Volcano           \strong{result from [VolcanoPlot()]} \cr
#'                              A list containing the data.frame with the transformed data and the ggplot object.
#' @param label_type            \strong{character} \cr
#'                              The label type.
#'                              Options are "FDR" (significant after correction) or "noFDR" (significant without correction) or "index" -> define indices to label
#' @param protein_name_column   \strong{character} \cr
#'                              The column name of the protein names in the RES data.frame.
#' @param ind                   \strong{integer} \cr
#'                              The index if the label_type is "index".
#' @param protein_names         \strong{character vector} \cr
#'                              The names of the proteins.
#'
#' @return A ggplot object with labels.
#' @export
#'
#' @seealso [VolcanoPlot()], [VolcanoPlot_ttest()], [VolcanoPlot_ANOVA()]
#'
#' @examples
#'

add_labels <- function(RES_Volcano,
                       label_type = "FDR",
                       ind = NULL,
                       protein_name_column = "Gene.names",
                       protein_names = NULL) {
  #### TODO: protein_name_column muss separat übergeben werden, weil dieSpalte nicht im Plot-Objekt vorhanden ist


  if (label_type == "FDR") {
    ind_label <- which(RES_Volcano$data$significance == "significant after FDR correction")
    if (length(ind_label) == 0) {
      warning("No significant proteins after FDR correction for labelling. Try changing label_type to noFDR or lower the fold change threshold.")
      return(RES_Volcano$plot)
    }
  }
  if (label_type == "noFDR") {
    ind_label <- which(RES_Volcano$data$significance %in% c("significant", "significant after FDR correction"))
    if (length(ind_label) == 0) {
      warning("No significant proteins for labelling. Try lowering the fold change threshold.")
      return(RES_Volcano$plot)
    }
  }
  if (label_type == "index") {
    ind_label <- ind
  }

  ind0 <<- ind_label

  print(RES_Volcano$data$transformed_FC[ind_label])

  ind_label_down <- ind_label[RES_Volcano$data$transformed_FC[ind_label] < 0]
  ind_label_up <- ind_label[RES_Volcano$data$transformed_FC[ind_label] > 0]

  ind1 <<- ind_label_up
  ind2 <<- ind_label_down

  labels_up <- rep(NA, nrow(RES_Volcano$data))
  labels_down <- rep(NA, nrow(RES_Volcano$data))
  if (is.null(protein_names)) {
    labels_up[ind_label_up] <- RES_Volcano$data[, protein_name_column][ind_label_up]
    labels_down[ind_label_down] <- RES_Volcano$data[, protein_name_column][ind_label_down]
  } else {
    labels_up[ind_label_up] <- protein_names[ind_label_up]
    labels_down[ind_label_down] <- protein_names[ind_label_down]
  }

  nudge_x = 0.2
  nudge_y = 0.2


  # change axis limits to have more space for labels
  xaxis_limits <- ggplot2::layer_scales(RES_Volcano)$x$get_limits()
  xaxis_limits <- xaxis_limits * 1.1

  yaxis_limits <- ggplot2::layer_scales(RES_Volcano)$y$get_limits()
  yaxis_limits[2] <- yaxis_limits[2] * 1.1 # only changge upper limit


  x1 <<- labels_up
  x2 <<- labels_down

  ### add labels
  plot <- RES_Volcano +
    ggplot2::ylim(yaxis_limits) + ggplot2::xlim(xaxis_limits) +
    # geom_point(data= RES_Volcano$RES[!is.na(labels_up) | !is.na(labels_down),],
    #             aes(x=transformed_FC,y=transformed_p)) +
    ggrepel::geom_label_repel(ggplot2::aes(label = labels_up), size = 2,
                              nudge_x = nudge_x, nudge_y = nudge_y, show.legend = FALSE,
                              max.overlaps = 20) +
    ggrepel::geom_label_repel(ggplot2::aes(label = labels_down), size = 2,
                              nudge_x = -nudge_x, nudge_y = nudge_y, show.legend = FALSE,
                              max.overlaps = 20)

  return(plot)
}



