
#' Calculation of significance categories for a volcano plot
#'
#' This function groups all proteins into the following three catrgories based
#' on the p-values (with and without FDR correction) and fold changes:
#' 1) not signifcant (p > thres_p or fc < thres_fc)
#' 2) significant (p < thres_p and fc > thres_fc, but p_adj > thres_p)
#' 3) significant after FDR-correction (p_adj < thres_p and fc > thres_fc)
#'
#' @param p vector of p-values before FDR-correction
#' @param p_adj vector of p-values after FDR-correction
#' @param fc vector of fold changes
#' @param thres_fc threshold for the fold changes
#' @param thres_p threshold for the p-values
#'
#' @return a factor with the three significance categories
#' @export
#'
#' @examples # TODO
calculate_significance_categories <- function(p, p_adj, fc, thres_fc=2, thres_p=0.05) {

  significance <- dplyr::case_when(
    p_adj <= thres_p & (fc >= thres_fc | fc <= 1/thres_fc) & !is.na(p) ~ "significant after FDR correction",
    p_adj > thres_p & p <= thres_p & (fc >= thres_fc | fc <= 1/thres_fc) & !is.na(p) ~ "significant",
    (p > thres_p | (fc < thres_fc & fc > 1/thres_fc)) & !is.na(p) ~ "not significant",
    is.na(p) ~ NA_character_
  )

  significance <- factor(significance, levels = c("not significant", "significant", "significant after FDR correction"))

  return(significance)
}








#' Volcano plot
#'
#' @param RES result table from t-test
#' @param columnname_p column name for p-value
#' @param columnname_padj columns name for adjusted p-value
#' @param columnname_FC column name for fold change
#' @param log_base_fc base for log transformation of fold change
#' @param log_base_p base for log transformation of p-value
#' @param is_FC_log TRUE if fold change is already log-transformed
#' @param is_p_log FALSE if p-value is already log-transformed
#' @param thres_fc threshold for fold change
#' @param thres_p threshold for p-value
#' @param show_thres_line TRUE if threshold lines should be shown
#' @param colour1 colour for not significant proteins
#' @param colour2 colour for significant proteins
#' @param colour3 colour for significant proteins after FDR correction
#' @param groupname1 name of first group
#' @param groupname2 name of second group
#' @param xlim limits for x-axis
#' @param ylim limits for y-axis
#' @param symmetric_x if TRUE, x-axis limits will be made symmetric (not used if xlim is defined)
#' @param legend_position position of the legend (e.g., "bottom" or "right")
#' @param plot_height height of plot
#' @param plot_width width of plot
#' @param plot_dpi resolution of plot
#' @param plot_device plot device
#' @param output_path path for output
#' @param suffix suffix for output file
#' @param base_size base size for theme
#' @param add_annotation TRUE if annotation should be added
#'
#' @return ggplot object
#' @export
#'
#' @examples # TODO
VolcanoPlot <- function(RES,
                        columnname_p = "p",
                        columnname_padj = "padj",
                        columnname_FC = "FC",
                        log_base_fc = 2,
                        log_base_p = 10,
                        is_FC_log = FALSE,
                        is_p_log = FALSE,
                        thres_fc = 2,
                        thres_p = 0.05,
                        show_thres_line = TRUE,
                        colour1 = "grey",
                        colour2 = "black",
                        colour3 = "orange",
                        groupname1 = "group1",
                        groupname2 = "group2",
                        xlim = NULL,
                        ylim = NULL,
                        symmetric_x = FALSE,
                        legend_position = "bottom",
                        plot_height = 15,
                        plot_width = 15,
                        plot_dpi = 300,
                        plot_device="pdf",
                        output_path = NULL,
                        suffix = NULL,
                        base_size = NULL,
                        add_annotation = TRUE){


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



  RES$significance <- calculate_significance_categories(p = p,
                                                        p_adj = padj,
                                                        fc = FC,
                                                        thres_fc = thres_fc,
                                                        thres_p = thres_p)

  ### transform p-values and fold changes and thresholds
  ### TODO: was wenn p-Wert exakt 0 ist?
  RES$transformed_FC = log(RES[, columnname_FC], base = log_base_fc) # default: log2(FC)
  RES$transformed_p = -log(RES[, columnname_p], base = log_base_p) # default: -log10(p)

  log_thres_p <- -log(thres_p, base = log_base_p)
  log_thres_fc <- log(thres_fc, base = log_base_fc)



  plot <- ggplot2::ggplot(data = RES, ggplot2::aes(x = transformed_FC, y = transformed_p, colour = significance)) +
    ggplot2::geom_point(alpha = 5/10) +
    ggplot2::scale_colour_manual(values = c("not significant" = colour1, "significant" = colour2, "significant after FDR correction" = colour3), drop = FALSE) +
    ### TODO: axis labels with expressions
    ggplot2::xlab(paste0("log",log_base_fc,"(FC)")) +
    ggplot2::ylab(paste0("-log",log_base_p,"(p)"))


  if (!is.null(base_size)) {
    plot <- plot + ggplot2::theme_bw(base_size = base_size)
  } else {
    plot <- plot + ggplot2::theme_bw()
  }

  # TODO: legend_position einstellbar machen
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


  ## draw threshold lines
  plot <- plot + ggplot2::geom_hline(yintercept = log_thres_p, linetype = "dotted") +
    ggplot2::geom_vline(xintercept =  log_thres_fc, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = -log_thres_fc, linetype = "dotted")


  #openxlsx::write.xlsx(x = RES, file=paste0(output_path,"transformedData", suffix, ".xlsx"), overwrite = TRUE,keepNA = TRUE)
  ggplot2::ggsave(paste0(output_path,"Volcano_Plot", suffix, ".",plot_device),
         plot = plot, device = plot_device,
         height = plot_height, width = plot_width, dpi = plot_dpi, units = "cm")
  return(list(RES = RES, plot = plot))
}



################################################################################
################################################################################
################################################################################


### which to label:
## significant after FDR correction
## significant
## top 10 on each side?




#' Add labels to a volcano plot
#'
#' @param RES_Volcano  result from volcanoPlot(): a list containing the data frame with the transformed data and the ggplot object
#' @param label_type "FDR" (significant after correction) or "noFDR" (significant without correction)
#' @param protein_name_column column name of the protein names in the RES data frame
#'
#' @return ggplot object with labels
#' @export
#'
#' @examples # TODO
add_labels <- function(RES_Volcano,
                       label_type = "FDR",
                       protein_name_column = "Gene.names") {


  if (label_type == "FDR") {
    ind_label <- which(RES_Volcano$RES$significance == "significant after FDR correction")
    if (length(ind_label) == 0) {
      warning("No significant proteins after FDR correction for labelling. Try changing label_type to noFDR or lower the fold change threshold.")
      return(RES_Volcano$plot)
    }
  }
  if (label_type == "noFDR") {
    ind_label <- which(RES_Volcano$RES$significance == "significant")
    if (length(ind_label) == 0) {
      warning("No significant proteins for labelling. Try lowering the fold change threshold.")
      return(RES_Volcano$plot)
    }
  }


  ind_label_down <- ind_label[RES_Volcano$RES$transformed_FC[ind_label] < 0]
  ind_label_up <- ind_label[RES_Volcano$RES$transformed_FC[ind_label] > 0]

  labels_up <- rep(NA, nrow(RES_Volcano$RES))
  labels_up[ind_label_up] <- RES_Volcano$RES[, protein_name_column][ind_label_up]
  labels_down <- rep(NA, nrow(RES_Volcano$RES))
  labels_down[ind_label_down] <- RES_Volcano$RES[, protein_name_column][ind_label_down]

  ### TODO: what if genenames are too long?
  ### TODO: what if genenames are not unique?
  ### TODO: what if genenames are not available?

  nudge_x = 0.2
  nudge_y = 0.2


  # change axis limits to have more space for labels
  xaxis_limits <- ggplot2::layer_scales(RES_Volcano$plot)$x$get_limits()
  xaxis_limits <- xaxis_limits * 1.1

  yaxis_limits <- ggplot2::layer_scales(RES_Volcano$plot)$y$get_limits()
  yaxis_limits[2] <- yaxis_limits[2] * 1.1 # only changge upper limit


  ### add labels
  plot <- RES_Volcano$plot +
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









# RES_Volcano: output from the Volcano_Standard function

#
# VolcanoPlot_labelled <- function(RES_Volcano, protein_name_column = "Gene.name",
#                                  columnname_p = "p",columnname_padj = "padj",columnname_FC = "FC",
#                                  thres_p = 0.05, thres_fc = 2,
#                                  plot_height = 15, plot_width = 15, plot_dpi = 300, plot_device="pdf",
#                                  output_path = NULL, suffix = NULL) {
#
#
#   fc <- RES_Volcano$RES[,columnname_FC]
#   p <- RES_Volcano$RES[,columnname_p]
#   padj <- RES_Volcano$RES[,columnname_padj]
#
#
#   ### which proteins to label
#   ind_label <- which(padj < thres_p & !is.na(padj) & (fc < 1/thres_fc | fc > thres_fc))
#   #D_signi <- D[ind_signi,]
#
#
#   ind_label_down <- ind_label[fc[ind_label] < 1]
#   ind_label_up <- ind_label[fc[ind_label] > 1]
#
#
#   labels_up <- rep(NA, nrow(RES_Volcano$RES))
#   labels_up[ind_label_up] <- RES_Volcano$RES[, protein_name_column][ind_label_up]
#
#
#   labels_down <- rep(NA, nrow(RES_Volcano$RES))
#   labels_down[ind_label_down] <- RES_Volcano$RES[, protein_name_column][ind_label_down]
#
#
#   nudge_x = 0.2
#   nudge_y = 0.2
#
#
#   xaxis_limits <- layer_scales(RES_Volcano$plot)$x$get_limits()
#   print(xaxis_limits)
#   xaxis_limits <- xaxis_limits * 1.3
#   print(xaxis_limits)
#
#
#   yaxis_limits <- layer_scales(RES_Volcano$plot)$y$get_limits()
#   print(yaxis_limits)
#   yaxis_limits[2] <- yaxis_limits[2] * 1.2
#   print(yaxis_limits)
#
#
#
#
#   ### add labels
#   plot <- RES_Volcano$plot +
#     ylim(yaxis_limits) + xlim(xaxis_limits) +
#    # geom_point(data= RES_Volcano$RES[!is.na(labels_up) | !is.na(labels_down),],
#   #             aes(x=transformed_FC,y=transformed_p)) +
#     ggrepel::geom_label_repel(aes(label = labels_up), size = 2,
#                               nudge_x = nudge_x, nudge_y = nudge_y, show.legend = FALSE,
#                               max.overlaps = 20) +
#     ggrepel::geom_label_repel(aes(label = labels_down), size = 2,
#                               nudge_x = -nudge_x, nudge_y = nudge_y, show.legend = FALSE,
#                               max.overlaps = 20)
#
#   ggsave(paste0(output_path,"Volcano_Plot_labelled", suffix, ".",plot_device),
#          plot = plot, device = plot_device,
#          height = plot_height, width = plot_width, dpi = plot_dpi, units = "cm")
#
# }
#









