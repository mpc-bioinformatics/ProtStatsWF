




#' Simple volcano plot from p-value, fold change and significance category
#'
#' @param p vector of p-values (before p-value correction)
#' @param FC vector of fold changes
#' @param significance_category vector of significance categories ("not significant", "significant", "significant after FDR correction")
#' @param log_base_fc base for logarithm of the fold change (default is 2).
#' @param log_base_p base for logarithm of the p-value (default is 10, the negative log-p-value is plotted, e.g. -log10(p)).
#' @param thres_p  threshold for p-value (default is 0.05)
#' @param thres_fc threshold for fold change (default is 2)
#' @param colour1 colour for not significant points, default is grey
#' @param colour2 colour for significant points, default is black
#' @param colour3 colour for significant points after FDR correction, default is orange
#' @param symmetric_x if TRUE, x-axis will be made symmetric, default is FALSE
#' @param legend_position position of the legen, default is "bottom"
#' @param base_size base size for ggplot2 theme, default is NULL
#' @param xlim x-axis limits, default is NULL, i.e. will be determined automatically
#' @param ylim y-axis limits, default is NULL, i.e. will be determined automatically
#'
#' @return Volcano plot as ggplot object
#' @export
#'
#' @examples # TODO
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
                        ylim = NULL) {


  ### transform p-values and fold changes and thresholds
  ### TODO: What if p-value is exactly zero?
  transformed_FC <- log(FC, base = log_base_fc) # default: log2(FC)
  transformed_p <- -log(p, base = log_base_p) # default: -log10(p)

  log_thres_p <- -log(thres_p, base = log_base_p)
  log_thres_fc <- log(thres_fc, base = log_base_fc)

  RES <- data.frame(transformed_FC = transformed_FC,
                    transformed_p = transformed_p,
                    significance = significance_category)


  significance <- NULL # silence notes when checking the package

  plot <- ggplot2::ggplot(data = RES, ggplot2::aes(x = transformed_FC, y = transformed_p, colour = significance)) +
    ggplot2::geom_point(alpha = 5/10) +
    ggplot2::scale_colour_manual(values = c("not significant" = colour1, "significant" = colour2, "significant after FDR correction" = colour3), drop = FALSE) +
    ### TODO: axis labels with expressions
    ### TODO: what if I do not have these categories?
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

  ## draw threshold lines
  plot <- plot + ggplot2::geom_hline(yintercept = log_thres_p, linetype = "dotted") +
    ggplot2::geom_vline(xintercept =  log_thres_fc, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = -log_thres_fc, linetype = "dotted")

  return(plot)

}

################################################################################
################################################################################
################################################################################



#' Calculation of significance categories for a volcano plot for the t-test
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
calculate_significance_categories_ttest <- function(p, p_adj, fc, thres_fc=2, thres_p=0.05) {

  significance <- dplyr::case_when(
    p_adj <= thres_p & p <= thres_p & (fc >= thres_fc | fc <= 1/thres_fc) & !is.na(p) ~ "significant after FDR correction",
    p_adj > thres_p & p <= thres_p & (fc >= thres_fc | fc <= 1/thres_fc) & !is.na(p) ~ "significant",
    (p > thres_p | (fc < thres_fc & fc > 1/thres_fc)) & !is.na(p) ~ "not significant",
    is.na(p) ~ NA_character_
  )

  significance <- factor(significance, levels = c("not significant", "significant", "significant after FDR correction"))

  return(significance)
}




#' Volcano plot for a t-test result
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
VolcanoPlot_ttest <- function(RES,
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



  RES$significance <- calculate_significance_categories_ttest(p = p,
                                                              p_adj = padj,
                                                              fc = FC,
                                                              thres_fc = thres_fc,
                                                              thres_p = thres_p)


  plot <- VolcanoPlot(p = RES$p,
                      FC = RES$FC,
                      significance_category = RES$significance,
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
  ggplot2::ggsave(paste0(output_path,"Volcano_Plot", suffix, ".",plot_device),
         plot = plot, device = plot_device,
         height = plot_height, width = plot_width, dpi = plot_dpi, units = "cm")
  return(plot = plot) # its not necessary to return the data also, they are already present in the plot object
}


################################################################################
################################################################################
################################################################################


#' Calculate significance categories for ANOVA
#'
#' @param p_posthoc vector of posthoc p-values
#' @param p_anova_adj vector of ANOVA p-values adjusted for multiple testing
#' @param p_anova vector of ANOVA p-values before p-value correction
#' @param fc vector of fold changes
#' @param thres_fc threshold for fold change
#' @param thres_p threshold for p-value
#'
#' @return vector of significance categories
#' @export
#'
#' @examples # TODO
calculate_significance_categories_ANOVA <- function(p_posthoc, p_anova_adj, p_anova, fc, thres_fc=2, thres_p=0.05) {

  significance <- dplyr::case_when(
    p_anova_adj <= thres_p & p_posthoc <= thres_p & (fc >= thres_fc | fc <= 1/thres_fc) & !is.na(p_posthoc) & !is.na(p_anova) ~ "significant after FDR correction", # ANOVA significant after FDR, posthoc also significant, fulfills FC threshold
    p_anova_adj > thres_p & p_anova <= thres_p & p_posthoc <= thres_p & (fc >= thres_fc | fc <= 1/thres_fc) & !is.na(p_posthoc) & !is.na(p_anova) ~ "significant", # ANOVA significant before FDR, posthoc also significant, fulfills FC threshold
    (p_anova > thres_p | p_posthoc > thres_p | (fc < thres_fc & fc > 1/thres_fc)) & !is.na(p_posthoc) & !is.na(p_anova) ~ "not significant", # ANOVA not significant or posthoc not significant or FC does not fulfill threshold
    is.na(p_posthoc) | is.na(p_anova) ~ NA_character_
  )

  significance <- factor(significance, levels = c("not significant", "significant", "significant after FDR correction"))

  return(significance)
}





#' Volcano Plots for ANOVA
#'
#' @param RES result of ProtStatsWF::ANOVA() or other ANOVA function
#' @param columnname_p_ANOVA column name of the p-values from the ANOVA
#' @param columnname_p_ANOVA_adj column name of the adjusted p-values from the ANOVA
#' @param columns_FC vector of column numbers for the fold changes
#' @param columns_p_posthoc vector of column numbers for the posthoc p-values
#' @param log_base_fc base for log2 transformation of fold changes
#' @param log_base_p base for log10 transformation of p-values
#' @param thres_p threshold for p-values
#' @param thres_fc threshold for fold changes
#' @param colour1 colour for not significant proteins
#' @param colour2 colour for significant proteins
#' @param colour3 colour for significant proteins after FDR correction
#' @param symmetric_x logical, if TRUE, x-axis will be made symmetric
#' @param legend_position position of the legend
#' @param base_size base size for the ggplot theme
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param add_labels add labels to interesting proteins in the volcano plot
#'
#' @return list of volcano plots as ggplot objects
#' @export
#'
#' @examples # TODO
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
  comp_names <- colnames(RES)[columns_p_posthoc]
  comp_names <- stringr::str_replace_all(comp_names, "p.posthoc.", "")
  comp_names <- stringr::str_replace_all(comp_names, "_", " ")

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




#' Add labels to a volcano plot
#'
#' @param RES_Volcano  result from volcanoPlot(): a list containing the data frame with the transformed data and the ggplot object
#' @param label_type "FDR" (significant after correction) or "noFDR" (significant without correction) or "index" -> define indizes to label
#' @param protein_name_column column name of the protein names in the RES data frame
#' @param ind index of the proteins to label
#' @param protein_names vector of protein names to label
#'
#' @return ggplot object with labels
#' @export
#'
#' @examples # TODO
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

  ind_label_down <- ind_label[RES_Volcano$data$transformed_FC[ind_label] < 0]
  ind_label_up <- ind_label[RES_Volcano$data$transformed_FC[ind_label] > 0]


  labels_up <- rep(NA, nrow(RES_Volcano$data))
  labels_down <- rep(NA, nrow(RES_Volcano$data))
  if (is.null(protein_names)) {
    labels_up[ind_label_up] <- RES_Volcano$data[, protein_name_column][ind_label_up]
    labels_down[ind_label_down] <- RES_Volcano$data[, protein_name_column][ind_label_down]
  } else {
    labels_up[ind_label_up] <- protein_names[ind_label_up]
    labels_down[ind_label_down] <- protein_names[ind_label_down]
  }

  ### TODO: what if genenames are too long?
  ### TODO: what if genenames are not unique?
  ### TODO: what if genenames are not available?

  nudge_x = 0.2
  nudge_y = 0.2


  # change axis limits to have more space for labels
  xaxis_limits <- ggplot2::layer_scales(RES_Volcano)$x$get_limits()
  xaxis_limits <- xaxis_limits * 1.1

  yaxis_limits <- ggplot2::layer_scales(RES_Volcano)$y$get_limits()
  yaxis_limits[2] <- yaxis_limits[2] * 1.1 # only change upper limit

  ### add labels
  plot <- RES_Volcano +
    ggplot2::ylim(yaxis_limits) + ggplot2::xlim(xaxis_limits) +
    ggrepel::geom_label_repel(ggplot2::aes(label = labels_up), size = 2,
                              nudge_x = nudge_x, nudge_y = nudge_y, show.legend = FALSE,
                              max.overlaps = 20) +
    ggrepel::geom_label_repel(ggplot2::aes(label = labels_down), size = 2,
                              nudge_x = -nudge_x, nudge_y = nudge_y, show.legend = FALSE,
                              max.overlaps = 20)

  return(plot)
}



