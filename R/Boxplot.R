

Boxplots <- function(D_long,
                     log_data = TRUE, log_base = 2,
                     method = "boxplot",
                     groupvar_name = "Group", group_colours = NULL,
                     plot_device = "png",
                     plot_height = 10, plot_width = 15, plot_dpi = 300,
                     base_size = 15){
  
  
  use_groups <- !all(is.na(D_long$group))
  
  #### log-transform data is necessary ####
  
  if(log_data) {
    D_long$value <- log(D_long$value, base = log_base)
  }
  
  name <- value <- group <- X <- NULL
  if (use_groups) {
    pl_boxplot <- ggplot2::ggplot(D_long, ggplot2::aes(x = name, y = value, fill = group)) +
      ggplot2::labs(fill = groupvar_name)
    if (!is.null(group_colours)) pl_boxplot <- pl_boxplot + ggplot2::scale_fill_manual(values = group_colours)
  } else {
    pl_boxplot <- ggplot2::ggplot(X, ggplot2::aes(x = name, y = value))
  }
  
  
  pl_boxplot <- pl_boxplot +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1)) +
    ggplot2::ylab("Log intensity") + ggplot2::xlab("Sample")
  
  if (method == "violinplot"){
    pl_boxplot <- pl_boxplot + ggplot2::geom_violin()
  }
  
  if (method == "boxplot") {
    pl_boxplot <- pl_boxplot + ggplot2::geom_boxplot()
  }
  
  return(pl_boxplot)
}