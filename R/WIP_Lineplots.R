

#### prerequisites:
### data with protein/peptide intensities, ideally already z-scored
### clustering (vector with cluster numbers, e.g. form hclust or extracted from a heatmap)

### id_columns?
### if data needs to be z-scored or not?
### if replicates should be aggregated or not (then groups have to be given??)


lineplots_of_clusters <- function(D, clustering, group = NULL) {

  ### TODO: einbauen, dass Daten hier noch ge-z-scored werden können, falls nötig

  nr_clusters <- length(unique(clustering))
  cluster_names <- sort(unique(clustering))

  cluster_colours <- scales::hue_pal()(nr_clusters)
  ### TODO: map cluster colours to colours in heatmap???



  line_plots <- list()

  for (i in 1:nr_clusters) {

    D_tmp <- D[clustering == cluster_names[i],]
    #D_tmp <- D_tmp[,-ncol(D_tmp)]  ### remove column with the cluster number

    ### TODO: einbauen, dass replicate einer Gruppe zusätzlich zusammengefasst werden!
    mean_profile <- colMeans(D_tmp, na.rm = TRUE)

    ### if group variable is given, merge replicates over group by median
    if (!is.null(group)) {
      mean_profile_tmp <- data.frame(mean_profile, group)
      aggregate(mean_profile ~ g, data = mean_profile_tmp, median)
      mean_profile_tmp_med <- mean_profile_tmp %>% group_by(group) %>% summarise(median = median(mean_profile, na.rm = TRUE), .groups = "keep")
      mean_profile <- mean_profile_tmp_med$median[match(group, mean_profile_tmp_med$group)]
    }

    ### calculate euclidean distance between the individual profiles and the mean profile
    ### TODO: What happens when replicates over group are merged in the mean profile?
    ### TODO: allow for different distance measures
    Dists_euclidean <- apply(D_tmp, 1, function(x) dist(rbind(x, mean_profile)))
    #Dists_spearman <- apply(D_tmp, 1, function(x) as.dist(1-cor(t(rbind(x, mean_profile)), method="spearman")))
    range(Dists_euclidean)


    X <- data.frame(D_tmp, Dists_euclidean, id = 1:nrow(D_tmp))

    X_long <- reshape2::melt(X, id.vars = c("id", "Dists_euclidean"))

    X_long <- rbind(X_long, data.frame(id = max(X_long$id)+1, Dists_euclidean = NA,
                                       variable = colnames(D_tmp),
                                       value = mean_profile))
    X_long <- dplyr::mutate(X_long, ClusterCenter = dplyr::case_when(is.na(Dists_euclidean) ~ "Cluster Center", TRUE ~ "Cluster Members"))



    pl <- ggplot2::ggplot(data = X_long, ggplot2::aes(x = variable, y = value, group = id,
                                    colour = Dists_euclidean, linetype = ClusterCenter, size = ClusterCenter)) +
      ggplot2::geom_line() + #geom_point() +
      ggplot2::scale_colour_gradient2(low = "red", mid = "yellow", high = "green", na.value = "black",
                             midpoint = 2, limits = c(0, max(X_long$Dists_euclidean, na.rm = TRUE)), name = "Distance \nto center") +
      ggplot2::scale_linetype_manual(values = c("dotted", "solid"), na.value = "solid", name = "") +
      ggplot2::scale_size_manual(values = c(1.3, 1), na.value = 1, guide = "none") +
      ggplot2::xlab("") + ggplot2::ylab("Z-Score") +
      ggplot2::scale_x_discrete(expand = c(0.03, 0.03)) +
      ggplot2::ggtitle(paste0("Cluster ", cluster_names[i], " (", nrow(X), " proteins)")) +
      ggplot2::theme_bw(base_size = 20) +
      ggplot2::theme(legend.key.width = ggplot2::unit(1.5,"cm"), axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggplot2::guides(linetype = ggplot2::guide_legend(override.aes = list(linewidth = 1.3), order = 1))+
      #guides(fill=guide_legend("Cluster center", override.aes=list(fill="black"))) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, colour = cluster_colours[i]), legend.position = "bottom")
    pl

    line_plots[[i]] <- pl
    print(pl)
  }

  return(line_plots)

}


