#' Clustering, Heatmap and Lineplots
#'
#' @param D dataframe with log-transformed protein intensities, e.g. filtered for significant proteins form the ANOVA or t-test results
#' @param id dataframe with id information e.g. protein names, gene names, accessions etc.
#' @param output_path path to save results
#' @param suffix suffix for file name
#' @param nr_clusters number of clusters. Default is NULL, meaning that the optimal number of clusters will be determined by dendextend::find_k().
#'
#' @return save heatmap and data frame with cluster information, as well as line plots
#' @export
#'
#' @examples # TODO
Clustering_heatmap_lineplots <- function(D, id, output_path, suffix = "", nr_clusters = NULL) {

  rownames(D) <- 1:nrow(D)  # reset rownames (important to match cluster information later)
  row_dend = stats::as.dendrogram(stats::hclust(amap::Dist(D, method = "pearson"))) # cluster the proteins

  ## generate global heatmap
  ## TODO: is this necessary?
  # ProtStatsWF::Heatmap_with_groups(D = D,
  #                                  id = id,
  #                                  #groups = timepoint_tmp,
  #                                  cluster_rows = TRUE,
  #                                  cluster_columns = FALSE,
  #                                  output_path = output_path,
  #                                  suffix = suffix)


  if (is.null(nr_clusters)) {
  # find optimal number of clusters based on silhouette values
  nr_clusters <- dendextend::find_k(row_dend)$k
  }

  ## define colours for each cluster
  cluster_colours <- scales::hue_pal()(nr_clusters)
  ## colour branches of the dendrogram to plot next to the heatmap
  row_dend_color = dendextend::color_branches(row_dend, k = nr_clusters, col = cluster_colours)

  ht <- ProtStatsWF::Heatmap_with_groups(D = D,
                                         id = id,
                                         # TODO: no filtering at the moment but it may be necessary/useful depending on the data
                                         filtermissings = ncol(D)+1,
                                         cluster_columns = FALSE,
                                         cluster_rows = row_dend_color,
                                         output_path = output_path,
                                         suffix = paste(suffix, nr_clusters, sep = "_"),
                                         ### TODO: allow omitting rows with missing values
                                         na_method = "impute",
                                         row_split = nr_clusters,
                                         row_gap = grid::unit(5, "mm"))

  ### get cluster for each protein (cluster number from heatmap doesn't correspond to apply cutree() on the dendrogram. This is why we need to get the cluster number from the heatmap directly).
  ht = ComplexHeatmap::draw(ht$Heatmap)
  x <- ComplexHeatmap::row_dend(ht)
  cluster <- integer(nrow(D_norm_tmp_signi))
  for (j in 1:nr_clusters) {
    cluster_members <- as.integer(names(stats::cutree(x[[j]],1))) ### get cluster members
    cluster[cluster_members] <- j
  }


  ### write table with cluster results
  RES_clustering <- cbind(id, cluster = cluster, D)
  openxlsx::write.xlsx(RES_clustering, paste(output_path, suffix, "_", nr_clusters, ".xlsx", sep = "_"))


  ###############
  ### draw lineplot for each cluster, coloured by distance to the cluster center

  D_zscore <- ht$data_wo_imputation
  id_columns <- 1:ncol(id)
  ### TODO: currently only plots without imputation. This should be changed to allow for imputation as well.

  grDevices::pdf(paste0(output_path, "/Lineplots_", suffix, "_", nr_clusters, ".pdf"), width = 10, height = 10)

  for (i in 1:nr_clusters) {

    ## choose only data points from the specific cluster
    D_tmp <- D_zscore[cluster == i,-id_columns]

    ## calculate mean profile of the cluster
    mean_profile <- colMeans(D_tmp, na.rm = TRUE)

    ## calculate euclidean distance of each protein to the cluster center
    Dists_euclidean <- apply(D_tmp, 1, function(x) stats::dist(rbind(x, mean_profile)))

    X <- data.frame(D_tmp, Dists_euclidean, id = 1:nrow(D_tmp))
    X_long <- reshape2::melt(X, id.vars = c("id", "Dists_euclidean"))

    X_long <- rbind(X_long, data.frame(id = max(X_long$id)+1, Dists_euclidean = NA,
                                       variable = colnames(D_tmp),
                                       value = mean_profile))
    X_long <- dplyr::mutate(X_long, ClusterCenter = dplyr::case_when(is.na(Dists_euclidean) ~ "Cluster Center", TRUE ~ "Cluster Members"))



    pl <- ggplot2::ggplot(data = X_long, ggplot2::aes(x = variable, y = value, group = id,
                                                      colour = Dists_euclidean, linetype = ClusterCenter, size = ClusterCenter)) +
      ggplot2::geom_line() +
      ggplot2::scale_colour_gradient2(low = "red", mid = "yellow", high = "green", na.value = "black",
                                      midpoint = 2, limits = c(0, max(X_long$Dists_euclidean, na.rm = TRUE)), name = "Distance \nto center") +
      ggplot2::scale_linetype_manual(values=c("dotted", "solid"), na.value = "solid", name = "") +
      ggplot2::scale_size_manual(values=c(1.3, 1), na.value = 1, guide = "none") +
      ggplot2::xlab("") + ggplot2::ylab("Z-Score") +
      ggplot2::scale_x_discrete(expand = c(0.03, 0.03)) +
      ggplot2::ggtitle(paste0("Cluster ", i, " (", nrow(X), " proteins)")) +
      ggplot2::theme_bw(base_size = 20) +
      ggplot2::theme(legend.key.width = ggplot2:unit(1.5,"cm"), axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggplot2::guides(linetype = ggplot2::guide_legend(override.aes = list(linewidth = 1.3)))+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, colour = cluster_colours[i]), legend.position = "bottom")

    print(pl)


  }

  grDevices::dev.off()

}



