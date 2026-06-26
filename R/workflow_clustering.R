

#' Workflow that clusters proteins for similar patterns over the samples and
#' produces a heatmap and lineplots
#'
#' @param data_path              \strong{character} \cr
#'                               The path to an .xlsx file containing the input data.
#' @param output_path            \strong{character} \cr
#'                               The path to the output folder.
#' @param intensity_columns      \strong{integer vector} \cr
#'                               The numbers of the intensity columns in the table.
#' @param nr_clusters \strong{integer(1)} \cr
#'        Number of clusters to cut the dendrogram into. If \code{NULL} (default),
#'        the optimal number of clusters is determined based on silhouette values
#'        using the \code{\link[dendextend]{find_k}} function.
#' @param cluster_colours \strong{character vector} \cr
#'       Colours to use for the different clusters. If \code{NULL} (default),
#'       the default ggplot color palette is used.
#' @param row_split \strong{logical(1)} \cr If TRUE, there will be space between row clusters in the heatmap.
#' @param dist_method \strong{character(1)} \cr
#'        Distance measure to use for the hierarchical clustering. In principle,
#'        all methods available in \code{\link[amap]{Dist}} are possible, however
#'        correlation-based metrics like "correlation",  "pearson" or "spearman"
#'        are recommended. The default is "correlation", which uses the centere
#'        Pearson correlation.
#' @param colour_dend \strong{logical(1)} \cr
#'       If \code{TRUE} (default), the branches of the dendrogram are coloured
#'       according to the clusters, using the defined cluster_colours.
#'
#' @param suffix \strong{character} \cr
#'         The suffix for the output files. It needs to start with an underscore.
#' @param plot_height_heatmap \strong{numeric} \cr
#'                               The height for the heatmap in cm.
#' @param plot_width_heatmap \strong{numeric} \cr
#'                               The width for the heatmap in cm.
#' @param plot_height_lineplot \strong{numeric} \cr
#'                               The height for the lineplots in cm.
#' @param plot_width_lineplot \strong{numeric} \cr
#'                               The width for the lineplots in cm.
#' @param plot_dpi \strong{numeric} \cr
#'                              The plot resolution for the heatmap.
#' @param column_name_protein \strong{character(1)} \cr
#'                              The name of the column containing the protein identifiers.
#' @param ... Additional parameters passed to \code{\link[ProtStatsWF]{Heatmap_with_groups}}.
#'
#' @returns Nothing, but saves a heatmap, a set of lineplots (one per cluster)
#' and a cluster table to the output folder.
#' @export
#'
#' @examples
workflow_clustering <- function(data_path,
                                output_path,
                                intensity_columns,

                                nr_clusters = NULL,
                                cluster_colours = NULL,
                                row_split = TRUE,
                                dist_method = "correlation",
                                colour_dend = TRUE,

                                suffix = "",
                                plot_height_heatmap = 15,
                                plot_width_heatmap = 15,
                                plot_height_lineplot= 20,
                                plot_width_lineplot = 25,
                                plot_dpi = 300,

                                column_name_protein = "Protein",
                                ...) {

  ### TODO: What to do with NAs
  ### TODO: what to do with constant rows (may happen for extremely low of high abundant proteins)

  ### TODO: option to aggregate data by group before clustering

  #### Prepare Data ####
  dataPrep <- prepareTtestData(data_path = data_path , intensity_columns = intensity_columns,
                               remove_missings = TRUE)

  Dprep2id <<- dataPrep$id

  clust <- clustering(dataPrep$D,
             dist_method = dist_method,
             nr_clusters = nr_clusters,
             cluster_colours = cluster_colours,
             colour_dend = colour_dend)


  ht <- ProtStatsWF::Heatmap_with_groups(D = dataPrep$D,
                                         id = dataPrep$id,
                                         # TODO: no filtering at the moment but it may be necessary/useful depending on the data
                                         #filtermissings = ncol(D),
                                         cluster_rows = clust$row_dend,
                                         cluster_columns = FALSE,
                                         log_data = FALSE,
                                         ### TODO: allow omitting rows with missing values
                                         na_method = "impute",
                                         row_split = clust$nr_clusters,
                                         row_gap = grid::unit(5, "mm"),
                                         ...)

  grDevices::png(paste0(output_path, "/heatmap", suffix, "_", clust$nr_clusters, ".png"),
                 height = plot_height_heatmap,
                 width = plot_width_heatmap, units = "cm", res = plot_dpi)
  graphics::plot(ht$heatmap)
  grDevices::dev.off()


  clusterInfo <- getClusterInfos(heatmap = ht$heatmap, nr_clusters = clust$nr_clusters, D = dataPrep$D, id = dataPrep$ID)
  openxlsx::write.xlsx(clusterInfo, paste0(output_path, "/cluster_table", suffix, "_", clust$nr_clusters, ".xlsx"))

  D_zscore <- data.frame(ht$heatmap@matrix, cluster = clusterInfo$cluster)

  D_zscore2 <<- D_zscore

  lineplots <- Lineplots(D_zscore = D_zscore, cluster_colours = clust$cluster_colours)

  grDevices::pdf(paste0(output_path, "/Lineplots", suffix, "_", clust$nr_clusters, ".pdf"),
                 width = plot_width_lineplot/2.54,
                 height = plot_height_lineplot/2.54)

  for(i in 1:clust$nr_clusters) {
    print(lineplots[i])
  }

  grDevices::dev.off()

  return(invisible(NULL))

}







