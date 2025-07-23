

#' Create a Heatmap.
#'
#' @param D                       \strong{data.frame} \cr
#'                                The data set containing only protein intensities, already filtered for interesting candidates.
#' @param id                      \strong{data.frame} \cr
#'                                The corresponding ID columns for the parameter D e.g. containing further columns like protein or gene names
#' @param protein_names_col       \strong{integer} \cr
#'                                The column with protein or gene names, if the names should be plotted.
#' @param na_method               \strong{character} \cr
#'                                The method with which missing values are handeled.
#'                                Options are "na.omit" (proteins with any missing values will be removed), "impute" (missing values will be imputed) and "keep" (missing values will be kept).
#'                                Note that clustering may not work when too many missing values are present.
#' @param filtermissings          \strong{integer} \cr
#'                                The threshold for missing values.
#'                                If a protein has more missing values, it will be filtered out.
#'                                Note that rows with only 1 or 2 valid values may cause problems with clustering.
#' @param groups                  \strong{character factor} \cr
#'                                The group of the data with two or more grouping variables.
#' @param group_colours           \strong{character} \cr
#'                                ? A named list of group colours (discrete vars) or colour functions (continuous vars)
#' @param column_split            \strong{character} \cr
#'                                ? The name of the column in groups, if the columns should be split.
#' @param cluster_column_slices   \strong{logical} \cr
#'                                ? If \code{TRUE}, column slices will be clustered.
#' @param cluster_rows            \strong{logical or dendrogram} \cr
#'                                If \code{TRUE}, the rows will be clustered. Can also be a dendrogram object which is used to cluster the rows.
#' @param cluster_cols            \strong{logical or dendrogram} \cr
#'                                If \code{TRUE}, the columns will be clustered. Can also be a dendrogram object which is used to cluster the columns.
#' @param dist_method             \strong{character} \cr
#'                                The distance metric for clustering. Options are "pearson", "spearman" and "euclidean".
#' @param clust_method            \strong{character} \cr
#'                                The linkage method for clustering. Options are "complete", "single" and "average".
#' @param symmetric_legend        \strong{logical} \cr
#'                                If \code{TRUE}, the colour code will be  symmetric.
#'                                Note that it only make sense for z-scored data.
#' @param scale_data              \strong{logical} \cr
#'                                If \code{TRUE}, the data will be scaled ( = z-scored).
#' @param legend_name             \strong{character} \cr
#'                                The name for legend.
#' @param title                   \strong{character} \cr
#'                                The title of the plot.
#' @param legend_colours          \strong{character} \cr
#'                                A vector of colours for colour gradient.
#' @param log_data                \strong{logical} \cr
#'                                If \code{TRUE}, the data will be log-transformed.
#' @param log_base                \strong{integer} \cr
#'                                The base for the log-transformation, if \code{log_data = TRUE}.
#' @param colour_scale_max        \strong{numeric} \cr
#'                                The cap value for which all greater values will receive the same color.
#' @param textsize                \strong{integer} \cr
#'                                The size of text in the plot.
#' @param ...                     Further arguments to Heatmap
#'
#' @return A heatmap of the given data.
#'
#' @export
#' @examples
#'

Heatmap_with_groups <- function(D,
                                id,
                                protein_names_col = NULL,
                                na_method = "na.omit",
                                filtermissings = 2,
                                groups = NULL,
                                group_colours = NULL,
                                column_split = NULL,
                                cluster_column_slices = FALSE,
                                cluster_rows = TRUE,
                                cluster_cols = TRUE,
                                dist_method = "pearson",
                                clust_method = "complete",
                                symmetric_legend = TRUE,
                                scale_data = TRUE,
                                #output_path = paste0(getwd(), "//"),
                                #suffix = NULL,
                                legend_name = "Legend",
                                title = "Heatmap",
                                legend_colours = c("blue", "white", "red"),
                                #plot_height = 20,
                                #plot_width = 20,
                                #plot_dpi = 300,
                                log_data = TRUE,
                                log_base = 2,
                                colour_scale_max = NULL,
                                textsize = 15,
                                top_annotation = NULL,
                                ...) {

  data.asmatrix <- as.matrix(D)

  ### log-transformation
  if (log_data) {
    data.asmatrix <- log(data.asmatrix, log_base)
  }

  ### calculation of z-scores
  if (scale_data) {
  data.asmatrix_scaled <- t(scale(t(data.asmatrix)))
  data.asmatrix <- data.asmatrix_scaled
  }


  ### remove rows with only missing values
  ind  <- rowSums(!is.na(data.asmatrix)) >= filtermissings
  data.asmatrix <- data.asmatrix[ind,]
  id <- id[ind,, drop = FALSE]


  ### cap colour gradient at a maximum value (may be valuable if there are extreme outliers)
  if (!is.null(colour_scale_max)) {
    data.asmatrix[data.asmatrix < -colour_scale_max] <- -colour_scale_max
    data.asmatrix[data.asmatrix > colour_scale_max] <- colour_scale_max
  }


  ### imputation or filter out rows with missing values
  if (na_method == "impute") {
    data.asmatrix[is.na(data.asmatrix)] <- 0
  }

  if (na_method == "na.omit") {
    #data.asmatrix_tmp <- data.asmatrix
    #rownames(data.asmatrix_tmp) <- 1:nrow(data.asmatrix_tmp)
    data.asmatrix <- stats::na.omit(data.asmatrix)
    ind <- as.numeric(rownames(data.asmatrix))

    ### TODO: geht kaputt weil keine rownames?
    #id <- id[ind,, drop = FALSE]
    #print(id)

    ### if there are no rows remaining after na.omit, throw error message
    if (nrow(data.asmatrix) == 0) {
      stop("All rows contain at least one missing value. Heatmap cannot be created.")
    }

  }

  if (na_method == "keep") {
    data.asmatrix <- data.asmatrix
  }




  ### make legend/colour gradient symmetric
  if (symmetric_legend) {
    minmax <- max(abs(c(min(data.asmatrix, na.rm = TRUE), max(data.asmatrix, na.rm = TRUE))))
    legend_colours <- circlize::colorRamp2(c(-minmax, 0, minmax), legend_colours)
  }

  ### get row labels from id dataframe
  if (!is.null(protein_names_col)) {
    row_labels <- id[, protein_names_col]
    #row_labels <- unlist(as.vector(row_labels))
    #print(row_labels)
  } else {
    row_labels <- rep("", nrow(data.asmatrix))
  }


  if (is.null(top_annotation)) {
    ### top annotation for groups
    if (!is.null(groups)) {
      top_annotation = ComplexHeatmap::HeatmapAnnotation(Group = groups,
                                                         col = group_colours,
                                                         annotation_name_gp = grid::gpar(fontsize = textsize), name = "Group",
                                                         annotation_legend_param = list(title_gp = grid::gpar(fontsize = textsize,
                                                                                                              fontface = "bold"),
                                                                                        labels_gp = grid::gpar(fontsize = textsize)))
    } else {
      top_annotation = NULL
    }
  }

  ### set column split
  if (!is.null(column_split)) {
    column_split = groups[, column_split]
  }


  ### set up hierarchical clustering
  if (is.logical(cluster_rows) && cluster_rows == TRUE) {
    cluster_rows = stats::as.dendrogram(stats::hclust(amap::Dist(data.asmatrix, method = dist_method), method = clust_method))
  }
  if (is.logical(cluster_cols) && cluster_cols == TRUE) {
    cluster_cols = stats::as.dendrogram(stats::hclust(amap::Dist(t(data.asmatrix), method = dist_method), method = clust_method))
  }

  #row_labels <<- row_labels
  #data.asmatrix <<- data.asmatrix

  print(row_labels)

  #row.names(data.asmatrix) <- row_labels


  ht <- ComplexHeatmap::Heatmap(data.asmatrix,
                column_title = title,
                #name = legend_name,
                cluster_rows = cluster_rows,
                cluster_columns = cluster_cols,
                cluster_column_slices = cluster_column_slices,
                top_annotation = top_annotation,
                column_split = column_split,
                row_labels = row_labels,
                col = legend_colours,
                heatmap_legend_param = list(direction = "vertical", title = legend_name, title_gp = grid::gpar(fontsize = textsize, fontface = "bold"),
                                            labels_gp = grid::gpar(fontsize = textsize)),
                column_names_gp = grid::gpar(fontsize = textsize),
                column_title_gp = grid::gpar(fontsize = textsize),
                ...)

  return(ht)
  #return(list("heatmap" = ht, "data_as_matrix" = cbind(id, data.asmatrix)))
}






