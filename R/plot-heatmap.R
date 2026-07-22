

#' Create a Heatmap.
#'
#' @param D                       \strong{data.frame} \cr
#'                                The data set containing only protein intensities, already filtered for interesting candidates.
#' @param id                      \strong{data.frame} \cr
#'                                The corresponding ID columns for the parameter D e.g. containing further columns like protein or gene names
#' @param proteinNamesCol        \strong{integer} \cr
#'                                The column with protein or gene names, if the names should be plotted.
#' @param naMethod               \strong{character} \cr
#'                                The method with which missing values are handeled.
#'                                Options are "na.omit" (proteins with any missing values will be removed), "impute" (missing values will be imputed) and "keep" (missing values will be kept).
#'                                Note that clustering may not work when too many missing values are present.
#' @param minValidValues         \strong{integer} \cr
#'                                The minimum number of valid values \string{per row}.
#'                                If a protein has less valid values, it will be filtered out.
#'                                Note that rows with only 1 or 2 valid values may cause problems with clustering.
#' @param groups                  \strong{character factor} \cr
#'                                The group of the data with two or more grouping variables.
#' @param groupColours            \strong{character} \cr
#'                                ? A named list of group colours (discrete vars) or colour functions (continuous vars)
#' @param columnSplit             \strong{character} \cr
#'                                ? The name of the column in groups, if the columns should be split.
#' @param clusterColumnSlices    \strong{logical} \cr
#'                                ? If \code{TRUE}, column slices will be clustered.
#' @param clusterRows             \strong{logical or dendrogram} \cr
#'                                If \code{TRUE}, the rows will be clustered. Can also be a dendrogram object which is used to cluster the rows.
#' @param clusterColumns          \strong{logical or dendrogram} \cr
#'                                If \code{TRUE}, the columns will be clustered. Can also be a dendrogram object which is used to cluster the columns.
#' @param distMethod              \strong{character} \cr
#'                                The distance metric for clustering. Options are "pearson", "spearman" and "euclidean".
#' @param clustMethod             \strong{character} \cr
#'                                The linkage method for clustering. Options are "complete", "single" and "average".
#' @param symmetricLegend         \strong{logical} \cr
#'                                If \code{TRUE}, the colour code will be  symmetric.
#'                                Note that it only make sense for z-scored data.
#' @param scaleData               \strong{logical} \cr
#'                                If \code{TRUE}, the data will be scaled ( = z-scored).
#' @param legendName              \strong{character} \cr
#'                                The name for legend (colour bar).
#' @param groupName               \strong{character} \cr
#'                                The name of the group variable, only used for the legend.
#' @param title                   \strong{character} \cr
#'                                The title of the plot.
#' @param legendColours           \strong{character} \cr
#'                                A vector of colours for colour gradient.
#' @param logData                 \strong{logical} \cr
#'                                If \code{TRUE}, the data will be log-transformed.
#' @param logBase                 \strong{integer} \cr
#'                                The base for the log-transformation, if \code{logData = TRUE}.
#' @param colourScaleMax          \strong{numeric} \cr
#'                                The cap value for which all greater values will receive the same color.
#' @param textSize                \strong{integer} \cr
#'                                The size of text in the plot.
#' @param verbose                 \strong{logical} \cr
#'                                If \code{TRUE}, messages will be printed.
#' @param ...                     Further arguments to Heatmap
#'
#' @return A heatmap of the given data.
#'
#' @export
#' @examples
#'

Heatmap_with_groups <- function(D,
                                id,
                                proteinNamesCol = NULL,
                                naMethod = "na.omit",
                                minValidValues = 2,
                                groups = NULL,
                                groupColours = NULL,
                                columnSplit = NULL,
                                clusterColumnSlices = FALSE,
                                clusterRows = TRUE,
                                clusterColumns = TRUE,
                                distMethod = "pearson",
                                clustMethod = "complete",
                                symmetricLegend = TRUE,
                                scaleData = TRUE,
                                #output_path = paste0(getwd(), "//"),
                                #suffix = NULL,
                                legendName = "Legend",
                                groupName = "Group",
                                title = "Heatmap",
                                legendColours = c("blue", "white", "red"),
                                #plot_height = 20,
                                #plot_width = 20,
                                #plot_dpi = 300,
                                logData = TRUE,
                                logBase = 2,
                                colourScaleMax = NULL,
                                textSize = 15,
                                top_annotation = NULL,
                                verbose = TRUE,
                                ...) {

  data.asmatrix <- as.matrix(D)

  ### log-transformation
  if (logData) {
    data.asmatrix <- log(data.asmatrix, logBase)
  }

  ### calculation of z-scores
  if (scaleData) {
  data.asmatrix_scaled <- t(scale(t(data.asmatrix)))
  data.asmatrix <- data.asmatrix_scaled
  }


  ### remove rows with too many missing values
  ind_row  <- rowSums(!is.na(data.asmatrix)) >= minValidValues
  data.asmatrix <- data.asmatrix[ind_row,]
  id <- id[ind_row,, drop = FALSE]
  if (verbose) message(sum(!ind_row), " rows with too many missing values removed.")

  ### cap colour gradient at a maximum value (may be valuable if there are extreme outliers)
  if (!is.null(colourScaleMax)) {
    data.asmatrix[data.asmatrix < -colourScaleMax] <- -colourScaleMax
    data.asmatrix[data.asmatrix > colourScaleMax] <- colourScaleMax
  }


  ### imputation or filter out rows with missing values
  if (naMethod == "impute") {
    data.asmatrix[is.na(data.asmatrix)] <- 0
  }

  if (naMethod == "na.omit") {
    #data.asmatrix_tmp <- data.asmatrix
    #rownames(data.asmatrix_tmp) <- 1:nrow(data.asmatrix_tmp)
    data.asmatrix <- stats::na.omit(data.asmatrix)
    # ind <- as.numeric(rownames(data.asmatrix))

    #id <- id[ind,, drop = FALSE]
    #print(id)

    ### if there are no rows remaining after na.omit, throw error message
    if (nrow(data.asmatrix) == 0) {
      message("All rows contain at least one missing value. Heatmap cannot be created.")
      return(NULL)
    }
    if (nrow(data.asmatrix) == 1) {
      message("Only one row remaining after removing rows with missing values. Heatmap cannot be created.")
      return(NULL)
    }

  }

  if (naMethod == "keep") {
    data.asmatrix <- data.asmatrix
  }




  ### make legend/colour gradient symmetric
  if (symmetricLegend) {
    minmax <- max(abs(c(min(data.asmatrix, na.rm = TRUE), max(data.asmatrix, na.rm = TRUE))))
    legendColours <- circlize::colorRamp2(c(-minmax, 0, minmax), legendColours)
  }

  ### get row labels from id dataframe
  if (!is.null(proteinNamesCol)) {
    row_labels <- id[, proteinNamesCol]
    #row_labels <- unlist(as.vector(row_labels))
    #print(row_labels)
  } else {
    row_labels <- rep("", nrow(data.asmatrix))
  }


  if (is.null(top_annotation)) {
    ### top annotation for groups
    if (!is.null(groups)) {
      top_annotation = ComplexHeatmap::HeatmapAnnotation(Group = groups,
                                                         col = groupColours,
                                                         annotation_name_gp = grid::gpar(fontsize = textSize), #name = groupName,
                                                         annotation_label = groupName,
                                                         annotation_legend_param = list(title_gp = grid::gpar(fontsize = textSize,
                                                                                                              fontface = "bold"),
                                                                                        labels_gp = grid::gpar(fontsize = textSize)))
    } else {
      top_annotation = NULL
    }
  }

  ### set column split
  if (!is.null(columnSplit)) {
    columnSplit = groups[, columnSplit]
  }


  X <<- data.asmatrix

  ### set up hierarchical clustering
  if (is.logical(clusterRows) && clusterRows == TRUE) {
    clusters_rows <- stats::hclust(amap::Dist(data.asmatrix, method = distMethod), method = clustMethod)
    clusterRows <- stats::as.dendrogram(clusters_rows)
  }
  if (is.logical(clusterColumns) && clusterColumns == TRUE) {
    clusters_columns <- stats::hclust(amap::Dist(t(data.asmatrix), method = distMethod), method = clustMethod)
    clusterColumns <- stats::as.dendrogram(clusters_columns)
  }

  #row_labels <<- row_labels
  #data.asmatrix <<- data.asmatrix

  #print(row_labels)

  #row.names(data.asmatrix) <- row_labels

  data.asmatrix2 <<- data.asmatrix

  ht <- ComplexHeatmap::Heatmap(data.asmatrix,
                column_title = title,
                #name = legendName,
                cluster_rows = clusterRows,
                cluster_columns = clusterColumns,
                cluster_column_slices = clusterColumnSlices,
                top_annotation = top_annotation,
                column_split = columnSplit,
                row_labels = row_labels,
                col = legendColours,
                heatmap_legend_param = list(direction = "vertical", title = legendName, title_gp = grid::gpar(fontsize = textSize, fontface = "bold"),
                                            labels_gp = grid::gpar(fontsize = textSize)),
                column_names_gp = grid::gpar(fontsize = textSize),
                column_title_gp = grid::gpar(fontsize = textSize),
                ...)

  return(list(heatmap = ht))#, clusters_columns = clusters_columns)) # clusters_rows = clusters_rows,
  #return(list("heatmap" = ht, "data_as_matrix" = cbind(id, data.asmatrix)))
}






