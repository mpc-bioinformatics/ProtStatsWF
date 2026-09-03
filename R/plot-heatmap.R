

#' Create a Heatmap.
#'
#' @param D **data.frame** \cr
#' The data set containing only protein intensities, already filtered for
#' interesting candidates.
#' @param id **data.frame** \cr
#' The corresponding ID columns for the parameter D e.g. containing further
#' columns like protein or gene names.
#' @param proteinNameColumn **character(1)** \cr
#' The name of the column in \code{id} with protein or gene names, if the
#' names should be plotted.
#' @param naMethod **character(1)** \cr
#' The method with which missing values are handled.
#' Options are "na.omit" (default, proteins with any missing values will be
#' removed), "impute" (missing values will be imputed) and "keep" (missing
#' values will be kept). Note that clustering may not work when too many missing
#' values are present.
#' @param minValidValues **integer(1)** \cr
#' The minimum number of valid values **per row**.
#' If a protein has less valid values, it will be filtered out.
#' Note that rows with only 1 or 2 valid values may cause problems with
#' clustering. Default is 2.
#' @param groups **data.frame** \cr
#' Dataframe containing grouping variables, which are used for column annotations
#' of the Heatmap. Default is NULL, meaning no column annotations are used.
#' @param groupColours **list** \cr
#' A named list of group colours (discrete vars) or colour functions
#' (continuous vars), used for the column annotation. Default is NULL, meaning
#' that random colours will be chosen. See [ComplexHeatmap::HeatmapAnnotation]
#' for details.
#' @param columnSplit **character** \cr  ## ????
#' The name of the column in the groups dataframe by which the columns should
#' be split. Default is NULL, then columns will not be splitted.
#' @param clusterColumnSlices **logical (1)** \cr
#' If \code{TRUE}, column slices will be clustered. Default is FALSE.
#' @param clusterRows **logical(1) or dendrogram** \cr
#' If \code{TRUE} (default), the rows will be clustered. Can also be a dendrogram
#' object which is used to cluster the rows.
#' @param clusterColumns **logical or dendrogram** \cr
#' If \code{TRUE} (default), the columns will be clustered. Can also be a dendrogram
#' object which is used to cluster the columns.
#' @param distMethod **character** \cr
#' The distance metric for clustering. Options are "pearson", "spearman" and
#' "euclidean". Default is "pearson".
#' @param clustMethod **character** \cr
#' The linkage method for clustering. Options are "complete", "single" and
#' "average". Default is "complete".
#' @param symmetricLegend **logical** \cr
#' If \code{TRUE} (default), the colour code will be symmetric.
#' Note that it only make sense for z-scored data.
#' @param scaleData **logical** \cr
#' If \code{TRUE} (default), the data will be scaled ( = z-scored).
#' @param legendName **character(1)** \cr
#' The name for legend (colour bar). Default is "Legend".
#' @param groupName **character** \cr
#' The name of the group variable, only used for the legend. Default is "Group".
#' @param title **character** \cr
#' The main title of the plot. Default is "Heatmap".
#' @param legendColours **character** \cr
#' A vector of colours for colour gradient.
#' @param logData **logical** \cr
#' If \code{TRUE}, the data will be log-transformed.
#' @param logBase **integer** \cr
#' The base for the log-transformation, if \code{logData = TRUE}.
#' @param colourScaleMax **numeric** \cr
#' The cap value for which all greater values will receive the same color.
#' @param textSize **integer** \cr
#' The size of text in the plot.
#' @param verbose **logical** \cr
#' If \code{TRUE}, messages will be printed.
#' @param ... Further arguments to Heatmap
#'
#' @return A heatmap of the given data.
#'
#' @export
#'
#' @importFrom checkmate assert assertCharacter assertDataFrame assertFlag
#' @importFrom checkmate assertInt assertList assertNumber assertSubset
#' @importFrom checkmate checkClass checkFlag
#'
#' @examples
#'

Heatmap_with_groups <- function(D,
                                id,
                                proteinNameColumn = NULL,
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
                                legendName = "Legend",
                                groupName = "Group",
                                title = "Heatmap",
                                legendColours = c("blue", "white", "red"),
                                logData = TRUE,
                                logBase = 2,
                                colourScaleMax = NULL,
                                textSize = 15,
                                top_annotation = NULL,
                                verbose = TRUE,
                                ...) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package \"ComplexHeatmap\" must be installed to plot a heatmap.",
         call. = FALSE)
  }


  checkmate::assertDataFrame(D, min.rows = 1, min.cols = 1)
  checkmate::assertDataFrame(id, nrows = nrow(D))
  checkmate::assertCharacter(proteinNameColumn, len = 1, null.ok = TRUE)
  checkmate::assertSubset(naMethod, c("na.omit", "impute", "keep"))
  checkmate::assertInt(minValidValues, lower = 0)
  checkmate::assertDataFrame(groups, nrows = ncol(D), null.ok = TRUE)
  checkmate::assertList(groupColours, null.ok = TRUE)
  checkmate::assertCharacter(columnSplit, len = 1, null.ok = TRUE)
  checkmate::assertFlag(clusterColumnSlices)
  checkmate::assert(checkmate::checkFlag(clusterRows),
                    checkmate::checkClass(clusterRows, "dendrogram"))
  checkmate::assert(checkmate::checkFlag(clusterColumns),
                    checkmate::checkClass(clusterColumns, "dendrogram"))
  checkmate::assertSubset(distMethod, c("pearson", "spearman", "euclidean"))
  checkmate::assertSubset(clustMethod, c("complete", "single", "average"))
  checkmate::assertFlag(symmetricLegend)
  checkmate::assertFlag(scaleData)
  checkmate::assertCharacter(legendName, len = 1)
  checkmate::assertCharacter(groupName, len = 1)
  checkmate::assertCharacter(title, len = 1)
  checkmate::assertCharacter(legendColours)
  checkmate::assertFlag(logData)
  checkmate::assertNumber(logBase, lower = 1)
  checkmate::assertNumber(colourScaleMax, null.ok = TRUE)
  checkmate::assertNumber(textSize, lower = 0)
  checkmate::assertFlag(verbose)

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


  ## TODO: do a manual clustering here. This allows also for more/better distance functions
  ## and control.

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
    data.asmatrix <- stats::na.omit(data.asmatrix)

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
  if (!is.null(proteinNameColumn)) {
    row_labels <- id[, proteinNameColumn]
    #row_labels <- unlist(as.vector(row_labels))
    #print(row_labels)
  } else {
    row_labels <- rep("", nrow(data.asmatrix))
  }


  if (is.null(top_annotation)) {
    ### top annotation for groups
    if (!is.null(groups)) {
      top_annotation = ComplexHeatmap::HeatmapAnnotation(df = groups,
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

  ### set up hierarchical clustering
  if (is.logical(clusterRows) && clusterRows == TRUE) {
    clusters_rows <- stats::hclust(amap::Dist(data.asmatrix, method = distMethod), method = clustMethod)
    clusterRows <- stats::as.dendrogram(clusters_rows)
  }
  if (is.logical(clusterColumns) && clusterColumns == TRUE) {
    clusters_columns <- stats::hclust(amap::Dist(t(data.asmatrix), method = distMethod), method = clustMethod)
    clusterColumns <- stats::as.dendrogram(clusters_columns)
  }

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
}






