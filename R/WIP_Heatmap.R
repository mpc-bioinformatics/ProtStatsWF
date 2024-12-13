################################################################################
#### Heatmap script using the ComplexHeatmap R package


### TODO: einbauen, dass er die Gruppenfarben richtig macht! (evtl. in getrennter Funktion?)



# D: data set with only intensities (should be log-transformed if necessary)
# id: data set with further columns, e.g. protein or gene names
# protein_names_col: column with protein or gene names (NULL if no protein names should be plotted)
# na.method: "na.omit" -> proteins with any missing values will be removed
#            "impute" -> imputation of missing values
#            "keep" -> keep missing values
# filtermissings: filter out proteins with more than X missing values
#                 (rows with only 1 or 2 valid values may cause problems with clustering)

# groups: data.frame with one or more grouping variables
# group_colours: named list of group colours (discrete vars) or colour functions (continuous vars)
# column_split = should columns be split? NULL if not, name of column in groups if yes
# cluster_column_slices: cluster the column slices?

# cluster_rows, cluster_columns: if TRUE, rows/columns will be clustered
# dist_method: distance metric for clustering, e.g. "pearson", "spearman", "euclidean"
# clust_method: linkage method for clustering, e.g. "complete", "single", "average"

# symmetric_legend: should colour code be made symmetric? (only make sense for z-scored data)
# scale_data: should data be scaled ( = z-scored)?
# output_path: path for exporting the plot
# suffix: suffix for file name

# legend_name: name for legend
# title: title
# legend_colours
# plot_height, plot_width, plot_dpi: setting for plot output (height/width in cm)
# ... further arguments to Heatmap






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
#' @param groups                  \strong{integer} \cr
#'                                ? A factor of group or data.frame with two or more grouping variables.
#' @param group_colours           named list of group colours (discrete vars) or colour functions (continuous vars)
#' @param column_split            \strong{character} \cr
#'                                ? The name of the column in groups, if the columns should be split.
#' @param cluster_column_slices   \strong{logical} \cr
#'                                ? If \code{TRUE}, column slices will be clustered.
#' @param cluster_rows            \strong{logical} \cr
#'                                If \code{TRUE}, the rows will be clustered.
#' @param cluster_columns         \strong{logical} \cr
#'                                If \code{TRUE}, the columns will be clustered.
#' @param dist_method             \strong{character} \cr
#'                                The distance metric for clustering. Options are "pearson", "spearman" and "euclidean".
#' @param clust_method            \strong{character} \cr
#'                                The linkage method for clustering. Options are "complete", "single" and "average".
#' @param symmetric_legend        \strong{logical} \cr
#'                                If \code{TRUE}, the colour code will be  symmetric.
#'                                Note that it only make sense for z-scored data.
#' @param scale_data              \strong{logical} \cr
#'                                If \code{TRUE}, the data will be scaled ( = z-scored).
#' @param output_path             \strong{character} \cr
#'                                The path to the output folder.
#' @param suffix                  \strong{character} \cr
#'                                The suffix of the file names should have one.
#' @param legend_name             \strong{character} \cr
#'                                The name for legend.
#' @param title                   \strong{character} \cr
#'                                The title of the plot.
#' @param legend_colours          \strong{character} \cr
#'                                A vector of colours for colour gradient.
#' @param plot_height             \strong{numeric} \cr
#'                                The plot height in cm.
#' @param plot_width              \strong{numeric} \cr
#'                                The plot width in cm.
#' @param plot_dpi                \strong{integer} \cr
#'                                The "dots per inch" of the plot aka. the plot resolution.
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
#' @export
#'
#' @examples 
#' 

Heatmap_with_groups <- function(D, id, protein_names_col = NULL,
                           na_method = "na.omit", filtermissings = 2,
                           groups = NULL, group_colours = NULL,
                           column_split = NULL, cluster_column_slices = FALSE,
                           cluster_rows = TRUE, cluster_columns = TRUE,
                           dist_method = "pearson", clust_method = "complete",
                           symmetric_legend = TRUE, scale_data = TRUE,
                           output_path = paste0(getwd(), "//"), suffix = NULL,
                           legend_name = "Legend", title = "Heatmap",
                           legend_colours = c("blue", "white", "red"),
                           plot_height = 20, plot_width = 20, plot_dpi = 300,
                           log_data = TRUE, log_base = 2,
                           colour_scale_max = NULL, textsize = 15,
                           ...){

  data.asmatrix <- as.matrix(D)

  if(log_data){
    data.asmatrix <- log(data.asmatrix, log_base)
  }



  if (scale_data) {
  ### calculation of z-scores
  data.asmatrix_scaled <- t(scale(t(data.asmatrix)))
  data.asmatrix <- data.asmatrix_scaled
  }

  ### remove rows with only missing values
  ind  <- rowSums(!is.na(data.asmatrix)) >= filtermissings
  data.asmatrix <- data.asmatrix[ind,]
  id <- id[ind,, drop = FALSE]

  ### imputation or filter out rows with missing values
  if(na_method == "impute"){
    data.asmatrix[is.na(data.asmatrix)] <- 0
  }
  if(na_method == "na.omit") {
    data.asmatrix_tmp <- data.asmatrix
    rownames(data.asmatrix_tmp) <- 1:nrow(data.asmatrix_tmp)
    data.asmatrix_tmp <- stats::na.omit(data.asmatrix_tmp)
    ind <- as.numeric(rownames(data.asmatrix_tmp))

    id <- id[ind,, drop = FALSE]

    ### if there are no rows remaining after na.omit, keep them and impute by 0
    if(nrow(data.asmatrix_tmp) == 0){
      print("All rows contain at least one missing value. Switched to imputation instead.")
      data.asmatrix[is.na(data.asmatrix)] <- 0
    } else {
      data.asmatrix <- data.asmatrix_tmp
    }

  }
  if (na_method == "keep") {
    data.asmatrix <- data.asmatrix
  }



  if (!is.null(colour_scale_max)) {
    data.asmatrix[data.asmatrix < -colour_scale_max] <- -colour_scale_max
    data.asmatrix[data.asmatrix > colour_scale_max] <- colour_scale_max

  }



  ### export data used in heatmap
  # openxlsx::write.xlsx(cbind(id, zscore = data.asmatrix), paste0(output_path,"Heatmap_data_", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)

  if (symmetric_legend) {
  minmax <- max(abs(c(min(data.asmatrix, na.rm = TRUE), max(data.asmatrix, na.rm = TRUE))))
  legend_colours <- circlize::colorRamp2(c(-minmax, 0, minmax), legend_colours)
  }


  if (!is.null(protein_names_col)) {
    row_labels <- id[, protein_names_col]
  } else {
    row_labels <- rep("", nrow(data.asmatrix))
  }

  ## to annotation for groups
  if (!is.null(groups)) {
    top_annotation = ComplexHeatmap::HeatmapAnnotation(Group = groups, col = group_colours, annotation_name_gp = grid::gpar(fontsize = textsize), name = "Group",
                                       annotation_legend_param = list(title_gp = grid::gpar(fontsize = textsize, fontface = "bold"),
                                                                      labels_gp = grid::gpar(fontsize = textsize)))
  } else {
    top_annotation = NULL
  }


  if (!is.null(column_split)) {
    column_split = groups[, column_split]
  }


  ht <- ComplexHeatmap::Heatmap(data.asmatrix,
                column_title = title , name = legend_name,
                cluster_rows = cluster_rows,
                clustering_distance_rows = dist_method, clustering_method_rows = clust_method,
                cluster_columns = cluster_columns,
                clustering_distance_columns = dist_method, clustering_method_columns = clust_method,
                cluster_column_slices = cluster_column_slices,
                top_annotation = top_annotation,
                column_split = column_split,
                row_labels = row_labels,
                col = legend_colours,
                heatmap_legend_param = list(direction = "vertical", title = "Legend", title_gp = grid::gpar(fontsize = textsize, fontface = "bold"),
                                            labels_gp = grid::gpar(fontsize = textsize)),
                column_names_gp = grid::gpar(fontsize = textsize),
                column_title_gp = grid::gpar(fontsize = textsize),
                ...)

  #png(paste0(output_path,"Heatmap_", suffix, ".png"), width = plot_width, height = plot_height, units = "cm", res = plot_dpi)
  #ComplexHeatmap::draw(ht, annotation_legend_side = "right", heatmap_legend_side = "right", merge_legend = TRUE)
  #dev.off()

  return(list("heatmap" = ht, "data_as_matrix" = data.asmatrix))
}






