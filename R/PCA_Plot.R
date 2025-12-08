#' Calculate a PCA plot for proteomics data.
#'
#' @param D                \strong{data.frame} \cr
#'                         The data set containing intensities of the sample.
#' @param id                \strong{data.frame} \cr
#'                         The corresponding ID columns for the parameter D e.g. containing further columns like protein or gene names
#' @param groupvar1        \strong{character vector} \cr
#'                         The variable used for colors.
#' @param groupvar2        \strong{character vector} \cr
#'                         The variable used for shapes.
#' @param impute           \strong{logical} \cr
#'                         If \code{TRUE}, missing values will be imputed.
#' @param impute_method    \strong{character} \cr
#'                         The imputation method. Options are "mean" or "median".
#' @param propNA           \strong{numeric} \cr
#'                         The proportion of allowed missing NAs for a protein, before it is discarded.
#' @param scale.           \strong{logical} \cr
#'                         If \code{TRUE}, the data will be scaled before computing the PCA.
#' @param PCx              \strong{integer} \cr
#'                         The principle component for the x-axis.
#' @param PCy              \strong{integer} \cr
#'                         The principle component for the y-axis.
#' @param groupvar1_name   \strong{character} \cr
#'                         The titles of legends for colour.
#' @param groupvar2_name   \strong{character} \cr
#'                         The titles of legends for shape.
#' @param group_colours    \strong{character vector} \cr
#'                         The colors for the groups.
#' @param alpha            \strong{logical} \cr
#'                         If \code{TRUE}, the data points will be transparent.
#' @param label            \strong{logical} \cr
#'                         If \code{TRUE}, the samples will be labeled.
#' @param label_seed       \strong{numeric} \cr
#'                         A numeric, which sets the seed for the label.
#' @param label_size       \strong{numeric} \cr
#'                         A numeric containing the size of the sample labels.
#' @param xlim             \strong{numeric} \cr
#'                         The limit of the x-axis.
#' @param ylim             \strong{numeric} \cr
#'                         The limit of the y-axis.
#' @param point.size       \strong{numeric} \cr
#'                         The size of the data points.
#' @param base_size        \strong{numeric} \cr
#'                         The base size for the plot.
#' @param ...              Additional arguments for ggplot2::plot.
#'
#' @return The PCA plot for a proteomics sample.
#' @export
#'
#' @examples
#'

PCA_Plot <- function(D,
                     id = NULL,
                     groupvar1 = NULL, groupvar2 = NULL,

                     impute = FALSE, impute_method = "mean", propNA = 0,
                     scale. = TRUE,
                     PCx = 1, PCy = 2,

                     groupvar1_name = "group", groupvar2_name = NULL,
                     group_colours = NULL, alpha = 1,
                     label = FALSE, label_seed = NA, label_size = 4,
                     xlim = NULL, ylim = NULL,

                     point.size = 4, base_size = 11,
                     ...

) {

  mess = ""

  filtered_data <- filter_PCA_data(D = D, id = id, impute = impute, impute_method = impute_method, propNA = propNA)
  filtered_D <- filtered_data$D
  filtered_id <- filtered_data$id

  if(is.null(filtered_D)){
    mess <- paste0(mess, "All rows were filtered out. \n")
    mess <- paste0(mess, "Try increasing the proportion of missing NAs allowed or imputing missing values. \n")
    message(mess)
    return(list("plot" = NULL, "D_PCA_plot" = NULL,
                "pca" = NULL, "message" = mess))
  }

  mess <- paste0(mess, nrow(filtered_D), " of ", nrow(D), " rows are used for PCA. \n")

  #### calculate PCA ####
  pca <- stats::prcomp(t(filtered_D), scale. = scale.)
  pred <- stats::predict(pca, t(filtered_D))
  summ <- summary(pca)

  var50 <- which(summ$importance[3,] >= 0.5)[1]
  mess <- paste0(mess, paste0("50% explained variance is reached with ", var50, " principle components.\n"))


  #### prepare data.frame ####
  # version with colour and shape
  if (!is.null(groupvar1) & !is.null(groupvar2)) {
    D_PCA <- data.frame(pred[,c(PCx,PCy)], groupvar1 = groupvar1, groupvar2 = groupvar2)
    ### more than 6 different shapes will otherwise give an error message:
    if (nlevels(D_PCA$groupvar2) > 6) pl <- pl + ggplot2::scale_shape_manual(values = 1:nlevels(D$groupvar2))
  }

  # version with only colour
  if (!is.null(groupvar1) & is.null(groupvar2)) {
    D_PCA <- data.frame(pred[,c(PCx,PCy)], groupvar1 = groupvar1, label = colnames(filtered_D))
  }

  # version without colour or shape
  if (is.null(groupvar1) & is.null(groupvar2)) {
    D_PCA <- data.frame(pred[,c(PCx,PCy)], label = colnames(filtered_D))
  }

  colnames(D_PCA)[1:2] <- c("PCx", "PCy")


  #### make plot ####
  pl <- ggplot2::ggplot(data = D_PCA, ggplot2::aes(x = PCx, y = PCy)) + # text = paste("sample:", "label")
    ggplot2::geom_point(ggplot2::aes(colour = groupvar1, shape = groupvar2), size = point.size, alpha = alpha)

  pl <- pl + ggplot2::labs(colour = groupvar1_name, shape = groupvar2_name)

  if(!is.null(groupvar1) & !is.null(group_colours)){
    pl <- pl + ggplot2::scale_colour_manual(values = group_colours)
  }

  if(label) {
    pl <- pl + ggrepel::geom_text_repel(ggplot2::aes(x=PCx, y=PCy, label = label, colour = groupvar1), size = label_size, seed = label_seed) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = ggplot2::aes(label = "")))
  }


  #### add % of explainable variance to the axis label ####
  pl <- pl + ggplot2::theme_bw(base_size = base_size) +
    ggplot2::xlab(paste0("PC", PCx, " (", round(100*summ$importance[2,PCx], 1), "%)")) +
    ggplot2::ylab(paste0("PC", PCy, " (", round(100*summ$importance[2,PCy], 1), "%)"))

  # allow more than 6 different shapes
  if (!is.null(groupvar1) & !is.null(groupvar2)) {
    if (nlevels(D_PCA$groupvar2) > 6) pl <- pl + ggplot2::scale_shape_manual(values = 1:nlevels(D$groupvar2))
  }


  #### add xlim and ylim if specified ####
  if (!is.null(xlim)) pl <- pl + xlim(xlim)
  if (!is.null(ylim)) pl <- pl + ylim(ylim)

  message(mess)

  Loadings <- as.data.frame(pca$rotation)
  if (!is.null(filtered_id)) {
    Loadings <- cbind(filtered_id, Loadings)
  }

  return(list("D_PCA_plot" = cbind(D_PCA, "Sample" = colnames(D)),
              "pca" = pca, "message" = mess, "filtered_D" = filtered_D, "loadings" = Loadings, "plot" = pl))
}

