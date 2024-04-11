#' Calculate a PCA plot for proteomics data
#'
#' @param D              A data.frame of the data set.
#' @param groupvar1      A variable used for colors.
#' @param groupvar2      A variable used for shapes.
#' @param impute         If \code{TRUE}, missing values will be imputed.
#' @param impute_method  A character containing the imputation method ("mean" or "median")
#' @param propNA         A numeric of the proportion of allowed missing NAs for a protein, before it is discarded.
#' @param scale.         If \code{TRUE}, the data will be scaled before computing the PCA.
#' @param PCx            The principle component for the x-axis (default: 1).
#' @param PCy            The principle component for the y-axis (default: 2).
#' @param groupvar1_name Titles of legends for colour and shape.
#' @param groupvar2_name Titles of legends for colour and shape.
#' @param group_colours  The colors for the groups.
#' @param alpha          If \code{TRUE}, the data points will be transparent.
#' @param label          If \code{TRUE}, the samples will be labeled.
#' @param label_seed     A numeric, which sets the seed for the label.
#' @param label_size     A numeric containing the size of the sample labels.
#' @param xlim           Limit of the x-axis.
#' @param ylim           Limit of the y-axis.
#' @param point.size     The size of the data points.
#' @param base_size      The base size for the plot.
#' @param ...            Additional arguments for ggplot2::plot.
#'
#' @return The PCA plot for a proteomics sample
#' @export
#'
#' @examples
#' \dontrun{
#' prepared_data <- prepareData(...)
#'
#' pca <- PCA_Plot(D = prepared_data[["D"]])
#'}
#'

PCA_Plot <- function(D,
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

  filtered_D <- filter_PCA_data(D = D, impute = impute, impute_method = impute_method, propNA = propNA)

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
  pl <- ggplot2::ggplot(data = D_PCA, ggplot2::aes(x = PCx, y = PCy), text = paste("sample:", "label")) +
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

  return(list("plot" = pl, "D_PCA_plot" = cbind(D_PCA, "Sample" = colnames(D)),
              "pca" = pca, "message" = mess))
}

