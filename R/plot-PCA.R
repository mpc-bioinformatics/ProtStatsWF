#' PCA plot for proteomics data.
#'
#' @param SE  **SummarizedExperiment object** \cr A SummarizedExperiment object containing the proteomics data (output element SE from [prepareDataSE]). The assay specified by the `assay` parameter will be used for PCA.
#' @param groupForColour  **character(1)** \cr The name of the column in the colData of the SummarizedExperiment object to be used for colouring the points in the PCA plot.
#' @param colourType **character(1)** \cr Type for colour, either "discrete" (default) or "continuous".
#' @param groupForShape  **character(1)** \cr The name of the column in the colData of the SummarizedExperiment object to be used for shaping the points in the PCA plot.
#' @param assay  **character(1)** \cr The name of the assay in the SummarizedExperiment object to be used for PCA.
#' @param imputeMethod  **character(1)** \cr The method to use for imputing missing values. Options are "mean" (default) and "median".
#' @param propNA  **numeric** \cr The proportion of missing values allowed before imputation is applied.
#' @param scale. **logical(1)** \cr Whether to scale data for the PCA. Default is TRUE.
#' @param PCx **integer(1)** \cr The principal component to be plotted on the x-axis. Default is 1.
#' @param PCy **integer(1)** \cr The principal component to be plotted on the y-axis. Default is 2.
#' @param groupColours **character** \cr Vector of colours for the groups. Default is NULL
#' (default ggplot2 colour palette will be used.)
#' @param alpha **numeric(1)** \cr The transparency level (between 0 and 1) of the points. Default is 1 (no transparency).
#' @param label **logical(1)** \cr Whether to label the points. Default is FALSE.
#' @param labelSeed **numeric(1)** \cr Seed for random number generator used for label placement.
#' Default is NA, which may lead to different label placement for multiple executions of this function.
#' @param labelSize **numeric(1)** \cr Size of the labels. Default is 4.
#' @param xlim **numeric(2)** \cr Limits for the x-axis. Default is NULL (automatic).
#' @param ylim **numeric(2)** \cr Limits for the y-axis. Default is NULL (automatic).
#' @param pointSize **numeric(1)** \cr Size of data points in the PCA. Default is 4.
#' @param baseSize **numeric(1)** \cr Base size for the plots. Default is 15.
#' @param NAValueColour **character(1)** \cr Colour for data points with missing values in the groupForColour variable. Default is grey.
#' @param NAValueShape **integer(1)** \cr Shape for data points with missing values in the groupForShape variable. Default is 0 (hollow square). Please see [graphics::pch()] for further information.
#' @param verbose **logical(1)** \cr Whether to print messages during the execution of the function. Default is TRUE.
#'
#' @return TODO
#' @export
#'
#' @examples
PCA_Plot <- function(SE,
                     groupForColour = NULL,
                     colourType = "discrete",
                     groupForShape = NULL,
                     assay = "intensity_norm",

                     imputeMethod = "mean",
                     propNA = 0,
                     scale. = TRUE,
                     PCx = 1,
                     PCy = 2,

                     groupColours = NULL,
                     alpha = 1,
                     label = FALSE,
                     labelSeed = NA,
                     labelSize = 4,
                     xlim = NULL,
                     ylim = NULL,

                     pointSize = 4,
                     baseSize = 11,
                     NAValueColour = "grey",
                     NAValueShape = 0,
                     verbose = TRUE
) {

  # warning if propNA = 0 and imputation method is chosen (no imp. will be done)

  filtered_data <- filter_PCA_data(SE = SE,
                                    assay = assay,
                                   imputeMethod = imputeMethod,
                                   propNA = propNA)

  if(nrow(SE) == 0) {
    if (verbose) message("All rows were filtered out. \n Try increasing the proportion of missing NAs allowed or imputing missing values.")
    return(list("plot" = NULL, "D_PCA_plot" = NULL,
                "pca" = NULL))
  }


  D <- SummarizedExperiment::assays(filtered_data)[[assay]]
  id <- SummarizedExperiment::rowData(filtered_data)

  if (!is.null(groupForColour)) {
    group1 <- SummarizedExperiment::colData(filtered_data)[,groupForColour]
    if (colourType == "discrete") {
      group1 <- factor(group1)
    }
  } else {
    group1 <- NULL
  }
  if (!is.null(groupForShape)) {
    group2 <- SummarizedExperiment::colData(filtered_data)[,groupForShape]
    group2 <- factor(group2)
  } else {
    group2 <- NULL
  }



  if (verbose)
    message(nrow(D), " of ", nrow(SE), " rows are used for PCA.")


  #### calculate PCA ####
  #print(D)

  pca <- stats::prcomp(t(D), scale. = scale.)
  pred <- stats::predict(pca, t(D))
  summ <- summary(pca)

  var50 <- which(summ$importance[3,] >= 0.5)[1]
  if (verbose)
    message("50% explained variance is reached with ", var50, " principle components.")

  pl <- ggplot2::ggplot(mapping = ggplot2::aes(x = PCx, y = PCy))

  #### prepare data.frame ####
  # version with colour and shape
  if (!is.null(groupForColour) & !is.null(groupForShape)) {
    D_PCA <<- data.frame(pred[,c(PCx,PCy)], group1 = group1, group2 = group2, label = colnames(D))
    colnames(D_PCA)[1:2] <- c("PCx", "PCy")
    ### more than 6 different shapes will otherwise give an error message:
    pl <- pl + ggplot2::geom_point(data = D_PCA, ggplot2::aes(colour = group1, shape = group2), size = pointSize, alpha = alpha)
  }

  # version with only colour
  if (!is.null(groupForColour) & is.null(groupForShape)) {
    D_PCA <- data.frame(pred[,c(PCx,PCy)], group1 = group1, label = colnames(D))
    colnames(D_PCA)[1:2] <- c("PCx", "PCy")
    pl <- pl + ggplot2::geom_point(data = D_PCA, ggplot2::aes(colour = group1), size = pointSize, alpha = alpha)
  }

  # version with only shape
  if (is.null(groupForColour) & !is.null(groupForShape)) {
    D_PCA <- data.frame(pred[,c(PCx,PCy)], group2 = group2, label = colnames(D))
    colnames(D_PCA)[1:2] <- c("PCx", "PCy")
    pl <- pl + ggplot2::geom_point(data = D_PCA, ggplot2::aes(shape = group2), size = pointSize, alpha = alpha)
  }

  # version without colour or shape
  if (is.null(groupForColour) & is.null(groupForShape)) {
    D_PCA <- data.frame(pred[,c(PCx,PCy)], label = colnames(D))
    colnames(D_PCA)[1:2] <- c("PCx", "PCy")
    pl <- pl + ggplot2::geom_point(data = D_PCA, size = pointSize, alpha = alpha)
  }

  if (!is.null(groupForColour)) pl <- pl + ggplot2::labs(colour = groupForColour)
  if (!is.null(groupForShape)) pl <- pl + ggplot2::labs(shape = groupForShape)

  if (!is.null(groupForColour) & (!is.null(groupColours) | is.numeric(D_PCA$group1)) & colourType == "discrete") {
    pl <- pl + ggplot2::scale_colour_manual(values = groupColours, na.value = NAValueColour)
  } else {
    if (is.numeric(D_PCA$group1)) {
      pl <- pl + ggplot2::scale_colour_continuous(na.value = NAValueColour)
    } else {
      pl <- pl + ggplot2::scale_colour_discrete(na.value = NAValueColour)
    }
  }

  if (label) {
    pl <- pl + ggrepel::geom_text_repel(data = D_PCA, ggplot2::aes(x=PCx, y=PCy, label = label, colour = group1),
                                        size = labelSize, seed = labelSeed) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = ggplot2::aes(label = "")))
  }

  #### add % of explainable variance to the axis label ####
  pl <- pl + ggplot2::theme_bw(base_size = baseSize) +
    ggplot2::xlab(paste0("PC", PCx, " (", round(100*summ$importance[2,PCx], 1), "%)")) +
    ggplot2::ylab(paste0("PC", PCy, " (", round(100*summ$importance[2,PCy], 1), "%)"))

  ### allow more than 6 different shapes
  if (!is.null(group2)) {
    if (nlevels(D_PCA$group2) > 6) {
      pl <- pl + ggplot2::scale_shape_manual(values = 1:nlevels(D$group2), na.value = NAValueShape)
    } else {
      pl <- pl + ggplot2::scale_shape_discrete(na.value = NAValueShape)
    }

  }


  #### add xlim and ylim if specified ####
  if (!is.null(xlim)) pl <- pl + xlim(xlim)
  if (!is.null(ylim)) pl <- pl + ylim(ylim)

  Loadings <- as.data.frame(pca$rotation)
  if (!is.null(id)) {
    Loadings <- cbind(id, Loadings)
  }

  return(list("D_PCA_plot" = cbind(D_PCA, "Sample" = colnames(D)),
              "pca" = pca, "filtered_data" = filtered_data, "loadings" = Loadings, "plot" = pl))
}

