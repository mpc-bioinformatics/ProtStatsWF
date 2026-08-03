


#' Batch normalization via linear regression
#'
#' @param SE
#' @param BatchColumn
#' @param pre_norm
#' @param pre_norm_method
#' @param return_coefs
#'
#' @returns
#' @export
#'
#' @examples
regression_normalization <- function(SE,
                                     BatchColumn = "Batch",
                                     pre_norm = TRUE,
                                     pre_norm_method = "loess",  ## TODO
                                     return_coefs = FALSE) {

  D <- as.data.frame(SummarizedExperiment::assay(SE, "intensity_norm"))
  batch <- factor(SummarizedExperiment::colData(SE)[[BatchColumn]])




  #### order data by batch
  ord <- order(batch)
  D <- D[, ord]
  batch <- batch[ord]

  if (pre_norm) {
    D_split <- split.default(D, batch)
    D_split_norm <- lapply(D_split, limma::normalizeBetweenArrays, method = "cyclicloess") ### TODO: Normalisierungsmethode auswählbar machen
    D_batchwise_norm <- do.call(cbind, D_split_norm)
  } else {
    D_batchwise_norm <- D
  }

  res <- rep(NA, nrow(D_batchwise_norm))
  p <- rep(NA, nrow(D_batchwise_norm))
  D_batchnorm <- data.frame()
  D_coef <- list()#matrix(nrow = nrow(DATA_batchwise_norm), ncol = length(levels(batch)))

  for (i in 1:nrow(D_batchwise_norm)) {

    y <- as.numeric(D_batchwise_norm[i,])
    if (all(is.na(y))) {D_batchnorm <- rbind(D_batchnorm, y);D_coef[[i]] <- rep(NA, length(levels(batch)));next} # wenn alles NA ist
    # falls für das Protein nur 1 Batch valide Werte hat, bleiben die Werte so:
    if (length(levels(droplevels(batch[!is.na(y)]))) <= 1) {D_batchnorm <- rbind(D_batchnorm, y);D_coef[[i]] <- rep(NA, length(levels(batch)));next}

    ind_noNA <- which(!is.na(y))
    y_tmp <- y[ind_noNA]
    batch_tmp <- droplevels(batch[ind_noNA])  # immer mit 1,2,3,... codiert

    mod <- lm(y_tmp ~ batch_tmp , contrasts = list(batch_tmp = contr.sum))
    coefs <- mod$coefficients[-1]
    coefs <- c(coefs, -sum(coefs))
    names(coefs) <- levels(batch_tmp)
    if (return_coefs) D_coef[[i]] <- coefs
    ### TODO: gibt irgendwie einen Fehler wenn dieses Proteinen nicht in allen batches vorkommt
    ### (z.B. Zeile 11 beim Normalisieren der 5 Batches vom Sepsis-Projekt)

    ynorm <- rep(NA, length(y))
    ynorm_ <- y_tmp - unname(coefs[as.numeric(batch_tmp)])
    ynorm[ind_noNA] <- ynorm_

    D_batchnorm <- rbind(D_batchnorm, ynorm)
  }


  colnames(D_batchnorm) <- colnames(D)
  D_batchnorm <- D_batchnorm[, order(ord)]  # restore previous order

  colnames(D_batchwise_norm) <- colnames(D)
  D_batchwise_norm <- D_batchwise_norm[, order(ord)] # restore previous order

  rownames(D_batchnorm) <- rownames(SE)
  rownames(D_batchwise_norm) <- rownames(SE)

  SummarizedExperiment::assay(SE, "intensity_batchwisenorm") <- D_batchwise_norm
  SummarizedExperiment::assay(SE, "intensity_batchnorm") <- D_batchnorm


  if (return_coefs) {
    return(list(SE_batchnorm = SE, D_coef = D_coef))
  } else {
    return(list(SE_batchnorm = SE, D_coef = NULL))
  }
}
