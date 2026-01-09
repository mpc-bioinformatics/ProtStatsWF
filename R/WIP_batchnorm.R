

# X <- DATA_TMT_Peptides[, 3:20]
# batch <- factor(rep(1:3, each = 6))
# log <-  TRUE
#
# library(limma)
#
# norm <- regression_normalization(X, batch, log = TRUE)


## TODO: Fortschrittsbalken einbauen!
## TODO: log-transform in den workflow auslagern
## TODO: id should be also an input and later cbinded to the normalized data
## TODO: output also data after prenorm and before batchnorm


#### X: Daten
#### batch: Faktor mit Batch-Informationen (codiert als factor)
#### log: sollen Daten vorher logarithmiert werden?
#### pre_norm: soll vorher jeder Batch einzeln mit loess normalisiert werden?
#### return_coefs: Soll data.frame mit Regressionskoeffizienten ausgegeben werden?

#' Regression normalization to remove batch effects
#'
#' @param D
#' @param batch
#' @param pre_norm
#' @param pre_norm_method
#' @param return_coefs
#'
#' @returns
#' @export
#'
#' @examples
regression_normalization <- function(D,
                                     batch,
                                     #log = TRUE,
                                     pre_norm = TRUE,
                                     pre_norm_method = "loess",
                                     return_coefs = FALSE) {
  require(limma)

  #### order data by batch
  ord <- order(batch)
  D <- D[, ord]
  batch <- batch[ord]

  #if(log) D <- log2(D)

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
  D_batchnorm <- D_batchnorm[, order(ord)]

  #colnames(D_coef) <- levels(batch)

  if (return_coefs) {
    return(list(D_batchnorm = D_batchnorm, D_coef = D_coef, D_batchwise_norm = D_batchwise_norm))
  } else {
    return(list(D_batchnorm = D_batchnorm, D_coef = NULL, D_batchwise_norm = D_batchwise_norm))
  }
}

