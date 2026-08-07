


#' Paired Wilcoxon-test for a single row of a data set.
#'
#' @param x                         \strong{numeric vector} \cr
#'                                  The abundances of the data.
#' @param group                     \strong{character factor} \cr
#'                                  The group membership of the data.
#' @param sample                    \strong{character factor} \cr
#'                                  The sample membership of the data.
#' @param delogForFC              \strong{logical} \cr
#'                                  If \code{TRUE}, the fold change will be calculated on the original scale.
#' @param minNrPairs              \strong{integer} \cr
#'                                  The minimum number of complete sample pairs.
#' @param logBase                  \strong{numeric} \cr
#'                                  The base of the logarithm for the log-transformation.
#' @param row                       \strong{integer} \cr
#'                                  The row number of the data for the function call.
#'
#' @return A vector with the following components: mean difference between groups, test statistics, p-value, free space fpr corrected p-value, fold changes (both directions), lower and upper limit of confidence interval, number of valid values per group.
#'
wilcoxon_single_row_paired <- function(x, group, sample, delogForFC = TRUE,
                             minNrPairs = NULL,
                             logBase = 2, row = NULL) {

  ## throw error if the two groups do not have the same length
  if (table(group)[1] != table(group)[2]) {
    stop("Groups don't have the same size, which is required for a paired t-test.")
  }

  ## log-transformation
  x <- unname(x)
  abundance = x

  tmp <- data.frame(abundance = abundance, group = group, sample = sample)
  groupnames <- levels(droplevels(group))

  if (length(groupnames) > 2) stop("More than 2 groups present, please use (repeated measurement) ANOVA.")

  ### column names for results:
  Y <- c("test_statistic",
         "p", "p.fdr",
         paste0("n_", groupnames[1]), paste0("n_", groupnames[2]))
  res <- rep(NA, 5)
  names(res) <- Y

  tmp_group1 <- tmp[tmp$group == groupnames[1],]
  tmp_group1 <- tmp_group1[order(tmp_group1$sample),]
  tmp_group2 <- tmp[tmp$group == groupnames[2],]
  tmp_group2 <- tmp_group2[order(tmp_group2$sample),]

  if (any(tmp_group1$sample != tmp_group2$sample)) {
    stop("Different samples present for both groups")
  }

  ## calculate differences between samples:
  diffs <- tmp_group1$abundance - tmp_group2$abundance
  ind_complete_pairs <- which(!is.na(diffs))

  if (sum(!is.na(diffs)) < minNrPairs) {return(res)}  # ensure that enough complete pairs are present

  ttest <- try({stats::wilcox.test(y = tmp_group1$abundance[ind_complete_pairs],  ## not possible to use tmp_na.omit as pairs may be shifted by omitting NAs
                       x = tmp_group2$abundance[ind_complete_pairs],
                       paired = TRUE)}, silent = TRUE)

  # it is still possible, that the ttest fails (e.g. if variance in one group is 0)
  if ("try-error" %in% class(ttest)) {warning(paste0("ttest failed for row ", row));return(res)}

  res[1] <- ttest$statistic
  res[2] <- ttest$p.value
  res[3] <- NA # free space for corrected p-value
  res[4] <- sum(!is.na(tmp_group1$abundance))
  res[5] <- sum(!is.na(tmp_group2$abundance))

  return(res)
}





################################################################################
################################################################################
################################################################################


#' Function to compute t-test (paired or unpaired)
#' @param D                         **data.frame** \cr
#'                                  The data set containing only protein intensities.
#' @param id                        **data.frame** \cr
#'                                  Data fram containing ID columns like protein names, gene names etc.
#' @param group                     **factor** \cr
#'                                  The group membership of the data.
#' @param sample                    **factor** \cr
#'                                  The sample numbers or IDs of the samples (only used to build pairs if paired = TRUE, otherwise it is ignored).
#' @param paired                    **logical** \cr
#'                                  If TRUE, the test will be paired, otherwise it will be unpaired. Default is FALSE.
#' @param varEqual                 **logical(1)** \cr
#'                                  If TRUE, the variances of the groups are expected to be equal. Default is FALSE.
#' @param logBeforeTest           **logical** \cr
#'                                  If TRUE, the data will be log-transformed before the test. Default is TRUE. Set it to FALSE if data is already log-transformed, e.g.
#'                                  if you have already normalized it during [prepareDataSE()].
#' @param delogForFC              **logical** \cr
#'                                  If TRUE, the fold change will be calculated on the original scale. Default is TRUE.
#' @param logBase                  **numeric** \cr
#'                                  The base of the logarithm for the log-transformation. Default is 2.
#' @param minObsPerGroup         **integer** \cr
#'                                  The minimum number of observations per group. For a paired ttest, this is counted as the minimum of complete pairs. Default is 3.
#' @param minObsPerGroupRatio   **numeric** \cr
#'                                  The minimum number of observations per group as a ratio (e.g, 0.8 = 80% valid values in each group needed).
#'
#' @return A data frame containing the results of the t-test.
#' @export
#'
#' @seealso [ttest_single_row()], [ttest_single_row_paired()]
#'
#' @examples
#'

wilcoxtest <- function(D, id = NULL, group, sample = NULL, paired = FALSE, varEqual = FALSE,
                  delogForFC = TRUE, logBase = 2,
                  minObsPerGroup = 3,
                  minObsPerGroupRatio = NULL) {

  if (paired) {
    RES <- pbapply::pbapply(D, 1, wilcoxon_single_row_paired, group = group,
                 delogForFC = delogForFC, minNrPairs = minObsPerGroup, sample = sample,
                 logBase = logBase)
  }

  RES <- t(RES)
  RES <- as.data.frame(RES)

  RES$p.fdr <- stats::p.adjust(RES$p, method = "fdr")

  ## sort columns in D by group and sample
  if(!is.null(sample)) {
    D <- D[,order(group, sample)]
  } else {
    D <- D[,order(group)]
  }

  if(!is.null(id)) {

    D <- cbind(id, D)

  }


  return(cbind(D, RES))
}




