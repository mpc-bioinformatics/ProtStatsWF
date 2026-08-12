#' Unpaired t-test for a single row of a data set.
#'
#' @description
#' Runs an unpaired t-test on a single row of peptide or protein intensity data
#' and returns p-values and fold changes.
#'
#' @param x **numeric** \cr
#' Vector of peptide or protein intensities.
#' @param group **factor** \cr
#' Vector containing the group membership of the samples. Must contain exactly
#' two factor levels and the same length as \code{x}.
#' @param varEqual **logical(1)** \cr
#' If TRUE, equal variances are assumed (Student's t-test). Default is FALSE
#' (Welch's t-test).
#' @param logBeforeTest **logical(1)** \cr
#' If TRUE, the data will be log-transformed before the t-test. Default is FALSE.
#' @param delogForFC **logical(1)** \cr
#' If TRUE, the fold change will be calculated on the original scale.
#' Default is TRUE.
#' @param logBase **numeric(1)** \cr
#' The base of the logarithm for the log-transformation. Default is 2.
#' @param minObsPerGroup **integer(1)** \cr
#' The minimum number of valid values per group to calculate the test for this
#' peptide or protein. If fewer valid values are present, NA is reported as
#' a result. Default is 3.
#' @param minObsPerGroupRatio **numeric(1)** \cr
#' The minimum proportion of valid values required per group (e.g. 0.8 = 80%
#' valid values in each group required). This is an alternative to
#' \code{minObsPerGroup}, especially when groups sizes are unbalanced. If you
#' use this setting, please set \code{minObsPerGroup = NULL}.
#' @param row **integer(1)** \cr
#' The current row number of the data, used for informative warnings.
#' Default is NULL.
#' @param verbose **logical(1)** \cr
#' If TRUE (default), messages are printed out.
#'
#' @return A named numeric vector with the following elements: mean of group 1,
#' mean of group 2, test statistic, p-value, placeholder for FDR-corrected
#' p-value, fold changes in both directions, confidence interval bounds, and
#' number of valid values per group.
#'
#' @seealso [ttest()], [.ttest_single_row_paired()]
#'
#' @examples
#'

.ttest_single_row <- function(x, group, varEqual = FALSE, logBeforeTest = FALSE,
                              delogForFC = TRUE, logBase = 2, minObsPerGroup = 3,
                              minObsPerGroupRatio = NULL, row = NULL, verbose = TRUE) {

  if (!is.null(minObsPerGroup) & !is.null(minObsPerGroupRatio)) {
    stop("Both minObsPerGroup and minObsPerGroupRatio are given, please define only one of them!")
  }



  x <- unlist(unname(x))
  if (logBeforeTest) {
    abundance = log(x, base = logBase)
  } else {
    abundance = x
  }

  tmp <- data.frame(abundance = abundance, group = group)
  groupnames <- levels(droplevels(group))

  if (length(groupnames) > 2) stop("More than 2 groups present, please use ANOVA.")

  ### column names for results:
  Y <- c(paste0("mean_", groupnames[1]),  # 1
         paste0("mean_", groupnames[2]),  # 2
         "test_statistic",                # 3
         "p", "p.fdr",                    # 4, 5
         paste0("FC_", groupnames[2], "_divided_by_", groupnames[1]),  # 6
         paste0("FC_", groupnames[1], "_divided_by_", groupnames[2]),  # 7
         "CI_lower", "CI_upper",                                       # 8, 9
         paste0("n_", groupnames[1]), paste0("n_", groupnames[2]), # 10, 11
         "NA_reason_code")                                               # 12
  res <- rep(NA, 12)
  names(res) <- Y

  ### drop missing values
  tmp_na.omit <- stats::na.omit(tmp)
  tmp_na.omit <- droplevels(tmp_na.omit)


  if(length(table(tmp_na.omit$group)) < 2) { # if less than 2 groups remain
    res[12] <- 2 ### reason: on/off protein with one group missing
    return(res)

  }              # ensure that both groups are still present

  if(!is.null(minObsPerGroup) & any(table(tmp_na.omit$group) < minObsPerGroup)) {  # ensure that each group has enough observations
    res[12] <- 1  ### reason: not enough observations in at least one group
    return(res)
  }
  if(!is.null(minObsPerGroupRatio) & any(table(tmp_na.omit$group)/table(tmp$group) < minObsPerGroupRatio)) {  # ensure that each group has enough observations
    res[12] <- 1  ### reason: not enough observations in at least one group
    return(res)
  }

  ttest <- try({stats::t.test(x = tmp_na.omit$abundance[tmp_na.omit$group == groupnames[2]],
                       y = tmp_na.omit$abundance[tmp_na.omit$group == groupnames[1]],
                       paired = FALSE, var.equal = varEqual)}, silent = TRUE)

  # it is still possible, that the ttest fails (e.g. if variance in one group is 0)
  if ("try-error" %in% class(ttest)) {warning(paste0("ttest failed for row ", row));res[12] <- 3;return(res)}  ### reason: other, e.g. var = 0

  res[1:2] <- ttest$estimate
  res[3] <- ttest$statistic
  res[4] <- ttest$p.value
  res[5] <- NA # free space for corrected p-value
  res[8:9] <- ttest$conf.int
  res[10] <- length(tmp_na.omit$abundance[tmp_na.omit$group == groupnames[1]])
  res[11] <- length(tmp_na.omit$abundance[tmp_na.omit$group == groupnames[2]])


  ### calculate fold changes
  if (delogForFC) {
    x2 <- logBase^tmp$abundance
  } else {
    x2 <- tmp$abundance
  }

  res[6] <- mean(x2[tmp$group == groupnames[2]], na.rm = TRUE) / mean(x2[tmp$group == groupnames[1]], na.rm = TRUE)
  res[7] <- 1/res[6]

  return(res)
}




###################################################################################################
###################################################################################################
###################################################################################################


#' Paired t-test for a single row of a data set.
#'
#' @description
#' Runs a paired t-test on a single row of intensity data and returns summary
#' statistics including fold changes.
#'
#' @param x **numeric** \cr
#' The abundances of the data.
#' @param group **factor** \cr
#' The group membership of the data. Must contain exactly two levels of equal
#' size.
#' @param sample **factor** \cr
#' The sample IDs used to match pairs across groups.
#' @param logBeforeTest **logical(1)** \cr
#' If TRUE, the data will be log-transformed before the test. Default is TRUE.
#' @param delogForFC **logical(1)** \cr
#' If TRUE, the fold change will be calculated on the original scale.
#' Default is TRUE.
#' @param minNrPairs **integer(1)** \cr
#' The minimum number of complete sample pairs required. Default is NULL.
#' @param logBase **numeric(1)** \cr
#' The base of the logarithm for the log-transformation. Default is 2.
#' @param row **integer(1)** \cr
#' The row number of the data, used for informative warnings. Default is NULL.
#'
#' @return A named numeric vector with the following elements: mean difference
#' between groups, test statistic, p-value, placeholder for FDR-corrected
#' p-value, fold changes in both directions, confidence interval bounds, and
#' number of valid values per group.
#'
#' @seealso [ttest()], [.ttest_single_row()]
#'
#' @examples
#'

.ttest_single_row_paired <- function(x, group, sample, logBeforeTest = TRUE, delogForFC = TRUE,
                             minNrPairs = NULL,
                             logBase = 2, row = NULL) {

  ## throw error if the two groups do not have the same length
  if (table(group)[1] != table(group)[2]) {
    stop("Groups don't have the same size, which is required for a paired t-test.")
  }

  ## log-transformation
  x <- unname(x)
  if (logBeforeTest) {
    abundance = log(x, base = logBase)
  } else {
    abundance = x
  }

  tmp <- data.frame(abundance = abundance, group = group, sample = sample)
  groupnames <- levels(droplevels(group))

  if (length(groupnames) > 2) stop("More than 2 groups present, please use (repeated measurement) ANOVA.")

  ### column names for results:
  Y <- c("mean_of_differences",
         "test_statistic",
         "p", "p.fdr",
         paste0("FC_", groupnames[2], "_divided_by_", groupnames[1]),
         paste0("FC_", groupnames[1], "_divided_by_", groupnames[2]),
         "CI_lower", "CI_upper",
         paste0("n_", groupnames[1]), paste0("n_", groupnames[2]))
  res <- rep(NA, 10)
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

  ttest <- try({stats::t.test(y = tmp_group1$abundance[ind_complete_pairs],  ## not possible to use tmp_na.omit as pairs may be shifted by omitting NAs
                       x = tmp_group2$abundance[ind_complete_pairs],
                       paired = TRUE)}, silent = TRUE)

  # it is still possible, that the ttest fails (e.g. if variance in one group is 0)
  if ("try-error" %in% class(ttest)) {warning(paste0("ttest failed for row ", row));return(res)}

  res[1] <- ttest$estimate
  res[2] <- ttest$statistic
  res[3] <- ttest$p.value
  res[4] <- NA # free space for corrected p-value
  res[7:8] <- ttest$conf.int
  res[9] <- sum(!is.na(tmp_group1$abundance))
  res[10] <- sum(!is.na(tmp_group2$abundance))

  ### calculate fold changes
  if (delogForFC) {
    x2 <- logBase^tmp$abundance
  } else {
    x2 <- tmp$abundance
  }

  res[5] <- mean(x2[tmp$group == groupnames[2]][ind_complete_pairs], na.rm = TRUE) /
    mean(x2[tmp$group == groupnames[1]][ind_complete_pairs], na.rm = TRUE)
  res[6] <- 1/res[5]

  return(res)
}





################################################################################
################################################################################
################################################################################


#' Compute t-test (paired or unpaired)
#'
#' @description
#' Applies an unpaired or paired t-test row-wise across a matrix of protein
#' intensities and returns a combined data frame of input data and test results
#' with FDR-adjusted p-values.
#'
#' @param D **data.frame** \cr
#' The data set containing only protein intensities (rows = proteins,
#' columns = samples).
#' @param id **data.frame** \cr
#' Data frame containing ID columns such as protein names or gene names.
#' Default is NULL.
#' @param group **factor** \cr
#' The group membership of the samples. Must contain exactly two levels.
#' @param sample **factor** \cr
#' Sample numbers or IDs used to build pairs when `paired = TRUE`. Ignored
#' when `paired = FALSE`. Default is NULL.
#' @param paired **logical(1)** \cr
#' If TRUE, a paired t-test is performed, otherwise unpaired. Default is FALSE.
#' @param varEqual **logical(1)** \cr
#' If TRUE, equal variances are assumed (Student's t-test). Default is FALSE
#' (Welch's t-test).
#' @param logBeforeTest **logical(1)** \cr
#' If TRUE, the data will be log-transformed before the test. Default is TRUE.
#' Set to FALSE if the data are already log-transformed, e.g. after
#' normalisation in [prepareDataSE()].
#' @param delogForFC **logical(1)** \cr
#' If TRUE, fold changes are calculated on the original (exponentiated) scale.
#' Default is TRUE.
#' @param logBase **numeric(1)** \cr
#' The base of the logarithm for the log-transformation. Default is 2.
#' @param minObsPerGroup **integer(1)** \cr
#' Minimum number of observations per group. For a paired t-test this is the
#' minimum number of complete pairs. Default is 3.
#' @param minObsPerGroupRatio **numeric(1)** \cr
#' Minimum proportion of valid values required per group (e.g. 0.8 = 80%).
#' Default is NULL.
#'
#' @return A data frame combining the (sorted) input intensities with the
#' per-protein t-test results including FDR-adjusted p-values.
#' @export
#'
#' @seealso [.ttest_single_row()], [.ttest_single_row_paired()]
#'
#' @examples
#'

ttest <- function(D, id = NULL, group, sample = NULL, paired = FALSE, varEqual = FALSE,
                  logBeforeTest = TRUE, delogForFC = TRUE, logBase = 2,
                  minObsPerGroup = 3,
                  minObsPerGroupRatio = NULL) {

  if (!paired) {
    RES <- pbapply::pbapply(D, 1, .ttest_single_row, group = group, logBeforeTest = logBeforeTest,
                 delogForFC = delogForFC, minObsPerGroup = minObsPerGroup,
                 logBase = logBase, varEqual = varEqual, minObsPerGroupRatio = minObsPerGroupRatio)
  }

  if (paired) {
    RES <- pbapply::pbapply(D, 1, .ttest_single_row_paired, group = group, logBeforeTest = logBeforeTest,
                 delogForFC = delogForFC, minNrPairs = minObsPerGroup, sample = sample,
                 logBase = logBase)
  }

  RES <- t(RES)
  RES <- as.data.frame(RES)

  RES$p.fdr <- stats::p.adjust(RES$p, method = "fdr")

  if(logBeforeTest){
    D_log <- log(D, base = logBase)
    colnames(D_log) <- paste0(colnames(D), "_log")
  } else {
    if (delogForFC) {
      D_delog <- logBase^D
      colnames(D_delog) <- paste0(colnames(D), "_delog")
    }
  }


  ## sort columns in D by group and sample
  if(!is.null(sample)) {
    D <- D[,order(group, sample)]
  } else {
    D <- D[,order(group)]
  }

  if(!is.null(id)) {
    if(logBeforeTest) {
      D <- cbind(id, D, D_log)
    } else {
      if (delogForFC) {
        D <- cbind(id, D, D_delog)
      } else {
        D <- cbind(id, D)
      }
    }
  }


  return(cbind(D, RES))
}




