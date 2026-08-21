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
#' @importFrom checkmate assertNumeric assertFactor assertFlag assertNumber
#' @importFrom checkmate assertIntegerish assertInt
#' @importFrom stats na.omit t.test
.ttest_single_row <- function(x, group, varEqual = FALSE, logBeforeTest = FALSE,
                              delogForFC = TRUE, logBase = 2,
                              minObs = 3, minObsRatio = NULL,
                              row = NULL, verbose = TRUE) {

  checkmate::assertNumeric(x)
  checkmate::assertFactor(group, len = length(x), n.levels = 2)
  checkmate::assertFlag(varEqual)
  checkmate::assertFlag(logBeforeTest)
  checkmate::assertFlag(delogForFC)
  checkmate::assertNumber(logBase, lower = 1)
  checkmate::assertIntegerish(minObs, null.ok = TRUE, lower = 2)
  checkmate::assertNumber(minObsRatio, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assertInt(row, null.ok = TRUE)
  checkmate::assertFlag(verbose)

  if (!is.null(minObs) & !is.null(minObsRatio)) {
    stop("Both minObs and minObsRatio are given, please define only one of them!")
  }

  x <- unlist(unname(x))
  if (logBeforeTest) intensity <- log(x, base = logBase) else intensity <- x
  tmp <- data.frame(intensity = intensity, group = group)
  groupnames <- levels(droplevels(group))
  if (length(groupnames) > 2) stop("More than 2 groups present, please use ANOVA.")

  ### column names for results:
  cnames <- c(paste0("mean_", groupnames[1]),  # 1
         paste0("mean_", groupnames[2]),  # 2
         "test_statistic",                # 3
         "p", "p.fdr",                    # 4, 5
         paste0("FC_", groupnames[2], "_divided_by_", groupnames[1]),  # 6
         paste0("FC_", groupnames[1], "_divided_by_", groupnames[2]),  # 7
         "CI_lower", "CI_upper",                                       # 8, 9
         paste0("n_", groupnames[1]), paste0("n_", groupnames[2]), # 10, 11
         "NA_reason_code")                                               # 12
  res <- rep(NA, length(cnames))
  names(res) <- cnames

  ### drop missing values
  tmp_na.omit <- droplevels(stats::na.omit(tmp))

  if (length(table(tmp_na.omit$group)) < 2) { # if less than 2 groups remain
    res[12] <- 2 ### reason: likely on/off protein with one group missing
    return(res)
  }              # ensure that both groups are still present
  if (!is.null(minObs) & any(table(tmp_na.omit$group) < minObs)) {  # ensure that each group has enough observations
    res[12]   <- 1  ### reason: not enough observations in at least one group
    return(res)
  }
  if (!is.null(minObsRatio) & any(table(tmp_na.omit$group)/table(tmp$group) < minObsRatio)) {  # ensure that each group has enough observations
    res[12]   <- 1  ### reason: not enough observations in at least one group
    return(res)
  }
  ttest <- try({stats::t.test(x = tmp_na.omit$intensity[tmp_na.omit$group == groupnames[2]],
                              y = tmp_na.omit$intensity[tmp_na.omit$group == groupnames[1]],
                              paired = FALSE, var.equal = varEqual)}, silent = TRUE)

  # it is still possible, that the ttest fails (e.g. if variance in one group is 0)
  if ("try-error" %in% class(ttest)) {
    print(ttest)
    #if (verbose) message("ttest failed for row ", row)
    res[12] <- 3
    return(res)
  }  ### reason: other, e.g. var = 0

  res[1:2] <- c(ttest$estimate[2], ttest$estimate[1])
  res[3] <- ttest$statistic
  res[4] <- ttest$p.value
  res[8:9] <- ttest$conf.int
  res[10] <- length(tmp_na.omit$intensity[tmp_na.omit$group == groupnames[1]])
  res[11] <- length(tmp_na.omit$intensity[tmp_na.omit$group == groupnames[2]])


  ### calculate fold changes
  #xFC <- ifelse(delogForFC, logBase^tmp$abundance, tmp$abundance)

  if (delogForFC) xFC <- logBase^tmp$intensity else xFC <- tmp$intensity

  res[5] <- mean(xFC[tmp$group == groupnames[2]], na.rm = TRUE) /
    mean(xFC[tmp$group == groupnames[1]], na.rm = TRUE)
  res[7] <- 1/res[6]

  return(res)
}




################################################################################
################################################################################
################################################################################


#' Paired t-test for a single row of a data set.
#'
#' @description
#' Runs an paired t-test on a single row of peptide or protein intensity data
#' and returns p-values and fold changes.
#'
#' @inheritParams .ttest_single_row
#' @param minPairs **integer(1)** \cr
#' The minimum number of complete sample pairs required. Default is 3.
#'
#' @return A named numeric vector with the following elements: mean difference
#' between groups, test statistic, p-value, placeholder for FDR-corrected
#' p-value, fold changes in both directions, confidence interval bounds, and
#' number of valid values per group.
#'
#' @seealso [ttest()], [.ttest_single_row()]
#'
#' @importFrom checkmate assertNumeric assertFactor assertFlag assertNumber
#' @importFrom checkmate assertIntegerish assertInt
#' @importFrom stats t.test
.ttest_single_row_paired <- function(x, group, sample, logBeforeTest = TRUE,
                              delogForFC = TRUE, minPairs = 3,
                             logBase = 2, row = NULL, verbose = TRUE) {

  checkmate::assertNumeric(x)
  checkmate::assertFactor(group, len = length(x), n.levels = 2)
  checkmate::assertFlag(logBeforeTest)
  checkmate::assertFlag(delogForFC)
  checkmate::assertNumber(logBase, lower = 1)
  checkmate::assertIntegerish(minPairs, null.ok = TRUE, lower = 2)
  checkmate::assertInt(row, null.ok = TRUE)
  checkmate::assertFlag(verbose)

  ## throw error if the two groups do not have the same length
  if (table(group)[1] != table(group)[2]) {
    stop("Groups don't have the same size, which is required for a paired t-test.")
  }

  ## log-transformation
  x <- unlist(unname(x))
  if (logBeforeTest) intensity <- log(x, base = logBase) else intensity <- x
  tmp <- data.frame(intensity = intensity, group = group, sample = sample)
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
  diffs <- tmp_group1$intensity - tmp_group2$intensity
  ind_complete_pairs <- which(!is.na(diffs))

  if (sum(!is.na(diffs)) < minPairs) {return(res)}  # ensure that enough complete pairs are present

  ttest <- try({stats::t.test(y = tmp_group1$intensity[ind_complete_pairs],  ## not possible to use tmp_na.omit as pairs may be shifted by omitting NAs
                       x = tmp_group2$intensity[ind_complete_pairs],
                       paired = TRUE)}, silent = TRUE)

  # it is still possible, that the ttest fails (e.g. if variance in one group is 0)
  if ("try-error" %in% class(ttest)) {message(paste0("ttest failed for row ", row));return(res)}

  res[1] <- ttest$estimate
  res[2] <- ttest$statistic
  res[3] <- ttest$p.value
  res[7:8] <- ttest$conf.int
  res[9] <- sum(!is.na(tmp_group1$intensity))
  res[10] <- sum(!is.na(tmp_group2$intensity))

  ### calculate fold changes
  if (delogForFC) xFC <- logBase^tmp$intensity else xFC <- tmp$intensity

  res[5] <- mean(xFC[tmp$group == groupnames[2]][ind_complete_pairs], na.rm = TRUE) /
    mean(xFC[tmp$group == groupnames[1]][ind_complete_pairs], na.rm = TRUE)
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
#' @param SE  **SummarizedExperiment object** \cr
#' Data in a SummarizedExperiment object, e.g. output SE from [prepareData]
#' function.
#' @param assay  **character(1)** \cr
#' The name of the assay in SE containing the protein intensities that should
#' be used for the PCA plot. Default is "intensity_norm", which uses the
#' normalized protein intensities if the SE object was generated by
#' [prepareData].
#' @param groupColumn **character(1)** \cr
#' The name of the column in the colData of SE containing the groups that will
#' be compared by the t-test. Please make sure that exactly two groups are
#' present.
#' @param sampleColumn **character(1)** \cr
#' The name of the column in the colData of SE containing the sample IDs. Only
#' necessary if `paired = TRUE`, then it will be used to match sample pairs.
#' @param paired **logical(1)** \cr
#' If TRUE, a paired t-test is performed, otherwise unpaired. Default is FALSE.
#' @inheritParams .ttest_single_row varEqual logBeforeTest delogForFC logbase
#' @param minObs **integer(1)** \cr
#' Minimum number of observations per group. For a paired t-test this is the
#' minimum number of complete pairs. Default is 3.
#' @param minObsRatio **numeric(1)** \cr
#' Minimum proportion of valid values required per group (e.g. 0.8 = 80%).
#' Currently not usable if `paired = TRUE`.
#'
#' @return A data frame containing the t-test results an used protein intensity
#' values.
#' @export
#'
#' @seealso For more details see [.ttest_single_row()] for the unpaired and
#'  [.ttest_single_row_paired()] for the paired t-test.
#'
#' @importFrom checkmate assertTRUE assertSubset assertFlag
#' @importFrom methods is
#' @importFrom pbapply pbapply
#' @importFrom stats p.adjust
#' @importFrom SummarizedExperiment assays colData rowData
#'
#'@examples
#'file_proteins <- system.file("extdata", "proteins_HCC.csv", package = "ProtStatsWF")
#'file_clinical  <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")
#'
#'D_hcc <- prepareData(dataPath = file_proteins, intensityColumns = 6:43,
#'                     proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
#'                     sampleNameColumn = "Sample", verbose = FALSE)
#'
#'ttest(SE = D_hcc$SE, assay = "intensity_norm",
#'                   groupColumn = "Group", sampleColumn = "PatientID",,
#'                   logBeforeTest = FALSE, delogForFC = TRUE, logBase = 2,
#'                   minObs = 3, paired = TRUE)
ttest <- function(SE, assay, groupColumn, sampleColumn = NULL, paired = FALSE,
                  varEqual = FALSE, logBeforeTest = TRUE, delogForFC = TRUE,
                  logBase = 2, minObs = 3, minObsRatio = NULL, verbose = TRUE) {

  checkmate::assertTRUE(methods::is(SE, "SummarizedExperiment"))
  checkmate::assertSubset(assay, names(SummarizedExperiment::assays(SE)))
  checkmate::assertSubset(groupColumn, colnames(SummarizedExperiment::colData(SE)))
  checkmate::assertSubset(sampleColumn, colnames(SummarizedExperiment::colData(SE)))
  checkmate::assertFlag(paired)

  D <- as.data.frame(SummarizedExperiment::assays(SE)[[assay]])
  id <- as.data.frame(SummarizedExperiment::rowData(SE))
  group <- SummarizedExperiment::colData(D_hcc$SE)[,groupColumn]
  if (!is.null(sampleColumn)) sample <- SummarizedExperiment::colData(D_hcc$SE)[,sampleColumn]

  if (!paired) {
    RES <- pbapply::pbapply(D, 1, .ttest_single_row, group = group,
                            logBeforeTest = logBeforeTest,
                 delogForFC = delogForFC, minObs = minObs,
                 logBase = logBase, varEqual = varEqual,
                 minObsRatio = minObsRatio)
  }

  if (paired) {
    RES <- pbapply::pbapply(D, 1, .ttest_single_row_paired, group = group,
                            logBeforeTest = logBeforeTest,
                 delogForFC = delogForFC, minPairs = minObs, sample = sample,
                 logBase = logBase)
  }

  RES <- t(RES)
  RES <- as.data.frame(RES)

  RES$p.fdr <- stats::p.adjust(RES$p, method = "fdr")

  if (logBeforeTest) {
    D_log <- log(D, base = logBase)
    colnames(D_log) <- paste0(colnames(D), "_log")
  } else {
    if (delogForFC) {
      D_delog <- logBase^D
      colnames(D_delog) <- paste0(colnames(D), "_delog")
    }
  }


  ## sort columns in D by group and sample
  if (!is.null(sampleColumn)) {
    D <- D[,order(group, sample)]
  } else {
    D <- D[,order(group)]
  }
  return(cbind(id, RES, D))
}




