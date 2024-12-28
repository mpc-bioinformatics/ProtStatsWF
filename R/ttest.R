


#' Unpaired t-test for a single row of a data set.
#'
#' @param x                         \strong{numeric vector} \cr
#'                                  The abundances of the data.
#' @param group                     \strong{character factor} \cr
#'                                  The group membership of the data.
#' @param log_before_test           \strong{logical} \cr
#'                                  If \code{TRUE}, the data will be log-transformed before the test.
#' @param delog_for_FC              \strong{logical} \cr
#'                                  If \code{TRUE}, the fold change will be calculated on the original scale.
#' @param min_obs_per_group         \strong{integer} \cr
#'                                  The minimum number of observations per group.
#' @param min_obs_per_group_ratio   \strong{numeric} \cr
#'                                  The minimum number of observations per group as a ratio (e.g, 0.8 = 80% valid values in each group needed).
#' @param log_base                  \strong{numeric} \cr
#'                                  The base of the logarithm for the log-transformation.
#' @param row                       \strong{integer} \cr
#'                                  The row number of the data for the function call.
#' @param var.equal                 \strong{logical} \cr
#'                                  If \code{TRUE}, the variances of the groups are expected to be equal.
#'
#' @return A vector with the following components: mean group 1, mean group 2, test statistics, p-value, free space fpr corrected p-value, fold changes (both directions), lower and upper limit of confidence interval, number of valid values per group.
#' @export
#' 
#' @seealso [ttest()], [ttest_single_row_paired()]
#'
#' @examples 
#' 

ttest_single_row <- function(x, group, log_before_test = TRUE, delog_for_FC = TRUE,
                             min_obs_per_group = NULL, min_obs_per_group_ratio = NULL,
                             log_base = 2, row = NULL, var.equal = FALSE) {

  if (!is.null(min_obs_per_group) & !is.null(min_obs_per_group_ratio)) {
    stop("Both min_obs_per_group and min_obs_per_group_ratio are given, please define only one of them!")
  }



  x <- unlist(unname(x))
  if (log_before_test) {
    abundance = log(x, base = log_base)
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

  if(!is.null(min_obs_per_group) & any(table(tmp_na.omit$group) < min_obs_per_group)) {  # ensure that each group has enough observations
    res[12] <- 1  ### reason: not enough observations in at least one group
    return(res)
  }
  if(!is.null(min_obs_per_group_ratio) & any(table(tmp_na.omit$group)/table(tmp$group) < min_obs_per_group_ratio)) {  # ensure that each group has enough observations
    res[12] <- 1  ### reason: not enough observations in at least one group
    return(res)
  }

  ttest <- try({stats::t.test(x = tmp_na.omit$abundance[tmp_na.omit$group == groupnames[2]],
                       y = tmp_na.omit$abundance[tmp_na.omit$group == groupnames[1]],
                       paired = FALSE, var.equal = var.equal)}, silent = TRUE)

  # it is still possible, that the ttest fails (e.g. if variance in one group is 0)
  if ("try-error" %in% class(ttest)) {warning(paste0("ttest failed for row ", row));res[12] <- 3;return(res)}  ### reason: other, e.g. var = 0
  ### TODO: error message mit ausgeben

  res[1:2] <- ttest$estimate
  res[3] <- ttest$statistic
  res[4] <- ttest$p.value
  res[5] <- NA # free space for corrected p-value
  res[8:9] <- ttest$conf.int
  res[10] <- length(tmp_na.omit$abundance[tmp_na.omit$group == groupnames[1]])
  res[11] <- length(tmp_na.omit$abundance[tmp_na.omit$group == groupnames[2]])


  ### calculate fold changes
  if (delog_for_FC) {
    x2 <- log_base^tmp$abundance
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

### TODO: add possibility to calculate FC pair-wise (diff of the values on log-scale, then mean, then delog)
### TODO: add possibility to set min_nr_pairs as a ratio!
### TODO: add error codes

#' Paired T-test for a single row of a data set.
#'
#' @param x                         \strong{numeric vector} \cr
#'                                  The abundances of the data.
#' @param group                     \strong{character factor} \cr
#'                                  The group membership of the data.
#' @param sample                    \strong{character factor} \cr
#'                                  The sample membership of the data.
#' @param log_before_test           \strong{logical} \cr
#'                                  If \code{TRUE}, the data will be log-transformed before the test.
#' @param delog_for_FC              \strong{logical} \cr
#'                                  If \code{TRUE}, the fold change will be calculated on the original scale.
#' @param min_nr_pairs              \strong{integer} \cr
#'                                  The minimum number of complete sample pairs.
#' @param log_base                  \strong{numeric} \cr
#'                                  The base of the logarithm for the log-transformation.
#' @param row                       \strong{integer} \cr
#'                                  The row number of the data for the function call.
#'
#' @return A vector with the following components: mean difference between groups, test statistics, p-value, free space fpr corrected p-value, fold changes (both directions), lower and upper limit of confidence interval, number of valid values per group.
#' @export
#' 
#' @seealso [ttest()], [ttest_single_row()]
#'
#' @examples 
#' 

ttest_single_row_paired <- function(x, group, sample, log_before_test = TRUE, delog_for_FC = TRUE,
                             min_nr_pairs = NULL,
                             log_base = 2, row = NULL) {

  ## throw error if the two groups do not have the same length
  if (table(group)[1] != table(group)[2]) {
    stop("Groups don't have the same size, which is required for a paired t-test.")
  }

  ## log-transformation
  x <- unname(x)
  if (log_before_test) {
    abundance = log(x, base = log_base)
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

  if (sum(!is.na(diffs)) < min_nr_pairs) {return(res)}  # ensure that enough complete pairs are present

  ttest <- try({stats::t.test(y = tmp_group1$abundance[ind_complete_pairs],  ## not possible to use tmp_na.omit as pairs may be shifted by omitting NAs
                       x = tmp_group2$abundance[ind_complete_pairs],
                       paired = TRUE)}, silent = TRUE)

  # it is still possible, that the ttest fails (e.g. if variance in one group is 0)
  if ("try-error" %in% class(ttest)) {warning(paste0("ttest failed for row ", row));return(res)}
  ### TODO: error message mit ausgeben

  res[1] <- ttest$estimate
  res[2] <- ttest$statistic
  res[3] <- ttest$p.value
  res[4] <- NA # free space for corrected p-value
  res[7:8] <- ttest$conf.int
  res[9] <- sum(!is.na(tmp_group1$abundance))
  res[10] <- sum(!is.na(tmp_group2$abundance))

  ### calculate fold changes
  if (delog_for_FC) {
    x2 <- log_base^tmp$abundance
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


#' Function to compute t-test (paired or unpaired)
#' @param D                         \strong{data.frame} \cr
#'                                  The data set containing only protein intensities, already filtered for interesting candidates.
#' @param id                        \strong{data.frame} \cr
#'                                  The corresponding ID columns for the parameter D e.g. containing further columns like protein or gene names.
#' @param group                     \strong{character factor} \cr
#'                                  The group membership of the data.
#' @param sample                    \strong{character factor} \cr
#'                                  The sample membership of the data.
#' @param paired                    \strong{logical} \cr
#'                                  If \code{TRUE}, the test will be paired, otherwise it will be unpaired.
#' @param var.equal                 \strong{logical} \cr
#'                                  If \code{TRUE}, the variances of the groups are expected to be equal.
#' @param log_before_test           \strong{logical} \cr
#'                                  If \code{TRUE}, the data will be log-transformed before the test.
#' @param delog_for_FC              \strong{logical} \cr
#'                                  If \code{TRUE}, the fold change will be calculated on the original scale.
#' @param log_base                  \strong{numeric} \cr
#'                                  The base of the logarithm for the log-transformation.
#' @param min_obs_per_group         \strong{integer} \cr
#'                                  The minimum number of observations per group.
#' @param filename                  \strong{characte} \cr
#'                                  The name of the output file.
#' @param min_obs_per_group_ratio   \strong{numeric} \cr
#'                                  The minimum number of observations per group as a ratio (e.g, 0.8 = 80% valid values in each group needed).
#'
#' @return A data frame containing the results of the t-test.
#' @export
#' 
#' @seealso [ttest_single_row()], [ttest_single_row_paired()]
#'
#' @examples 
#' 

ttest <- function(D, id = NULL, group, sample = NULL, paired = FALSE, var.equal = FALSE,
                  log_before_test = TRUE, delog_for_FC = TRUE, log_base = 2,
                  min_obs_per_group = 3,
                  filename = "results_ttest.xlsx", min_obs_per_group_ratio = NULL) {

  if (!paired) {
    RES <- pbapply::pbapply(D, 1,ttest_single_row, group = group, log_before_test = log_before_test,
                 delog_for_FC = delog_for_FC, min_obs_per_group = min_obs_per_group,
                 log_base = log_base, var.equal = var.equal, min_obs_per_group_ratio = min_obs_per_group_ratio)
  }

  if (paired) {
    RES <- pbapply::pbapply(D, 1, ttest_single_row_paired, group = group, log_before_test = log_before_test,
                 delog_for_FC = delog_for_FC, min_nr_pairs = min_obs_per_group, sample = sample,
                 log_base = log_base)
  }

  RES <- t(RES)
  RES <- as.data.frame(RES)

  RES$p.fdr <- stats::p.adjust(RES$p, method = "fdr")

  if(log_before_test){
    D_log <- log(D, base = log_base)
    colnames(D_log) <- paste0(colnames(D), "_log")
  } else {
    if (delog_for_FC) {
      D_delog <- log_base^D
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
    if(log_before_test) {
      D <- cbind(id, D, D_log)
    } else {
      if (delog_for_FC) {
        D <- cbind(id, D, D_delog)
      } else {
        D <- cbind(id, D)
      }
    }
  }

  openxlsx::write.xlsx(cbind(D, RES), filename, keepNA = TRUE, overwrite = TRUE)

  return(cbind(D, RES))
}




