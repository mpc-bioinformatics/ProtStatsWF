

#' ANOVA to compare three or more experimental groups
#'
#' @param D                    \strong{data.frame} \cr
#'                             The data set containing only protein intensities of the sample.
#' @param id                   \strong{data.frame} \cr
#'                             The corresponding ID columns for the parameter D.
#' @param group                \strong{character factor} \cr
#'                             The groups of the data.
#' @param sample               \strong{character factor} \cr
#'                             The sample membership in the data.
#' @param paired               \strong{logical} \cr
#'                             If \code{TRUE}, a paired test will be done, otherwise an unpaired test.
#' @param var.equal            \strong{logical} \cr
#'                             If \code{TRUE}, the variances are assumed to be equal.
#' @param log_before_test      \strong{logical} \cr
#'                             If \code{TRUE}, the data will be log-transformed.
#' @param delog_for_FC         \strong{logical} \cr
#'                             If \code{TRUE}, the fold change will be calculated without the log-transformation.
#' @param log_base             \strong{integer} \cr
#'                             The base for the log-transformation, if \code{log_before_test = TRUE}.
#' @param min_obs_per_group    \strong{integer} \cr
#'                             The minimum number of observations per group.
#' @param min_perc_per_group   \strong{integer} \cr
#'                             The minimum ratio of observations per group as an alternative to min_obs_per_group.
#'
#' @return A data.frame with p-values and fold changes
#' @export
#'
#' @examples
#'

ANOVA <- function(D,
                  id = NULL,
                  group,
                  sample = NULL,
                  paired = FALSE,
                  var.equal = TRUE,
                  log_before_test = TRUE,
                  delog_for_FC = TRUE,
                  log_base = 2,
                  min_obs_per_group = 3,
                  min_perc_per_group = NULL) {

  if (!is.null(min_obs_per_group) & !is.null(min_perc_per_group)) stop("Please specify only one of min_obs_per_group or min_perc_per_group.")


  if (paired) {
    print("Repeated Measures ANOVA")
    ### repeated measures ANOVA
    RES <- pbapply::pbapply(D, 1, ANOVA_repeatedMeasurements_single_row, group = group, sample = sample,
                   log_before_test = log_before_test,
                   delog_for_FC = delog_for_FC, min_obs_per_group = min_obs_per_group,
                   min_perc_per_group = min_perc_per_group,
                   log_base = log_base)
  } else {
    if (var.equal) {
      #### Standard ANOVA (equeal variances
      print("Standard ANOVA")
      RES <- pbapply::pbapply(D, 1, ANOVA_standard_single_row, group = group, log_before_test = log_before_test,
                   delog_for_FC = delog_for_FC, min_obs_per_group = min_obs_per_group,
                   min_perc_per_group = min_perc_per_group,
                   log_base = log_base)
    } else {
      ### Welch ANOVA (unequal variances)
      print("Welch ANOVA")
      i <<- 0
      RES <- pbapply::pbapply(D, 1, ANOVA_Welch_single_row, group = group, log_before_test = log_before_test,
                     delog_for_FC = delog_for_FC, min_obs_per_group = min_obs_per_group,
                     min_perc_per_group = min_perc_per_group,
                     log_base = log_base)
    }
  }

  RES <- t(RES)
  RES <- as.data.frame(RES)

  # adjust the p-values of the ANOVA
  RES$p.anova.fdr <- stats::p.adjust(RES$p.anova, method = "fdr")

  if (!is.null(id)) {
    D <- cbind(id, D)
  }
  #openxlsx::write.xlsx(cbind(D, RES), filename, keepNA = TRUE, overwrite = TRUE)

  return(cbind(D, RES))
}

################################################################################
################################################################################
################################################################################



#' Standard ANOVA (equal variances)
#'
#' @param x Numeric vector of protein intensities.
#' @param group Factor vector indicating the group membership.
#' @param min_obs_per_group Integer indicating the minimum number of observations per group.
#' @param log_before_test Logical indicating whether the data should be log-transformed before the test.
#' @param log_base Numeric indicating the base of the logarithm.
#' @param delog_for_FC Logical indicating whether the fold change should be calculated on the original scale.
#' @param min_perc_per_group Numeric indicating the minimum ratio of observations per group.
#'
#' @return Vector with p-values and fold changes.
#' @export
#'
#' @examples # TODO
ANOVA_standard_single_row <- function(x,
                                      group,
                                      min_obs_per_group = 3,
                                      log_before_test = TRUE,
                                      log_base = 2,
                                      delog_for_FC = TRUE,
                                      min_perc_per_group = NULL) {

  x <- as.numeric(unname(x))

  nr_groups <- length(levels(group))
  groups <- factor(levels(group), levels = levels(group))

  comparisons <- gtools::combinations(n = nr_groups, r = 2)
  comparisons_char <- matrix(data = NA_character_, nrow = nrow(comparisons), ncol = ncol(comparisons))
  for (i in 1:nr_groups) {
    comparisons_char[comparisons == i] <- as.character(groups[i])
  }
  comparisons <- comparisons_char
  nr_comparisons <- nrow(comparisons)

  if (log_before_test) x <- log(x, base = log_base)

  D_tmp <- data.frame(intensity = x, group)
  D_tmp_naomit <- stats::na.omit(D_tmp)

  ## remove groups with less than min_obs_per_group values
  tab_lev <- table(D_tmp_naomit$group) # table of groups without NAs
  tab_lev_perc <- tab_lev/table(group)
  if (!is.null(min_obs_per_group) & any(tab_lev < min_obs_per_group)) {
    group_remove <- names(tab_lev)[which(tab_lev < min_obs_per_group)]
    D_tmp_naomit <- D_tmp_naomit[!D_tmp_naomit$group %in% group_remove, ]
    D_tmp_naomit <- droplevels(D_tmp_naomit)
  }
  if (!is.null(min_perc_per_group) & any(tab_lev_perc < min_perc_per_group)) {
    group_remove <- names(tab_lev_perc)[which(tab_lev_perc < min_perc_per_group)]
    D_tmp_naomit <- D_tmp_naomit[!D_tmp_naomit$group %in% group_remove, ]
    D_tmp_naomit <- droplevels(D_tmp_naomit)
  }


  ### TODO: die Benennung der Spalten in die aeussere Funktion verlagern
  namesfc1 <- paste0("FC_", comparisons[,1], "_divided_by_", comparisons[,2])
  namesfc2 <- paste0("FC_", comparisons[,2], "_divided_by_", comparisons[,1])
  namesfc <- character(2*nr_comparisons)
  namesfc[seq(1,2*nr_comparisons, by = 2)] <- namesfc1
  namesfc[seq(2,2*nr_comparisons, by = 2)] <- namesfc2
  cnames <- c("p.anova", "p.anova.fdr",
              paste0("p.posthoc.", comparisons[,1], "_vs_", comparisons[,2]),
              namesfc)


  if (length(levels(D_tmp_naomit$group)) <= 1) {
    res <- rep(NA, 3 * nr_comparisons + 2)
    names(res) <- cnames
    return(res)
  } # ANOVA + adj + posthocs + fold changes

  aov.obj  <- stats::aov(intensity ~ group, data = D_tmp_naomit)
  sum.aov  <- summary(aov.obj)
  p.anova  <- sum.aov[[1]]['group','Pr(>F)']
  p.posthoc_tmp  <- stats::TukeyHSD(aov.obj )$group[,'p adj']
  p.anova.fdr <- NA

  p.posthoc <- rep(NA, nr_comparisons)

  if (length(levels(D_tmp_naomit$group)) == 2) {
    for (j in 1:nr_comparisons) {
      if (all(levels(D_tmp_naomit$group) %in% comparisons[j,])) {
        p.posthoc[j] <- p.anova
      } else {
        p.posthoc[j] <- NA
      }
    }
  } else {

    for (j in 1:nr_comparisons) {
      name_tmp <- paste0(comparisons[j,2], "-", comparisons[j,1])
      if(any(names(p.posthoc_tmp) == name_tmp)) {
        p.posthoc[j] <- p.posthoc_tmp[name_tmp]
      } else {
        p.posthoc[j] <- NA
      }
    }
  }

  if (delog_for_FC) {
    x2 <- log_base^x
  } else {
    x2 <- x
  }

  fcs <- NULL
  name.fcs <- NULL
  for (j in 1:(nr_groups-1)) {
    for (k in (j+1):nr_groups) {
      fc1 <- mean(x2[group == groups[j]], na.rm = TRUE) / mean(x2[group == groups[k]], na.rm = TRUE)
      fc2 <- 1/fc1
      fcs <- cbind(fcs, fc1, fc2)
    }
  }

  res <- c(p.anova = p.anova, p.anova.fdr = p.anova.fdr, p.posthoc, fcs)
  names(res) <- cnames

  return(res)
}




################################################################################
################################################################################


#' Repeated measures ANOVA (paired samples)
#'
#' @param x Numeric vector of protein intensities.
#' @param group Factor vector indicating the group membership.
#' @param sample vector indicating the sample.
#' @param min_obs_per_group Integer indicating the minimum number of observations per group.
#' @param log_before_test Logical indicating whether the data should be log-transformed before the test.
#' @param log_base Numeric indicating the base of the logarithm.
#' @param delog_for_FC Logical indicating whether the fold change should be calculated on the original scale.
#' @param min_perc_per_group Numeric indicating the minimum ratio of observations per group.
#'
#' @return Vector with p-values and fold changes.
#' @export
#'
#' @examples # TODO
ANOVA_repeatedMeasurements_single_row <- function(x,
                                                  group,
                                                  sample,
                                                  min_obs_per_group = 3,
                                                  log_before_test = TRUE,
                                                  log_base = 2,
                                                  delog_for_FC = TRUE,
                                                  min_perc_per_group = NULL) {
  x <- as.numeric(unname(x))

  nr_groups <- length(levels(group))
  groups <- factor(levels(group), levels = levels(group))

  h <- 2 + nr_groups + 2*choose(nr_groups, 2)

  comparisons <- gtools::combinations(n = nr_groups, r = 2)
  comparisons_char <- matrix(data = NA_character_, nrow = nrow(comparisons), ncol = ncol(comparisons))
  for (i in 1:nr_groups) {
    comparisons_char[comparisons == i] <- as.character(groups[i])
  }
  comparisons <- comparisons_char
  nr_comparisons <- nrow(comparisons)

  if (log_before_test) x <- log(x, base = log_base)

  D_tmp <- data.frame(intensity = x, group = group, sample = factor(sample))
  D_tmp_naomit <- stats::na.omit(D_tmp)


  ## remove groups with less than min_obs_per_group values
  tab_lev <- table(D_tmp_naomit$group) # table of groups without missing values
  tab_lev_perc <- tab_lev/table(group)
  if (!is.null(min_obs_per_group) & any(tab_lev < min_obs_per_group)) {
    group_remove <- names(tab_lev)[which(tab_lev < min_obs_per_group)]
    D_tmp_naomit <- D_tmp_naomit[!D_tmp_naomit$group %in% group_remove, ]
    D_tmp_naomit <- droplevels(D_tmp_naomit)
  }
  if (!is.null(min_perc_per_group) & any(tab_lev_perc < min_perc_per_group)) {
    group_remove <- names(tab_lev_perc)[which(tab_lev_perc < min_perc_per_group)]
    D_tmp_naomit <- D_tmp_naomit[!D_tmp_naomit$group %in% group_remove, ]
    D_tmp_naomit <- droplevels(D_tmp_naomit)
  }

  ### TODO: die Benennung der Spalten in die aeussere Funktion verlagern
  namesfc1 <- paste0("FC_", comparisons[,1], "_divided_by_", comparisons[,2])
  namesfc2 <- paste0("FC_", comparisons[,2], "_divided_by_", comparisons[,1])
  namesfc <- character(2*nr_comparisons)
  namesfc[seq(1,2*nr_comparisons, by = 2)] <- namesfc1
  namesfc[seq(2,2*nr_comparisons, by = 2)] <- namesfc2
  cnames <- c("p.anova", "p.anova.fdr",
              paste0("p.posthoc.", comparisons[,1], "_vs_", comparisons[,2]),
              namesfc)

  if (length(levels(D_tmp_naomit$group)) <= 1) {  ## if <2 groups are left
    res <- rep(NA, 3 * nr_comparisons + 2) # ANOVA + adj + posthocs + fold changes
    names(res) <- cnames
    return(res)
  }


  Lme.mod <- try({nlme::lme(intensity ~ group, random = ~1 | sample, data = D_tmp_naomit, na.action = stats::na.omit)})
  if ("try-error" %in% class(Lme.mod)) {
    res <- rep(NA, 3 * nr_comparisons + 2)
    names(res) <- cnames
    return(res)
    }  # Modellberechnung fehlgeschlagen
  lme.aov <- stats::anova(Lme.mod)  # Wald test
  p.anova <- lme.aov$`p-value`[2]  # p-Wert der ANOVA fuer den Faktor group
  if (is.nan(p.anova)) {
    z <- rep(NA, h)
    names(z) <- cnames
    return(z)
  }
  p.anova.fdr <- NA

  posthoc <- summary(multcomp::glht(Lme.mod, linfct = multcomp::mcp(group = "Tukey")))  # Tukey posthoc Tests
  p.posthoc_tmp <- as.vector(posthoc$test$pvalues)
  names_p.posthoc_tmp <- unname(dimnames(posthoc$linfct)[[1]])


  p.posthoc <- rep(NA, nr_comparisons)
  for (j in 1:nr_comparisons) {
    name_tmp <- paste0(comparisons[j,2], " - ", comparisons[j,1])
    if (any(names_p.posthoc_tmp == name_tmp)) {
      p.posthoc[j] <- p.posthoc_tmp[names_p.posthoc_tmp == name_tmp]
    } else {
      p.posthoc[j] <- NA
    }
  }


  if (delog_for_FC) {
    x2 <- log_base^x
  } else {
    x2 <- x
  }

  fcs <- NULL
  name.fcs <- NULL
  for (j in 1:(nr_groups - 1)) {
    for (k in (j + 1):nr_groups) {
      fc1 <- mean(x2[group == groups[j]], na.rm = TRUE) / mean(x2[group == groups[k]], na.rm = TRUE)
      fc2 <- 1/fc1
      fcs <- cbind(fcs, fc1, fc2)
    }
  }
  res <- c(p.anova = p.anova, p.anova.fdr = p.anova.fdr, p.posthoc, fcs)
  names(res) <- cnames

  return(res)
}




################################################################################
################################################################################
################################################################################



#' Welch ANOVA (unequal variances)
#'
#' @param x Numeric vector of protein intensities.
#' @param group Factor vector indicating the group membership.
#' @param min_obs_per_group Integer indicating the minimum number of observations per group.
#' @param log_before_test Logical indicating whether the data should be log-transformed before the test.
#' @param log_base Numeric indicating the base of the logarithm.
#' @param delog_for_FC Logical indicating whether the fold change should be calculated on the original scale.
#' @param min_perc_per_group Numeric indicating the minimum ratio of observations per group.
#'
#' @return Vector with p-values and fold changes.
#' @export
#'
#' @examples # TODO
ANOVA_Welch_single_row <- function(x,
                                   group,
                                   min_obs_per_group,
                                   log_before_test = TRUE,
                                   log_base = 2,
                                   delog_for_FC = TRUE,
                                   min_perc_per_group = NULL) {

  x <- as.numeric(unname(x))

  nr_groups <- length(levels(group))
  groups <- factor(levels(group), levels = levels(group))

  comparisons <- gtools::combinations(n = nr_groups, r = 2)
  comparisons_char <- matrix(data = NA_character_, nrow = nrow(comparisons), ncol = ncol(comparisons))
  for (i in 1:nr_groups) {
    comparisons_char[comparisons == i] <- as.character(groups[i])
  }
  comparisons <- comparisons_char
  nr_comparisons <- nrow(comparisons)

  if (log_before_test) x <- log(x, base = log_base)

  D_tmp <- data.frame(intensity = x, group = group)
  D_tmp_naomit <- stats::na.omit(D_tmp)

  ## remove groups with less than min_obs_per_group values
  tab_lev <- table(D_tmp_naomit$group) # table of groups without NAs
  tab_lev_perc <- tab_lev/table(group)
  if (!is.null(min_obs_per_group) & any(tab_lev < min_obs_per_group)) {
    group_remove <- names(tab_lev)[which(tab_lev < min_obs_per_group)]
    D_tmp_naomit <- D_tmp_naomit[!D_tmp_naomit$group %in% group_remove, ]
    D_tmp_naomit <- droplevels(D_tmp_naomit)
  }
  if (!is.null(min_perc_per_group) & any(tab_lev_perc < min_perc_per_group)) {
    group_remove <- names(tab_lev_perc)[which(tab_lev_perc < min_perc_per_group)]
    D_tmp_naomit <- D_tmp_naomit[!D_tmp_naomit$group %in% group_remove, ]
    D_tmp_naomit <- droplevels(D_tmp_naomit)
  }

  ### TODO: die Benennung der Spalten in die aeussere Funktion verlagern
  namesfc1 <- paste0("FC_", comparisons[,1], "_divided_by_", comparisons[,2])
  namesfc2 <- paste0("FC_", comparisons[,2], "_divided_by_", comparisons[,1])
  namesfc <- character(2*nr_comparisons)
  namesfc[seq(1,2*nr_comparisons, by = 2)] <- namesfc1
  namesfc[seq(2,2*nr_comparisons, by = 2)] <- namesfc2
  cnames <- c("p.anova", "p.anova.fdr",
              paste0("p.posthoc.", comparisons[,1], "_vs_", comparisons[,2]),
              namesfc)

  if (length(levels(D_tmp_naomit$group)) <= 1) {  # if <= 1 group is left
    res <- rep(NA, 3 * nr_comparisons + 2) # ANOVA + adj + posthocs + fold changes
    names(res) <- cnames
    return(res)
  }

  model <- stats::lm(intensity ~ group, data = D_tmp_naomit)
  ANOVA <- try({suppressMessages(car::Anova(model, white.adjust = TRUE))})

  if ("try-error" %in% class(ANOVA)) {
    res <- rep(NA, 3 * nr_comparisons + 2)
    names(res) <- cnames
    return(res)
  }

  p.anova <- ANOVA$`Pr(>F)`[1]
  p.anova.fdr <- NA

  ### single Welch-Tests as post-hoc tests
  p.posthoc <- rep(NA, nrow(comparisons))
  for (j in 1:nrow(comparisons)) {
    D_tmp1 <- D_tmp_naomit[D_tmp_naomit$group %in% comparisons[j, ],]
    D_tmp1$group <- droplevels(D_tmp1$group)
    if (length(levels(D_tmp1$group)) != 2) {  # if there are not enough values, skip
      next
    } else {
      p.posthoc[j] <- stats::t.test(intensity ~ group, data = D_tmp1)$p.value
    }
  }
  p.posthoc <- stats::p.adjust(p.posthoc, method = "holm")


  if (delog_for_FC) {
    x2 <- log_base^x
  } else {
    x2 <- x
  }

  fcs <- NULL
  name.fcs <- NULL
  for (j in 1:(nr_groups-1)) {
    for (k in (j+1):nr_groups) {
      fc1 <- mean(x2[group == groups[j]], na.rm = TRUE) / mean(x2[group == groups[k]], na.rm = TRUE)
      fc2 <- 1/fc1
      fcs <- cbind(fcs, fc1, fc2)
    }
  }

  res <- c(p.anova = p.anova, p.anova.fdr = p.anova.fdr, p.posthoc, fcs)
  names(res) <- cnames

  return(res)
}



