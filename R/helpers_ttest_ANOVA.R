#' Calculation of significance categories for a volcano plot for the t-test.
#'
#' This function groups all proteins into the following three catrgories based
#' on the p-values (with and without FDR correction) and fold changes:
#' 1) not signifcant (p > thresP or fc < thresFC)
#' 2) significant (p < thresP and fc > thresFC, but pAdj > thresP)
#' 3) significant after FDR-correction (pAdj < thresP and fc > thresFC)
#'
#' @param p          \strong{numeric vector} \cr
#'                   The p-values before FDR-correction.
#' @param pAdj       \strong{numeric vector} \cr
#'                   The p-values after FDR-correction.
#' @param fc         \strong{numeric vector} \cr
#'                   The values of the fold changes.
#' @param thresFC    \strong{numeric} \cr
#'                   The threshold for the fold changes.
#' @param thresP     \strong{numeric} \cr
#'                   The threshold for the p-values.
#'
#' @return A factor with the three significance categories.
#' @export
#'
#' @examples
#'

calculate_significance_categories_ttest <- function(p, pAdj, fc, thresFC = 2, thresP = 0.05) {

  significance <- dplyr::case_when(
    pAdj <= thresP & p <= thresP & (fc >= thresFC | fc <= 1/thresFC) & !is.na(p) ~ "significant after FDR correction",
    pAdj > thresP & p <= thresP & (fc >= thresFC | fc <= 1/thresFC) & !is.na(p) ~ "significant",
    (p > thresP | (fc < thresFC & fc > 1/thresFC)) & !is.na(p) ~ "not significant",
    is.na(p) ~ NA_character_
  )

  significance <- factor(significance, levels = c("not significant", "significant", "significant after FDR correction"))

  return(significance)
}



#' Calculate significance categories for ANOVA.
#'
#' @param p_posthoc     \strong{numeric vector} \cr
#'                      The posthoc p-values .
#' @param p_anova_adj   \strong{numeric vector} \cr
#'                      The p-values after FDR-correction.
#' @param p_anova       \strong{numeric vector} \cr
#'                      The p-values before FDR-correction
#' @param fc            \strong{numeric vector} \cr
#'                      The values of the fold changes.
#' @param thres_fc      \strong{numeric} \cr
#'                      The threshold for the fold changes.
#' @param thres_p       \strong{numeric} \cr
#'                      The threshold for the p-values.
#'
#' @return A factor containing the significances.
#' @export
#'
#' @examples
#'

calculate_significance_categories_ANOVA <- function(p_posthoc, p_anova_adj, p_anova, fc, thres_fc=2, thres_p=0.05) {

  significance <- dplyr::case_when(
    p_anova_adj <= thres_p & p_posthoc <= thres_p & (fc >= thres_fc | fc <= 1/thres_fc) & !is.na(p_posthoc) & !is.na(p_anova) ~ "significant after FDR correction", # ANOVA significant after FDR, posthoc also significant, fulfills FC threshold
    p_anova_adj > thres_p & p_anova <= thres_p & p_posthoc <= thres_p & (fc >= thres_fc | fc <= 1/thres_fc) & !is.na(p_posthoc) & !is.na(p_anova) ~ "significant", # ANOVA significant before FDR, posthoc also significant, fulfills FC threshold
    (p_anova > thres_p | p_posthoc > thres_p | (fc < thres_fc & fc > 1/thres_fc)) & !is.na(p_posthoc) & !is.na(p_anova) ~ "not significant", # ANOVA not significant or posthoc not significant or FC does not fulfill threshold
    is.na(p_posthoc) | is.na(p_anova) ~ NA_character_
  )

  significance <- factor(significance, levels = c("not significant", "significant", "significant after FDR correction"))

  return(significance)
}
