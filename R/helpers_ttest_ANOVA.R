#' Calculation of significance categories for a volcano plot for the t-test.
#'
#' This function groups all proteins into the following three categories based
#' on the p-values (with and without FDR correction) and fold changes:
#' 1) not significant (p > thresP or fc < thresFC)
#' 2) significant (p < thresP and fc > thresFC, but pAdj > thresP)
#' 3) significant after FDR-correction (pAdj < thresP and fc > thresFC)
#'
#' If a p-value is NA, the corresponding significance category will also be NA.
#'
#' @param p **numeric**\cr
#' Vector of p-values before FDR-correction.
#' @param pAdj **numeric** \cr
#' Vector of p-values after FDR-correction.
#' @param fc **numeric** \cr
#' Vector of fold changes.
#' @param thresFC **numeric(1)** \cr
#' Threshold for the fold changes.
#' @param thresP **numeric(1)** \cr
#' Threshold for the p-values.
#'
#' @return A factor vector with the three significance categories.
#'
#' @importFrom checkmate assertNumeric
#' @importFrom dplyr case_when
.calcSignCat_ttest <- function(p, pAdj, fc, thresFC = 2, thresP = 0.05) {

  checkmate::assertNumeric(p, lower = 0, upper = 1)
  checkmate::assertNumeric(pAdj, lower = 0, upper = 1, len = length(p))
  checkmate::assertNumeric(fc, len = length(p))
  checkmate::assertNumeric(thresFC, len = 1)
  checkmate::assertNumeric(thresP, len = 1, lower = 0, upper = 1)

  ind_sign_FDR <- (pAdj <= thresP & p <= thresP &
                (fc >= thresFC | fc <= 1/thresFC) & !is.na(p))
  ind_sign <- (pAdj > thresP & p <= thresP &
                 (fc >= thresFC | fc <= 1/thresFC) & !is.na(p))
  ind_nosign <- (p > thresP | (fc < thresFC & fc > 1/thresFC)) & !is.na(p)

  significance <- dplyr::case_when(
    ind_sign_FDR ~ "significant after FDR correction",
    ind_sign ~ "significant",
    ind_nosign ~ "not significant",
    is.na(p) ~ NA_character_
  )

  significance <- factor(significance,
   levels = c("not significant", "significant",
              "significant after FDR correction"))

  return(significance)
}



#' Calculate significance categories for ANOVA.
#'
#'
#' This function groups all proteins into the following three categories based
#' on the p-values (with and without FDR correction, ANOVA and posthoc) and
#' fold changes:
#' 1) not significant (p.anova > thresP or p_posthoc > thresPfc < thresFC)
#' 2) significant (p < thresP and fc > thresFC, but pAdj > thresP)
#' 3) significant after FDR-correction (pAdj < thresP and fc > thresFC)
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
