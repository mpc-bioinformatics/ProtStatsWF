#' Prepare data from an xlsx sheet
#'
#'
#' @param data_path              A character containing the path to an .xlsx file.
#' @param intensity_columns      An integer vector containing the intensity columns of the table.
#' 
#' @return A list containing an intensity data.frame, an IDs data frame and a factor of the sample groups
#'
#' @examples
#' \dontrun{
#' in_path <- "/Users/thisuser/Documents/dataFolder/data.xlsx"
#' int_cols <- 3:8
#'
#' result <- workflow_ttest(data_path = in_path, intensity_columns = int_cols)
#'}

prepareTtestData <- function(data_path,
                             intensity_columns
){
  
  D <- openxlsx::read.xlsx(data_path, na.strings = c("NA", "NaN", "Filtered","#NV"))
  
  id <- D[, -intensity_columns]
  D <- D[, intensity_columns]
  
  D[D == 0] <- NA
  
  group <- factor(limma::strsplit2(colnames(D), "_")[,1])
  number_of_groups <- length(levels(group))
  
  sample <- factor(limma::strsplit2(colnames(D), "_")[,2])
  
  return(list("D" = D, "ID" = id, "group" = group, "number_of_groups" = number_of_groups, "sample" = sample))
}



#' Calculation of significance categories for a volcano plot for the t-test
#'
#' This function groups all proteins into the following three catrgories based
#' on the p-values (with and without FDR correction) and fold changes:
#' 1) not signifcant (p > thres_p or fc < thres_fc)
#' 2) significant (p < thres_p and fc > thres_fc, but p_adj > thres_p)
#' 3) significant after FDR-correction (p_adj < thres_p and fc > thres_fc)
#'
#' @param p vector of p-values before FDR-correction
#' @param p_adj vector of p-values after FDR-correction
#' @param fc vector of fold changes
#' @param thres_fc threshold for the fold changes
#' @param thres_p threshold for the p-values
#'
#' @return a factor with the three significance categories
#' @export
#'
#' @examples # TODO
calculate_significance_categories_ttest <- function(p, p_adj, fc, thres_fc=2, thres_p=0.05) {
  
  significance <- dplyr::case_when(
    p_adj <= thres_p & p <= thres_p & (fc >= thres_fc | fc <= 1/thres_fc) & !is.na(p) ~ "significant after FDR correction",
    p_adj > thres_p & p <= thres_p & (fc >= thres_fc | fc <= 1/thres_fc) & !is.na(p) ~ "significant",
    (p > thres_p | (fc < thres_fc & fc > 1/thres_fc)) & !is.na(p) ~ "not significant",
    is.na(p) ~ NA_character_
  )
  
  significance <- factor(significance, levels = c("not significant", "significant", "significant after FDR correction"))
  
  return(significance)
}



#' Calculate significance categories for ANOVA
#'
#' @param p_posthoc vector of posthoc p-values 
#' @param p_anova_adj vector of p-values after FDR-correction
#' @param p_anova vector of p-values before FDR-correction
#' @param fc vector of fold changes
#' @param thres_fc threshold for the fold changes
#' @param thres_p threshold for the p-values
#'
#' @return A factor containing the significances
#' @export
#'
#' @examples
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