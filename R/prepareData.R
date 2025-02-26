#' Prepare proteomics data for analysis.
#'
#' @param data_path               \strong{character} \cr
#'                                The path to an .xlsx file containing the input data.
#' @param intensity_columns       \strong{integer vector} \cr
#'                                The numbers of the intensity columns in the table.
#' @param na_strings              \strong{character} \cr
#'                                A vector containing the symbols to be recognized as missing values (with the exception of 0).
#' @param zero_to_NA              \strong{logical} \cr
#'                                If \code{TRUE}, 0 will be treated as missing value.
#' @param do_log_transformation   \strong{logical} \cr
#'                                If \code{TRUE}, the data will be log-transformed.
#' @param log_base                \strong{numeric} \cr
#'                                The base used, if \code{do_log_transformation = TRUE}.
#' @param use_groups              \strong{logical} \cr
#'                                If \code{TRUE}, group information encoded in the column names is used.
#' @param group_colours           \strong{character vector} \cr
#'                                The hex codes for the group colors, if the data has groups. If \code{NULL}, a default color scale will be used.
#' @param normalization           \strong{character} \cr
#'                                The method of normalization. Options are "nonorm" (no normalization), "median", "loess", "quantile" or "lts" normalization.
#' @param lts_quantile            \strong{numeric} \cr
#'                                The quantile for the lts normalization if \code{normalization = "lts"}.
#'
#' @return A list containing the prepared data and the ids of the data as data.frames as well as the groups, number of groups and group colors.
#' @export
#'
#' @examples
#'

prepareData <- function (data_path,
                         intensity_columns,
                         na_strings = c("NA", "NaN", "Filtered","#NV"),
                         zero_to_NA = TRUE,
                         do_log_transformation = TRUE, log_base = 2,
                         use_groups = FALSE, group_colours = NULL,
                         normalization = "loess", lts_quantile = 0.8){


  #### read and prepare data file ####

  D <- openxlsx::read.xlsx(data_path, na.strings = na_strings)
  mess = ""

  id <- D[, -intensity_columns]
  D <- D[, intensity_columns]

  if(zero_to_NA) {
    D[D == 0] <- NA
    mess <- paste0(mess, "Zeros set to NA. \n")
  }

  if(do_log_transformation) {
    D <- log(D, base = log_base)
    mess <- paste0(mess, "Log-transformation with base ", log_base ,". \n")
  }


  #### make data groups ####

  if (use_groups) {
    group <- factor(limma::strsplit2(colnames(D), "_")[,1])
    mess <- paste0(mess, "Groups used. \n")
  } else {
    group <- NULL
    mess <- paste0(mess, "Groups not used. \n")
  }

  nr_groups <- length(levels(group))

  if (is.null(group_colours) & nr_groups >= 1) group_colours <- scales::hue_pal()(nr_groups)


  #### normalize the data ####

  D <- automatedNormalization(DATA = D, method = normalization, is_log_transformed = do_log_transformation, log_base = log_base, lts.quantile = lts_quantile)

  mess <- paste0(mess, D$message)
  D <- D$data


  #### calculate long form ####

  D_long <- tidyr::pivot_longer(data = D, cols = 1:ncol(D))
  if (use_groups) {
    D_long$group <- factor(limma::strsplit2(D_long$name, "_")[,1])
  } else {
    D_long$group <- NA
  }

  ### add column with sample number
  if (use_groups) {
    D_long$sample <- limma::strsplit2(D_long$name, "_")[,2]
  } else {
    D_long$sample <- NA
  }

  message(mess)

  return (list("D" = D, "ID" = id, "D_long" = D_long, "group" = group, "number_groups" = nr_groups, "group_colors" = group_colours, "message" = mess))
}
