#' Prepare proteomics data for analysis.
#'
#' @template param_dataPath
#' @param intensityColumns **integer** \cr The columns numbers containing protein intensities in the data set.
#' @param proteinNameColumn **character(1)** \cr The name of the column in the data file that contains the protein names. Default is "Protein".
#' @param sampleInfoPath **character(1)** \cr The path to the sample information file with group information.
#'    This file must contain a column corresponding to the column names of the data file (only intensity columns).
#'    Further columns are optional and may give information on groups, timepoints or patient information (patient ID, age, gender etc.).
#' @param sampleNameColumn **character(1)** \cr The name of the column in the sample information file that contains the sample names.
#'    The sample names have to correspond to the columns names of the data file (only for intensity columns),
#'    but they do not necessarily have to be in the same order as the columns. They will be used to match the data file and the sample information file.
# @param groupColumn **character(1)** \cr The name of the column in the sample information file that contains the group information.
#    Will be used for colouring in the plots. If NULL, no colouring is done.
# @param group2Column **character(1)** \cr The name of the column in the sample information file that contains a second grouping variable, e.g. timepoints.
#    This will be used in the PCA plot for "shape" of the data points.
# @param sampleIDColumn **character(1)** \cr The name of the column in the sample information file that contains the sample IDs.
#    This is especially important if a paired analysis is going to be used (paired t-test, repeated measures ANOVA),
#    as the sample IDs will be used to link the paired samples.
#'
#' @param doLogTrans **logical(1)** \cr If \code{TRUE}, the data will be log-transformed.
#' @param logBase **numeric(1)** \cr The base for the logarithm, if \code{do_log_transformation = TRUE}. Default is 2.
#' @param normMethod **character(1)** \cr The method of normalization. Options are "nonorm" (no normalization), "median", "loess", "quantile" or "lts" normalization.
#' @param ltsQuantile **numeric(1)** The quantile for the lts normalization if \code{normalization = "lts"}.
#'
#' @param fileType **character(1)** \cr Type of input file: "csv" or "tsv" or "txt" or "xlsx".
#' @param sep **character(1)** \cr The field separator, e.g. " " for blanks, "," for comma or "\\t" for tab. Default is ",".
#' @param dec **character(1)** \cr Decimal separator, e.g. "," for comma or "." for dot. Default is ".".
#' @param header **logical(1)** \cr If TRUE, first line is counted as column names.
#' @param sheet **integer(1)** \cr Sheet number (only needed for xlsx files, default is to use the first sheet).
#' @param zeroToNA **logical(1)** \cr If \code{TRUE}, 0 will be treated as missing value.
#' @param NAStrings **character** \cr A vector containing the symbols to be recognized as missing values (except 0).
#' @param verbose **logical(1)** \cr If \code{TRUE}, messages about the performed steps will be printed in the console.
#'
#' @return A list containing the prepared data and the ids of the data as data.frames as well as the groups, number of groups and group colors.
#' @export
#'
#' @examples
#'
prepareDataSE <- function(dataPath,
                          intensityColumns,
                          proteinNameColumn = "Protein",
                          sampleInfoPath = NULL,
                          sampleNameColumn = "SampleName",
                          doLogTrans = TRUE,
                          logBase = 2,
                          normMethod = "loess",
                          ltsQuantile = 0.8,

                          fileType = "xlsx",
                          sep = ",",
                          dec = ".",
                          header = TRUE,
                          sheet = 1,
                          zeroToNA = TRUE,
                          NAStrings = c("NA", "NaN", "Filtered","#NV"),
                          verbose = TRUE) {


  ### import data file

  if (fileType == "csv" | fileType == "txt" | fileType == "tsv") {
    if (fileType == "csv") {
      sep <- ","
    } else if (fileType == "tsv") {
      sep <- "\t"
    }
    D_complete <- utils::read.table(dataPath, sep = sep, header = header, dec = dec,
                                    quote = "\"")
    if (!is.null(sampleInfoPath)) {
      sampleInfo <- utils::read.table(sampleInfoPath, sep = sep, header = header, dec = dec,
                                      quote = "\"")
    }
  }
  if (fileType == "xlsx") {
    D_complete <- openxlsx::read.xlsx(dataPath, colNames = header, sheet = sheet)
    if (!is.null(sampleInfoPath)) {
      sampleInfo <- openxlsx::read.xlsx(sampleInfoPath, colNames = TRUE, sheet = 1)
    }
  }

  id <- D_complete[, -intensityColumns]
  D <- D_complete[, intensityColumns] # only intensity columns
  rownames(id) <- id[, proteinNameColumn]
  rownames(D) <- id[, proteinNameColumn]

  #ID <<- id

  ## data preprocessing (NAs, log, normalization)
  if (zeroToNA) {
    D[D == 0] <- NA
    if(verbose) message("Zeros set to NA.")
  }

  if (doLogTrans) {
    D <- log(D, base = logBase)
    if (verbose) message("Log-transformation with base ", logBase, ".")
  }

  D_norm <- automatedNormalization(DATA = D, method = normMethod,
                              is_log_transformed = doLogTrans,
                              log_base = logBase, lts.quantile = ltsQuantile,
                              verbose = verbose)

  ### TODO: add option to save normalized data as xlsx file
  # if (outType == "xlsx") {
  #   exportSE(prepared_data$SE, file = file.path(output_path, paste0("D_norm", suffix, ".xlsx")))
  # }



  ### TODO: check if all samples can be found in sampleInfo
  ### TODO: remove additional samples from sampleInfo

  ## read in sample information file, if given
  if (!is.null(sampleInfoPath)) {
    #sampleInfo <- openxlsx::read.xlsx(sampleInfoPath, colNames = TRUE)
    ind <- match(sampleInfo[, sampleNameColumn], colnames(D))
    D_norm <- D_norm[,ind] # sort columns of D like sampleInfo

    D <- D[, ind]
    if (verbose) message("Sample information file read in.")
  } else {
    sampleInfo <- data.frame(SampleName = colnames(D))
    rownames(sampleInfo) <- colnames(D)
    if(verbose) message("No sample information file given.")
  }

  #SI <<- sampleInfo

  #if (normMethod != "nonorm") {
    assays <- list(intensity_norm = as.matrix(D_norm), intensity = as.matrix(D))
  #} else {
  #  assays <- list(intensity = as.matrix(D))
  #}

  # TODO: include metadata, e.g. about normalization method?
  #if (is.null(sampleInfo)) {
    SE <- SummarizedExperiment::SummarizedExperiment(assays = assays,
                                                     rowData = id,
                                                     colData = sampleInfo)
  #} else {
  #  SE <- SummarizedExperiment::SummarizedExperiment(assays = assays,
  #                                                   colData = sampleInfo,
  #                                                   rowData = id)
  #}

  ### long format:
  # TODO: this is a generic function from tidySummarizedExperiments package
  # TODO: D_long is a tibble, will that cause problems?
  suppressMessages({
  D_long <- tidySummarizedExperiment:::pivot_longer.SummarizedExperiment(SE,
                  cols = proteinNameColumn)
  })
  D_long <- dplyr::select(D_long, -c("name", "value"))
  # bring factor levels in same order as in sampleInfo
  D_long$.sample <- factor(D_long$.sample, levels = sampleInfo[, sampleNameColumn])

  #print(sampleInfo[, sampleNameColumn])
  #print(levels(D_long$.sample))

  return(list(SE = SE, D_long = D_long))
}




### helper function to extract assay data and pivot to long format.
### Also add the information from colData
pivot_longer_SE <- function(SE, cols) {

  # Extract assay
  D <- as.data.frame(SummarizedExperiment::assay(SE))
  D$id. <- rownames(D)

  D_long <- tidyr::pivot_longer(D, -id., names_to = "sample.", values_to = "value.")

  # Add colData
  col_data <- as.data.frame(SummarizedExperiment::colData(SE))
  col_data$sample. <- rownames(col_data)
  D_long <- dplyr::left_join(D_long, col_data, by = "sample.")

  # Add rowData
  row_data <- as.data.frame(SummarizedExperiment::rowData(SE))
  row_data$id. <- rownames(row_data)
  D_long <- dplyr::left_join(D_long, row_data, by = "id.")

}




























#' Prepare proteomics data for analysis.
#'
#' @param data_path               \strong{character} \cr
#'                                The path to an .xlsx file containing the input data.
#' @param filetype **character(1)** \cr Type of input file: "csv" or "tsv" or "txt" or "xlsx".
#' @param sep **character(1)** \cr The field separator, e.g. " " for blanks, "," for comma or "\\t" for tab. Default is ",".
#' @param dec **character(1)** \cr Decimal separator, e.g. "," for comma or "." for dot. Default is ".".
#' @param header **logical(1)** \cr If TRUE, first line is counted as column names.
#' @param sheet **integer(1)** \cr Sheet number (only needed for xlsx files, default is to use the first sheet).
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
                         filetype = "xlsx",
                         sep = ",",
                         dec = ".",
                         header = TRUE,
                         sheet = 1,
                         intensity_columns,
                         na_strings = c("NA", "NaN", "Filtered","#NV"),
                         zero_to_NA = TRUE,
                         do_log_transformation = TRUE, log_base = 2,
                         use_groups = FALSE, group_colours = NULL,
                         normalization = "loess", lts_quantile = 0.8){


  #### read and prepare data file ####


  if (filetype == "csv" | filetype == "txt" | filetype == "tsv") {
    if (filetype == "csv") {
      sep <- ","
    } else if (filetype == "tsv") {
      sep <- "\t"
    }
    D <- utils::read.table(data_path,
                           sep = sep,
                           header = header,
                           dec = dec,
                           quote = "\"")
  }
  if (filetype == "xlsx") {
    D <- openxlsx::read.xlsx(data_path, colNames = header, sheet = sheet)
  }


  #D <- openxlsx::read.xlsx(data_path, na.strings = na_strings)
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
