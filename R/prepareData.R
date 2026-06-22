
#' Prepare quantitative proteomics data and combine it with sample information
#' (e.g., clinical data) in a SummarizedExperiment object.
#'
#' @param dataPath **character(1)** \cr Path to the data file (xlsx, csv, tsv, or txt) containing the quantitative proteomics data. The file should have a column with protein names and columns with intensity values for each sample.
#' @param intensityColumns **integer** \cr Column numbers that contain the intensity values. E.g. 1:12 if the first 12 columns contain intensity values for 12 samples.
#' @param proteinNameColumn **character(1)** \cr Column name that contains the protein names. Default is "Protein".
#' @param sampleInfoPath **character(1)** \cr Path to the file containing sample information. Default is NULL (no sample info file).
#' @param sampleNameColumn **character(1)** \cr Column name of the sampleInfo that contains the sample names. Default is "SampleName".
#' @param doLogTrans **logical(1)** \cr Whether to perform log transformation. Default is TRUE.
#' @param logBase **numeric(1)** \cr Base for log transformation. Default is 2.
#' @param normMethod **character(1)** \cr Normalization method to use. Default is "loess". See [automatedNormalization()] for options.
#' @param ltsQuantile **numeric(1)** \cr Quantile to use for the least trimmed squares regression in the "lts" normalization method. Default is 0.8.
#' @param fileType **character(1)** \cr Type of the data file. One of "xlsx", "csv", "tsv", or "txt". Default is "xlsx".
#' @param sep **character(1)** \cr Separator for csv, tsv, or txt files. Default is "," (comma). Ignored for xlsx files.
#' @param dec **character(1)** \cr Decimal point character for csv, tsv, or txt files. Default is "." (dot). Ignored for xlsx files.
#' @param header **logical(1)** \cr Whether the data file has a header row. Default is TRUE. Ignored for xlsx files, where the header is always read.
#' @param sheet **integer(1)** \cr Sheet number to read from an xlsx file. Default is 1. Ignored for csv, tsv, or txt files.
#' @param zeroToNA **logical(1)** \cr Whether to convert zero intensity values to NA. Default is TRUE.
#' @param NAStrings **character** \cr Character vector of strings to interpret as NA when reading the data file. Default is c("NA", "NaN", "Filtered","#NV", "").
#' @param verbose **logical(1)** \cr Whether to print messages about the data processing steps. Default is TRUE.
#'
#' @returns List with two elements:
#' \itemize{
#'  \item SE: A SummarizedExperiment object containing the normalized intensity data in
#'  the assay "intensity_norm", the original intensity data in the assay "intensity",
#'  the protein information in rowData (all other columns in data that are not the proteinName or
#'  intensity columns), and the sample information in colData..
#'  \item D_long: A data frame in long format containing the normalized intensity values
#'  for each protein and sample, along with the corresponding sample information from colData.
#'  }
#' @export
#'
#' @examples
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
                          NAStrings = c("NA", "NaN", "Filtered","#NV", ""),
                          verbose = TRUE) {


  ### import data file

  if (fileType == "csv" | fileType == "txt" | fileType == "tsv") {
    if (fileType == "csv") {
      sep <- ","
    } else if (fileType == "tsv") {
      sep <- "\t"
    }
    D_complete <- utils::read.table(dataPath, sep = sep, header = header, dec = dec,
                                    quote = "\"", na.strings = NAStrings)
    if (!is.null(sampleInfoPath)) {
      sampleInfo <- utils::read.table(sampleInfoPath, sep = sep, header = header, dec = dec,
                                      quote = "\"", na.strings = NAStrings)
    }
  }
  if (fileType == "xlsx") {
    D_complete <- openxlsx::read.xlsx(dataPath, colNames = header, sheet = sheet,
                                      na.strings = NAStrings)
    if (!is.null(sampleInfoPath)) {
      sampleInfo <- openxlsx::read.xlsx(sampleInfoPath, colNames = TRUE,
                                        sheet = 1, na.strings = NAStrings)
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
  suppressMessages({
  D_long <- tidySummarizedExperiment:::pivot_longer.SummarizedExperiment(SE,
                  cols = proteinNameColumn)
  })
  D_long <- dplyr::select(D_long, -c("name", "value"))
  # bring factor levels in same order as in sampleInfo
  D_long$.sample <- factor(D_long$.sample, levels = sampleInfo[, sampleNameColumn])
  D_long <- as.data.frame(D_long)

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









