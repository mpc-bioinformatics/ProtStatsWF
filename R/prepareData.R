#' Data preparation
#'
#' Prepare quantitative proteomics data and combine it with sample information
#' (optional, e.g., clinical data) in a SummarizedExperiment object.
#'
#' @param dataPath **character(1)** \cr
#' Path to the data file (xlsx, csv, tsv, or txt format) containing the
#' quantitative proteomics data. The file should have at least a column with
#' protein names and columns with intensity values for each sample.
#' @param intensityColumns **integer** \cr
#' Column numbers in data that contain the intensity values. E.g. 1:12 if
#' the first 12 columns contain intensity values for the 12 samples.
#' @param proteinNameColumn **character(1)** \cr
#' Column name in data that contains the protein names. Default is "Protein".
#' @param sampleInfoPath **character(1)** \cr
#' Path to the file containing sample information (optional). Default is NULL
#' (no sample info file).
#' @param sampleNameColumn **character(1)** \cr
#' Column name of the sampleInfo that contains the sample names. Default is
#' "Sample". Those names must correspond to column names in data.
#' @param zeroToNA **logical(1)** \cr
#' If TRUE (default), zero intensity values are converted to NA.
#' @param doLogTrans **logical(1)** \cr
#' If TRUE (default), intensity values will be log-transformed. As it is highly
#' recommended to log-transform proteomics intensities, only set to FALSE if the
#' data are already log-transformed or you have a special case like TMT ratios.
#' @param logBase **numeric(1)** \cr
#' Base for log transformation. Default is 2.
#' @param normMethod **character(1)** \cr
#' Normalization method, default is "loess". Options are "nonorm", "median",
#' "quantile", "loess" and "lts".
#' @param ltsQuantile **numeric(1)** \cr
#' Quantile to use for the least trimmed squares regression in the "lts"
#' normalization method. Default is 0.8.
#' @param sep **character(1)** \cr
#' Separator for csv, tsv, or txt files. Default is "," (comma).
#' Ignored for xlsx files.
#' @param dec **character(1)** \cr
#' Decimal point character for csv, tsv, or txt files. Default is "." (dot).
#' Ignored for xlsx files.
#' @param header **logical(1)** \cr
#' If TRUE (default), the first row is treated as column names. Ignored
#' for xlsx files, where the header is always read.
#' @param sheet **integer(1)** \cr
#' Excel sheet number to read in for dataPath. Only used if fileType = "xlsx".
#' @param NAStrings **character** \cr
#' Character vector of strings to interpret as NA when reading the data file.
#' Default is c("NA", "NaN", "Filtered","#NV", "").
#' @param verbose **logical(1)** \cr
#' If TRUE (default), messages about the data processing steps are printed.
#'
#' @returns List with two elements:
#' \itemize{
#'  \item SE: A SummarizedExperiment object containing the normalized intensity
#'  data in the assay "intensity_norm", the original intensity data in the assay
#'  "intensity", the protein information in rowData (all other columns in data
#'  that are not the proteinName or intensity columns), and the sample
#'  information in colData.
#'  \item D_long: A data frame in long format containing the normalized
#'  intensity values for each protein and sample, along with the corresponding
#'  sample information from colData.
#'  }
#' @export
#'
#' @importFrom checkmate assert assertCharacter assertFileExists assertFlag
#' @importFrom checkmate assertIntegerish assertNumeric assertSubset
#' @importFrom checkmate checkFileExists
#' @importFrom dplyr left_join select
#' @importFrom openxlsx read.xlsx
#' @importFrom SummarizedExperiment assay colData rowData SummarizedExperiment
#' @importFrom tidyr pivot_longer
#' @importFrom tools file_ext
#' @importFrom utils read.table
#'
#'
#' @examples
#' file_proteins <- system.file("extdata", "proteins_HCC.csv",
#'   package = "ProtStatsWF")
#' file_clinical <- system.file("extdata", "clinical_data.csv",
#'   package = "ProtStatsWF")
#' prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
#'   proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
#'   sampleNameColumn = "Sample", verbose = FALSE)
prepareDataSE <- function(dataPath,
                          intensityColumns,
                          proteinNameColumn = "Protein",
                          sampleInfoPath = NULL,
                          sampleNameColumn = "Sample",
                          zeroToNA = TRUE,
                          doLogTrans = TRUE,
                          logBase = 2,
                          normMethod = "loess",
                          ltsQuantile = 0.8,
                          sep = ",",
                          dec = ".",
                          header = TRUE,
                          sheet = 1,
                          NAStrings = c("NA", "NaN", "Filtered","#NV", ""),
                          verbose = TRUE) {
  checkmate::assertFileExists(dataPath)
  if (!is.null(sampleInfoPath))
    checkmate::assert(checkmate::checkFileExists(sampleInfoPath))
  checkmate::assertFlag(zeroToNA)
  checkmate::assertFlag(doLogTrans)
  checkmate::assertNumeric(logBase, len = 1)
  checkmate::assertCharacter(sep, len = 1)
  checkmate::assertCharacter(dec, len = 1)
  checkmate::assertFlag(header)
  checkmate::assertIntegerish(sheet)
  checkmate::assertCharacter(NAStrings)
  checkmate::assertFlag(verbose)

  filetype_data <- tools::file_ext(dataPath)
  checkmate::assertSubset(filetype_data,
                          choices = c("csv", "tsv", "txt", "xlsx"))
  if (filetype_data %in% c("csv", "txt", "tsv")) {
    D_complete <- utils::read.table(dataPath, sep = sep, header = header,
                                    dec = dec, quote = "\"",
                                    na.strings = NAStrings)
  }
  if (filetype_data == "xlsx") {
    D_complete <- openxlsx::read.xlsx(dataPath, colNames = header,
                                      sheet = sheet, na.strings = NAStrings)
  }

  if (!is.null(sampleInfoPath)) {
    filetype_sampleInfo <- tools::file_ext(sampleInfoPath)
    checkmate::assertSubset(filetype_sampleInfo,
                            choices = c("csv", "tsv", "txt", "xlsx"))
    if (filetype_sampleInfo %in% c("csv", "txt", "tsv")) {
      sampleInfo <- utils::read.table(sampleInfoPath, sep = sep,
                                      header = header, dec = dec,
                                      quote = "\"", na.strings = NAStrings)
    }
    if (filetype_sampleInfo == "xlsx") {
      sampleInfo <- openxlsx::read.xlsx(sampleInfoPath, colNames = TRUE,
                                        sheet = 1, na.strings = NAStrings)
    }
    checkmate::assertSubset(sampleNameColumn, choices = colnames(sampleInfo))
    if (verbose) message("Sample information file read in.")
  }

  checkmate::assertIntegerish(intensityColumns, any.missing = FALSE, lower = 1,
                              upper = ncol(D_complete),
                              max.len = ncol(D_complete),
                              min.len = 1, unique = TRUE)
  checkmate::assertSubset(proteinNameColumn, choices = colnames(D_complete))

  id <- D_complete[, -intensityColumns]
  D <- D_complete[, intensityColumns]
  rownames(id) <- id[, proteinNameColumn]
  rownames(D) <- id[, proteinNameColumn]

  if (zeroToNA) {
    D[D == 0] <- NA
    if (verbose) message("Zeros set to NA.")
  }

  if (doLogTrans) {
    D <- log(D, base = logBase)
    if (verbose) message("Log-transformation with base ", logBase, ".")
  }

  D_norm <- .normalization(DATA = D, method = normMethod,
                              is_log_transformed = TRUE, ## we assume that data is log-transformed!
                              log_base = logBase, lts.quantile = ltsQuantile,
                              verbose = verbose)

  if (!is.null(sampleInfoPath)) {
    ind <- match(sampleInfo[, sampleNameColumn], colnames(D))
    D_norm <- D_norm[,ind] # sort columns of D like sampleInfo
    D <- D[, ind]
  } else {
    sampleInfo <- data.frame(SampleName = colnames(D))
    rownames(sampleInfo) <- colnames(D)
  }

  assays <- list(intensity_norm = as.matrix(D_norm), intensity = as.matrix(D))
  SE <- SummarizedExperiment::SummarizedExperiment(assays = assays,
                                                   rowData = id,
                                                   colData = sampleInfo)
  D_norm_long <- as.data.frame(SummarizedExperiment::assay(SE,
                                                           "intensity_norm"))
  D_norm_long$.feature <- rownames(D_norm_long)
  D_long <- tidyr::pivot_longer(D_norm_long, cols = -".feature",
    names_to = ".sample", values_to = "intensity_norm")
  col_data <- as.data.frame(SummarizedExperiment::colData(SE))
  col_data$.sample <- rownames(col_data)
  D_long <- dplyr::left_join(D_long, col_data, by = ".sample")
  row_data <- as.data.frame(SummarizedExperiment::rowData(SE))
  row_data$.feature <- rownames(row_data)
  D_long <- dplyr::left_join(D_long, row_data, by = ".feature")
  # bring factor levels in same order as in sampleInfo
  D_long$.sample <- factor(D_long$.sample,
                           levels = sampleInfo[, sampleNameColumn])
  D_long <- as.data.frame(D_long)

  return(list(SE = SE, D_long = D_long))
}







#' Pivot SE to long format
#'
#'helper function to extract assay data and pivot to long format.
#'Also adds the information from colData and rowData.
#'
#' @param SE A SummarizedExperiment object.
#' @param cols Columns to pivot.
#'
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom tidyr pivot_longer
#'
#' @returns A data frame in long format.
#'
#' @examples
.pivot_longer_SE <- function(SE, cols) {

  # Extract assay
  D <- as.data.frame(SummarizedExperiment::assay(SE))
  D$id. <- rownames(D)

  D_long <- tidyr::pivot_longer(D, -id., names_to = "sample.",
                                values_to = "value.")

  # Add colData
  col_data <- as.data.frame(SummarizedExperiment::colData(SE))
  col_data$sample. <- rownames(col_data)
  D_long <- dplyr::left_join(D_long, col_data, by = "sample.")

  # Add rowData
  row_data <- as.data.frame(SummarizedExperiment::rowData(SE))
  row_data$id. <- rownames(row_data)
  D_long <- dplyr::left_join(D_long, row_data, by = "id.")

}









