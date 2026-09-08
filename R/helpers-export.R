#' Export SummarizedExperiment object to an Excel file, with one assay per sheet.
#'
#' @param SE \strong{SummarizedExperiment object} \cr Object to be exported.
#' @param file \strong{character} \cr file path to the output Excel file.
#'      If the file already exists, it will be overwritten.
#'
#' @returns nothing, but an Excel file is written containing the assay matrices.
#' Each assay is written to a separate sheet, and the sheet name corresponds to
#' the assay name. Row names of the assay matrices are included in the Excel file.
#' @export
#'
#' @examples
#' file_proteins <- system.file("extdata", "proteins_HCC.csv",
#'   package = "ProtStatsWF")
#' file_clinical <- system.file("extdata", "clinical_data.csv",
#'   package = "ProtStatsWF")
#' D_hcc <- prepareData(file_proteins, intensityColumns = 6:43,
#'   proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
#'   sampleNameColumn = "Sample", verbose = FALSE)
#' exportSE(D_hcc$SE, file = file.path(tempdir(), "HCC_data.xlsx"))
exportSE <- function(SE, file) {
  
  wb <- openxlsx::createWorkbook()
  
  rData <- as.data.frame(SummarizedExperiment::rowData(SE))
  
  for (assay_name in SummarizedExperiment::assayNames(SE)) {
    openxlsx::addWorksheet(wb, assay_name)
    mat <- as.data.frame(SummarizedExperiment::assay(SE, assay_name))
    mat <- cbind(rData, mat)
    openxlsx::writeData(wb, sheet = assay_name, x = mat, rowNames = TRUE,
                        keepNA = TRUE)
  }

  openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
  return(invisible(NULL))
}

