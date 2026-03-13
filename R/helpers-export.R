#' Export SummarizedExperiment object to an Excel file, with one assay per sheet.
#'
#' @param SE \strong{SummarizedExperiment object} \cr Object to be exported.
#'      E.g. the result of \code{\link{combineComparisons}}.
#' @param file \strong{character} \cr file path to the output Excel file.
#'      If the file already exists, it will be overwritten.
#'
#' @returns nothing, but an Excel file is written containing the assay matrices.
#' Each assay is written to a separate sheet, and the sheet name corresponds to
#' the assay name. Row names of the assay matrices are included in the Excel file.
#' @export
#'
#' @examples
#' # TODO
# TODO: allow csv, tsv etc as output
exportSE <- function(SE, file) {
  
  wb <- openxlsx::createWorkbook()
  
  cData <- SummarizedExperiment::colData(SE)
  
  for (assay_name in SummarizedExperiment::assayNames(SE)) {
    openxlsx::addWorksheet(wb, assay_name)
    mat <- as.data.frame(SummarizedExperiment::assay(SE, assay_name))
    mat <- cbind(cData, mat)
    openxlsx::writeData(wb, sheet = assay_name, x = mat, rowNames = TRUE,
                        keepNA = TRUE)
  }

  openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
  return(invisible(NULL))
}

