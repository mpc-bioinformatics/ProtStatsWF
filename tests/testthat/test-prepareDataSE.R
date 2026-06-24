
test_that("Prepare data from csv file", {
  file_proteins <- system.file("extdata", "proteins_HCC.csv", package = "ProtStatsWF")
  file_clinical <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")
  D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
                proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
                sampleNameColumn = "Sample", fileType = "csv")

  expect_snapshot(D$D_long)
  expect_snapshot(D$SE)
  expect_snapshot(SummarizedExperiment::assays(D$SE)$intensity)
  expect_snapshot(SummarizedExperiment::assays(D$SE)$intensity_norm)
  expect_snapshot(SummarizedExperiment::colData(D$SE))
  expect_snapshot(SummarizedExperiment::rowData(D$SE))

})

# TODO: test different file formats and settings


