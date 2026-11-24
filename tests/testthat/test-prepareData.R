
test_that("Prepare data from csv file", {
  file_proteins <- system.file("extdata", "proteins_HCC.csv", package = "ProtStatsWF")
  file_clinical <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")
  SE <- prepareData(dataPath = file_proteins, intensityColumns = 6:43,
                proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
                sampleNameColumn = "Sample", verbose = FALSE)

  #expect_snapshot(D$D_long)
  expect_snapshot(SE)
  expect_snapshot(SummarizedExperiment::assays(SE)$intensity)
  expect_snapshot(SummarizedExperiment::assays(SE)$intensity_norm)
  expect_snapshot(SummarizedExperiment::colData(SE))
  expect_snapshot(SummarizedExperiment::rowData(SE))

})

# TODO: test different file formats and settings



test_that("Test prepareData if sampleInfo describes data points not present in dataset", {
  file_proteins <- system.file("extdata", "proteins_HCC_small_missingSamples.csv", package = "ProtStatsWF")
  file_clinical <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")

  expect_warning({
  SE <- prepareData(dataPath = file_proteins, intensityColumns = 6:41,
                   proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
                   sampleNameColumn = "Sample", verbose = FALSE)
  })

  #expect_snapshot(D$D_long)
  expect_snapshot(SE)
  expect_snapshot(SummarizedExperiment::assays(SE)$intensity)
  expect_snapshot(SummarizedExperiment::assays(SE)$intensity_norm)
  expect_snapshot(SummarizedExperiment::colData(SE))
  expect_snapshot(SummarizedExperiment::rowData(SE))

})


