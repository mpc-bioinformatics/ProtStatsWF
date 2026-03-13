
test_that("Data preparation with groups (median normalization) for test_file_1", {
  dataPath <- system.file("extdata", "test_file_2.xlsx", package = "ProtStatsWF")
  sampleInfoPath <- system.file("extdata", "test_file_1_sampleInfo.xlsx", package = "ProtStatsWF")

  prepData <- prepareDataSE(dataPath = dataPath, 
                          intensityColumns = 3:11,
                          proteinNameColumn = "peptide",
                          normMethod = "median", 
                          sampleInfoPath = sampleInfoPath, 
                          verbose = FALSE)
  
  expect_snapshot(prepData$SE)
  expect_snapshot(SummarizedExperiment::colData(prepData$SE))
  expect_snapshot(SummarizedExperiment::rowData(prepData$SE))
  expect_snapshot(SummarizedExperiment::assays(prepData$SE)$intensity)
  expect_snapshot(SummarizedExperiment::assays(prepData$SE)$intensity_norm)
  expect_snapshot(print(prepData$D_long, n = 200))
})

test_that("Data preparation without groups (loess normalization) for test_file_1", {
  dataPath <- system.file("extdata", "test_file_1.xlsx", package = "ProtStatsWF")
 
  prepData <- prepareDataSE(dataPath = dataPath, 
                          intensityColumns = 3:11, 
                          proteinNameColumn = "peptide",
                          normMethod = "quantile",
                          sampleInfoPath = NULL, 
                          verbose = FALSE)
  
  expect_snapshot(prepData$SE)
  expect_snapshot(SummarizedExperiment::colData(prepData$SE))
  expect_snapshot(SummarizedExperiment::rowData(prepData$SE))
  expect_snapshot(SummarizedExperiment::assays(prepData$SE)$intensity)
  expect_snapshot(SummarizedExperiment::assays(prepData$SE)$intensity_norm)
  expect_snapshot(print(prepData$D_long, n = 200))
})