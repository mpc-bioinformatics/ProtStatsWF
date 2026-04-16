
test_that("PCA plot for test_file_2 (no imputation, no shape, no colour, no labels)", {
  dataPath <- system.file("extdata", "test_file_2.xlsx", package = "ProtStatsWF")
  sampleInfoPath <- system.file("extdata", "test_file_1_sampleInfo.xlsx", package = "ProtStatsWF")

  pData <- prepareDataSE(dataPath = dataPath, intensityColumns = 3:11,
                         proteinNameColumn = "peptide",
                         normMethod = "median",
                         sampleInfoPath = sampleInfoPath,
                         verbose = FALSE)
  pResult <- PCA_Plot(SE = pData$SE, imputationMethod = "none",
                      verbose = FALSE)

  vdiffr::expect_doppelganger("PCA_test_file_1", pResult$plot)
  expect_snapshot(pResult$D_PCA_plot)
  expect_snapshot(pResult$pca)
})


test_that("PCA plot for test_file_2 (no imputation, shape, colour, labels)", {
  dataPath <- system.file("extdata", "test_file_2.xlsx", package = "ProtStatsWF")
  sampleInfoPath <- system.file("extdata", "test_file_1_sampleInfo.xlsx",
                                package = "ProtStatsWF")

  pData <- prepareDataSE(dataPath = dataPath, intensityColumns = 3:11,
                         proteinNameColumn = "peptide",
                         normMethod = "median",
                         sampleInfoPath = sampleInfoPath,
                         verbose = FALSE)
  pResult <- PCA_Plot(SE = pData$SE, imputationMethod = "none",
                      groupForColour = "group",
                      groupForShape = "replicate",
                      label = TRUE,
                      label_seed = 123,
                      verbose = FALSE)

  vdiffr::expect_doppelganger("PCA_test_file_2", pResult$plot)
  expect_snapshot(pResult$D_PCA_plot)
  expect_snapshot(pResult$pca)
})



test_that("PCA plot for test_file_2 (imputation, shape, no colour, labels)", {
  dataPath <- system.file("extdata", "test_file_2.xlsx", package = "ProtStatsWF")
  sampleInfoPath <- system.file("extdata", "test_file_1_sampleInfo.xlsx",
                                package = "ProtStatsWF")

  pData <- prepareDataSE(dataPath = dataPath, intensityColumns = 3:11,
                         proteinNameColumn = "peptide",
                         normMethod = "median",
                         sampleInfoPath = sampleInfoPath,
                         verbose = FALSE)
  pResult <- PCA_Plot(SE = pData$SE, imputationMethod = "mean", propNA = 0.4,
                      groupForColour = NULL,
                      groupForShape = "replicate",
                      label = TRUE,
                      label_seed = 123,
                      verbose = FALSE)

  vdiffr::expect_doppelganger("PCA_test_file_3", pResult$plot)
  expect_snapshot(pResult$D_PCA_plot)
  expect_snapshot(pResult$pca)
})


test_that("PCA plot for test_file_2 (imputation, no shape, colour, no labels, specified group colours)", {
  dataPath <- system.file("extdata", "test_file_2.xlsx", package = "ProtStatsWF")
  sampleInfoPath <- system.file("extdata", "test_file_1_sampleInfo.xlsx",
                                package = "ProtStatsWF")

  pData <- prepareDataSE(dataPath = dataPath, intensityColumns = 3:11,
                         proteinNameColumn = "peptide",
                         normMethod = "median",
                         sampleInfoPath = sampleInfoPath,
                         verbose = FALSE)
  pResult <- PCA_Plot(SE = pData$SE, imputationMethod = "mean", propNA = 0.4,
                      groupForColour = "group",
                      groupForShape = NULL,
                      label = FALSE,
                      labelSeed = 123,
                      groupColours =  c("yellow", "red", "blue"),
                      verbose = FALSE)

  vdiffr::expect_doppelganger("PCA_test_file_4", pResult$plot)
  expect_snapshot(pResult$D_PCA_plot)
  expect_snapshot(pResult$pca)
})



