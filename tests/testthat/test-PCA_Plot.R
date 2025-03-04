
test_that("PCA plot with no NAs for test_file_1", {
  pData <- prepareData(data_path = system.file("extdata", "test_file_1.xlsx", package = "ProtStatsWF"), intensity_columns = 3:11)
  pResult <- PCA_Plot(D = pData[["D"]])
  expect_equal(pResult[["plot"]], NULL)
})

test_that("PCA plot with NAs for test_file_1", {
  pData <- prepareData(data_path = system.file("extdata", "test_file_1.xlsx", package = "ProtStatsWF"), intensity_columns = 3:11)
  pResult <- PCA_Plot(D = pData[["D"]], propNA = 0.2, impute = TRUE)
  vdiffr::expect_doppelganger("PCA_plot_test_file_1", pResult[["plot"]])
})

test_that("PCA plot with NAs and label for test_file_1", {
  pData <- prepareData(data_path = system.file("extdata", "test_file_1.xlsx", package = "ProtStatsWF"), intensity_columns = 3:11)
  pResult <- PCA_Plot(D = pData[["D"]], propNA = 0.2, impute = TRUE, seed = 42)
  vdiffr::expect_doppelganger("PCA_plot_labeled_test_file_1", pResult[["plot"]])
})