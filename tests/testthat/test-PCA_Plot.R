# test_that("PCA plot with no NAs", {
#   pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11)
#   pResult <- PCA_Plot(D = pData[["D"]])
#   expect_equal(pResult[["plot"]], NULL)
# })

# test_that("PCA plot with NAs", {
#   pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11)
#   pResult <- PCA_Plot(D = pData[["D"]], propNA = 0.2, impute = TRUE)
#   vdiffr::expect_doppelganger("PCA plot", pResult[["plot"]])
# })

# test_that("PCA plot with NAs and label", {
#   pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11)
#   pResult <- PCA_Plot(D = pData[["D"]], propNA = 0.2, impute = TRUE, seed = 42)
#   vdiffr::expect_doppelganger("PCA plot labeled", pResult[["plot"]])
# })



# NEW #

test_that("PCA plot with no NAs for test_file_1", {
  pData <- prepareData(data_path = test_path("testdata", "test_file_1.xlsx"), intensity_columns = 3:11)
  pResult <- PCA_Plot(D = pData[["D"]])
  expect_equal(pResult[["plot"]], NULL)
})

test_that("PCA plot with NAs for test_file_1", {
  pData <- prepareData(data_path = test_path("testdata", "test_file_1.xlsx"), intensity_columns = 3:11)
  pResult <- PCA_Plot(D = pData[["D"]], propNA = 0.2, impute = TRUE)
  vdiffr::expect_doppelganger("PCA_plot_test_file_1", pResult[["plot"]])
})

test_that("PCA plot with NAs and label for test_file_1", {
  pData <- prepareData(data_path = test_path("testdata", "test_file_1.xlsx"), intensity_columns = 3:11)
  pResult <- PCA_Plot(D = pData[["D"]], propNA = 0.2, impute = TRUE, seed = 42)
  vdiffr::expect_doppelganger("PCA_plot_labeled_test_file_1", pResult[["plot"]])
})