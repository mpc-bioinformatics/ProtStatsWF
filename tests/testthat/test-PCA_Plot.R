test_that("PCA plot with no NAs", {
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11)
  pResult <- PCA_Plot(D = pData[["D"]])
  expect_equal(pResult[["plot"]], NULL)
})

test_that("PCA plot with NAs", {
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11)
  pResult <- PCA_Plot(D = pData[["D"]], propNA = 0.2, impute = TRUE)
  vdiffr::expect_doppelganger("PCA plot", pResult[["plot"]])
})

test_that("PCA plot with NAs and label", {
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11)
  pResult <- PCA_Plot(D = pData[["D"]], propNA = 0.2, impute = TRUE, seed = 42)
  vdiffr::expect_doppelganger("PCA plot labeled", pResult[["plot"]])
})
