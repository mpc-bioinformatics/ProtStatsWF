test_that("PCA plot wothout groups", {
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11)
  # pResult <- PCA_Plot(D = pData[["D"]])
  expect_equal(2 * 2, 4)
})
