# test_that("Regular ValidValue Plot", {
#   pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = TRUE)
#   pResult <- ValidValuePlot(D_long = pData[["D_long"]])
#   expect_snapshot(pResult[["table"]])
#   vdiffr::expect_doppelganger("ValidValuePlot", pResult[["plot"]])
# })



# NEW #

test_that("Regular ValidValue Plot for test_file_1", {
  pData <- prepareData(data_path = test_path("testdata", "test_file_1.xlsx"), intensity_columns = 3:11, use_groups = TRUE)
  pResult <- ValidValuePlot(D_long = pData[["D_long"]])
  expect_snapshot(pResult[["table"]])
  vdiffr::expect_doppelganger("ValidValuePlot_test_file_1", pResult[["plot"]])
})