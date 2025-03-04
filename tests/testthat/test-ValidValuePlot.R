
test_that("Regular ValidValue Plot for test_file_1", {
  pData <- prepareData(data_path = system.file("extdata", "test_file_1.xlsx", package = "ProtStatsWF"), intensity_columns = 3:11, use_groups = TRUE)
  pResult <- ValidValuePlot(D_long = pData[["D_long"]])
  expect_snapshot(pResult[["table"]])
  vdiffr::expect_doppelganger("ValidValuePlot_test_file_1", pResult[["plot"]])
})