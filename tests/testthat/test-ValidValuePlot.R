test_that("the spreadsheet is correct", {
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = TRUE)
  pGraph <- ValidValuePlot(D_long = pData[[6]])
  expect_snapshot(pGraph[[1]])
})

test_that("the plot is correct", {
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = TRUE)
  pGraph <- ValidValuePlot(D_long = pData[[6]])
  vdiffr::expect_doppelganger("ValidValuePlot", pGraph[2])
})
