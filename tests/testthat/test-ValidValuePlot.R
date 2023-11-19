test_that("Regular ValidValue Plot", {
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = TRUE)
  pVVP <- ValidValuePlot(D_long = pData[[3]])
  expect_snapshot(pVVP[[1]])
  vdiffr::expect_doppelganger("ValidValuePlot", pVVP[2])
})