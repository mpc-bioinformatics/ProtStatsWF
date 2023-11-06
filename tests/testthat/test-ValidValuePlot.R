test_that("the spreadsheet is correct", {
  pData <- prepareData(data_path = "/Users/kalar/Documents/0_Studium/SHK-WHK/Testdata/preprocessed_peptide_data_D1.xlsx", intensity_columns = 3:17, use_groups = TRUE)
  pGraph <- ValidValuePlot(D_long = D1[[6]])
  expect_snapshot(pGraph[[1]])
})

test_that("the plot is correct", {
  pData <- prepareData(data_path = "/Users/kalar/Documents/0_Studium/SHK-WHK/Testdata/preprocessed_peptide_data_D1.xlsx", intensity_columns = 3:17, use_groups = TRUE)
  pGraph <- ValidValuePlot(D_long = D1[[6]])
  vdiffr::expect_doppelganger("ValidValuePlot", pGraph[2])
})
