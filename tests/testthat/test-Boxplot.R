test_that("Regular Boxplot", {
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = TRUE)
  plot <- Boxplots(D_long = pData[[6]])
  vdiffr::expect_doppelganger("Boxplot", plot)
})
