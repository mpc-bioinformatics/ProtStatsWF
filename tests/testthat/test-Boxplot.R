test_that("Regular Boxplot", {
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = TRUE)
  pResult <- Boxplots(D_long = pData[[3]])
  plot <- pResult[[1]]
  vdiffr::expect_doppelganger("Boxplot", plot)
})
