test_that("Boxplot with groups", {
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = TRUE)
  pResult <- Boxplots(D_long = pData[[3]], method = "boxplot")
  plot <- pResult[[1]]
  vdiffr::expect_doppelganger("Boxplot", plot)
})

test_that("Violinplot without groups", {
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = FALSE)
  pResult <- Boxplots(D_long = pData[[3]], method = "violinplot")
  plot <- pResult[[1]]
  vdiffr::expect_doppelganger("Violinplot", plot)
})
