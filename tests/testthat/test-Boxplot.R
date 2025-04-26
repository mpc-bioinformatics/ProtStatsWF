
test_that("Boxplot with groups for test_file_1", {
  pData <- prepareData(data_path = system.file("extdata", "test_file_1.xlsx", package = "ProtStatsWF"), intensity_columns = 3:11, use_groups = TRUE)
  pResult <- Boxplots(D_long = pData[["D_long"]], method = "boxplot")
  plot <- pResult[["plot"]]
  vdiffr::expect_doppelganger("Boxplot_test_file_1", plot)
})

test_that("Violinplot without groups for test_file_1", {
  pData <- prepareData(data_path = system.file("extdata", "test_file_1.xlsx", package = "ProtStatsWF"), intensity_columns = 3:11, use_groups = FALSE)
  pResult <- Boxplots(D_long = pData[["D_long"]], method = "violinplot")
  plot <- pResult[["plot"]]
  vdiffr::expect_doppelganger("Violinplot_test_file_1", plot)
})