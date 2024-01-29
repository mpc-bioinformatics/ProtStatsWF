test_that("Single MA plot", {
  
  pData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11)
  D = pData[["D"]]
  s1 = D[,1]
  s2 = D[,2]
  
  grDevices::png(file = test_path("testdata", "MA_plot_snapshot.png"))
  MA_Plot_single(sample_1 = s1, sample_2 = s2, do_log_transformation = FALSE, alpha = FALSE)
  grDevices::dev.off()
  
  expect_snapshot_file(path = test_path("testdata", "MA_plot_snapshot.png"), name = "MA_Plot")
})
