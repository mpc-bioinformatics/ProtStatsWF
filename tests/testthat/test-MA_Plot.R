test_that("Single MA plot", {

  pData <- prepareData(data_path = test_path("testdata", "test_file_MA_plots.xlsx"), intensity_columns = 3:6)
  D = pData[["D"]]
  s1 = D[,1]
  s2 = D[,2]

  grDevices::png(file = test_path("testdata", "MA_plot_snapshot.png"))
  MA_Plot_single(sample_1 = s1, sample_2 = s2, do_log_transformation = FALSE, alpha = FALSE)
  grDevices::dev.off()

  expect_snapshot_file(path = test_path("testdata", "MA_plot_snapshot.png"), name = "MA_Plot")
})

test_that("plot", {

  pData <- prepareData(data_path = test_path("testdata", "test_file_MA_plots.xlsx"), intensity_columns = 3:6)
  D = pData[["D"]]

  pResult <- MA_Plots(D = D, do_log_transformation = FALSE, output_path = test_path("testdata"))

  expect_snapshot(pResult)
})

