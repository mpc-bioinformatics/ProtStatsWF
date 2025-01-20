test_that("Single MA plot", {

  # Skip this test on continuous integration systems like GitHub Actions
  # The function expect_snapshot_file is otherwise too strict 
  # and there is no way to get a few pixel of tolerance
  # testthat::skip_on_ci()
  
  # Create a temporary directory so no permanent files are put on a package users directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE)) 
  
  png_file_path <- file.path(temp_dir, "result_MA_plot_snapshot.png")
  
  pData <- prepareData(data_path = test_path("testdata", "test_file_MA_plots.xlsx"), intensity_columns = 3:6)
  D = pData[["D"]]
  s1 = D[,1]
  s2 = D[,2]

  grDevices::png(file = png_file_path)
  MA_Plot_single(sample_1 = s1, sample_2 = s2, do_log_transformation = FALSE, alpha = FALSE)
  grDevices::dev.off()

  expect_snapshot_file(path = png_file_path, name = "MA_Plot")
})



test_that("plot", {

  pData <- prepareData(data_path = test_path("testdata", "test_file_MA_plots.xlsx"), intensity_columns = 3:6)
  D = pData[["D"]]

  pResult <- MA_Plots(D = D, do_log_transformation = FALSE, output_path = test_path("testdata"), suffix = "_result")

  expect_snapshot(pResult)
})

