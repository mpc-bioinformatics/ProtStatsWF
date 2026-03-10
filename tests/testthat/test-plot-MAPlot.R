test_that("Single MA plot", {

  # Skip this test on continuous integration systems like GitHub Actions
  # The function expect_snapshot_file is otherwise too strict
  # and there is no way to get a few pixel of tolerance
  testthat::skip_on_ci()

  # Create a temporary directory so no permanent files are put on a package users directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  png_file_path <- file.path(temp_dir, "result_MA_plot_snapshot.png")

  dataPath <- system.file("extdata", "test_file_MA_plots.xlsx", package = "ProtStatsWF")
  pData <- prepareDataSE(dataPath = dataPath, intensityColumns = 3:6, 
                         proteinNameColumn = "peptides", verbose = FALSE)
  D = SummarizedExperiment::assays(pData$SE)$intensity_norm
  s1 = D[,1]
  s2 = D[,2]

  grDevices::png(file = png_file_path)
  MA_Plot_single(sample_1 = s1, sample_2 = s2, alpha = FALSE)
  grDevices::dev.off()

  expect_snapshot_file(path = png_file_path, name = "MA_Plot.png", 
                       variant = Sys.info()[["sysname"]])
})



test_that("MA-plots", {

  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  #on.exit(unlink(temp_dir, recursive = TRUE))

  dataPath <- system.file("extdata", "test_file_MA_plots.xlsx", package = "ProtStatsWF")
  pData <- prepareDataSE(dataPath = dataPath, intensityColumns = 3:6, 
                         proteinNameColumn = "peptides", verbose = FALSE)
  D = SummarizedExperiment::assays(pData$SE)$intensity_norm

  expect_message({
    MA_Plots(D = D, output_path = NULL, suffix = "_result")
  }, "6 MA plots generated.")
    
  #expect_snapshot_file(path = paste0(temp_dir, "\\MA_Plots_result.pdf"), 
  #                     name = "MA_Plot_result.pdf", 
  #                     variant = Sys.info()[["sysname"]])
})

