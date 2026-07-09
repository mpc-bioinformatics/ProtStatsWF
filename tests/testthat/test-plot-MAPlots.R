
test_that("Test MA plots", {
  # Create a temporary directory so no permanent files are put on a package users directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  MA_Plots(D = SummarizedExperiment::assay(D_hcc$SE),
           outPath = temp_dir, suffix = "",
           maxPlots = 2, alpha = 1,
           plotHeight = 15,
           plotWidth = 15,
           sampling = 1, verbose = FALSE)

  expect_true(file.exists(file.path(temp_dir, "MA_Plots.pdf")))

})

## TODO: test different other settings
## TODO: test single MA plot with snapshot



