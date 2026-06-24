
test_that("Test MA plots", {
  file_proteins <- system.file("extdata", "proteins_HCC.csv", package = "ProtStatsWF")
  file_clinical <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")
  D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
                     proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
                     sampleNameColumn = "Sample", fileType = "csv", verbose = FALSE)


  # Create a temporary directory so no permanent files are put on a package users directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  MA_Plots(D = SummarizedExperiment::assay(D$SE),
           outPath = temp_dir, suffix = "",
           maxPlots = 10, alpha = 1,
           plotHeight = 15,
           plotWidth = 15,
           sampling = 1, verbose = FALSE)

  expect_true(file.exists(file.path(temp_dir, "MA_Plots.pdf")))

})

## TODO: test different other settings
## TODO: test single MA plot



