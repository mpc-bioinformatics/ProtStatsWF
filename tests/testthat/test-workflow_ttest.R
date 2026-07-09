
test_that("workflow_ttest runs without error on HCC data", {
  file_proteins <- system.file("extdata", "proteins_HCC.csv", package = "ProtStatsWF")
  file_clinical <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")

  D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
                     proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
                     sampleNameColumn = "Sample", fileType = "csv",
                     verbose = FALSE)

  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))


  result <- workflow_ttest(D = D, groupColumn = "Group", outputPath = temp_dir,
                           verbose = FALSE)

  # Returns a message log
  expect_type(result, "list")
  expect_named(result, "message")
  expect_type(result$message, "character")

  # Core output files are created
  expect_true(file.exists(file.path(temp_dir, "results_ttest.xlsx")))
  expect_true(file.exists(file.path(temp_dir, "results_onoff.xlsx")))
  expect_true(file.exists(file.path(temp_dir, "volcano_plot.pdf")))
  expect_true(file.exists(file.path(temp_dir, "histogram_p_value.pdf")))
  expect_true(file.exists(file.path(temp_dir, "message_log_ttest.txt")))
})

test_that("workflow_ttest errors with more than 2 groups", {
  file_proteins <- system.file("extdata", "proteins_HCC.csv", package = "ProtStatsWF")
  file_clinical <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")

  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))


  D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
                     proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
                     sampleNameColumn = "Sample", fileType = "csv",
                     verbose = FALSE)

  # Inject a fake third group into colData to trigger the guard
  SummarizedExperiment::colData(D$SE)$FakeGroup <- factor(
    c(rep("A", 10), rep("B", 14), rep("C", 14))
  )

  expect_error(
    workflow_ttest(D = D, groupColumn = "FakeGroup", outputPath = temp_dir,
                   verbose = FALSE),
    "exactly 2 groups"
  )
})
