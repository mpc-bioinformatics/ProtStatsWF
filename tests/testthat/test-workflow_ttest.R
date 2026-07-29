
test_that("workflow_ttest runs without error on HCC data", {
  file_proteins <- system.file("extdata", "proteins_HCC_small.csv", package = "ProtStatsWF")
  file_clinical  <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")

  D_hcc_small <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
                         proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
                         sampleNameColumn = "Sample", verbose = FALSE)

  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))


  result <- workflow_ttest(D = D_hcc_small, groupColumn = "Group", outputPath = temp_dir,
                           verbose = FALSE, thresFC = 1, thresP = 0.6,
                           significantAfterFDR = FALSE)

  #
  # result <- workflow_ttest(D = D_hcc_small, groupColumn = "Group", outputPath = "tmp",
  #                          verbose = TRUE, thresFC = 1, thresP = 0.6,
  #                          significantAfterFDR = FALSE)


  # Returns a message log
  expect_type(result, "list")
  expect_named(result, c("test_results", "significance"))


  # Core output files are created
  expect_true(file.exists(file.path(temp_dir, "results_ttest.xlsx")))
  expect_gt(file.info(file.path(temp_dir, "results_ttest.xlsx"))$size, 0) # test that file size > 0
  expect_true(file.exists(file.path(temp_dir, "results_onoff.xlsx")))
  expect_gt(file.info(file.path(temp_dir, "results_onoff.xlsx"))$size, 0) # test that file size > 0
  expect_true(file.exists(file.path(temp_dir, "volcano_plot.pdf")))
  expect_gt(file.info(file.path(temp_dir, "volcano_plot.pdf"))$size, 0) # test that file size > 0
  expect_true(file.exists(file.path(temp_dir, "histogram_p_value.pdf")))
  expect_gt(file.info(file.path(temp_dir, "histogram_p_value.pdf"))$size, 0) # test that file size > 0
  expect_true(file.exists(file.path(temp_dir, "boxplots_candidates.pdf")))
  expect_gt(file.info(file.path(temp_dir, "boxplots_candidates.pdf"))$size, 0) # test that file size > 0
})

#### TODO: test existence of heatmap (currently only 1 candidate without missing values)


test_that("workflow_ttest errors with more than 2 groups", {
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  D <- D_hcc

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
