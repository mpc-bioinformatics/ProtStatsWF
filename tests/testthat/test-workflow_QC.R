test_that("Test MA plots", {

  # Create a temporary directory so no permanent files are put on a package users directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  workflow_QC(D_hcc,
              groupColumn = "Group",
              group2Column = NULL,
              groupColours = NULL,

              outPath = temp_dir,
              outType = "xlsx",
              suffix = "",
              NAOut = "NA",
              verbose = FALSE,

              baseSize = 15,
              plotDevice = "pdf",
              plotHeight_BP_VV = 10,
              plotWidth_BP_VV = 15,
              plotHeight_PCA_MA = 15,
              plotWidth_PCA_MA = 15,
              plotDPI = 300,

              boxplotMethod = "boxplot",

              MAMaxPlots = 5,
              MAAlpha = 1,

              PCAImputeMethod = "mean",
              PCAPropNA = 0,
              PCAScale = TRUE,
              PCAAlpha = 1,
              PCALabel = FALSE,
              PCALabelSeed = NA,
              PCALabelSize = 4,
              PCAXlim = NULL,
              PCAYlim = NULL,
              PCAPointSize = 4
  )

  expect_true(file.exists(file.path(temp_dir, "boxplot.pdf")))
  expect_gt(file.info(file.path(temp_dir, "boxplot.pdf"))$size, 0) # test that file size > 0
  expect_true(file.exists(file.path(temp_dir, "D_PCA.csv")))
  expect_gt(file.info(file.path(temp_dir, "D_PCA.csv"))$size, 0) # test that file size > 0
  expect_true(file.exists(file.path(temp_dir, "D_validvalues.csv")))
  expect_gt(file.info(file.path(temp_dir, "D_validvalues.csv"))$size, 0) # test that file size > 0
  expect_true(file.exists(file.path(temp_dir, "MA_Plots.pdf")))
  expect_gt(file.info(file.path(temp_dir, "MA_Plots.pdf"))$size, 0) # test that file size > 0
  expect_true(file.exists(file.path(temp_dir, "PCA_data_after_imputation.csv")))
  expect_gt(file.info(file.path(temp_dir, "PCA_data_after_imputation.csv"))$size, 0) # test that file size > 0
  expect_true(file.exists(file.path(temp_dir, "PCA_plot.pdf")))
  expect_gt(file.info(file.path(temp_dir, "PCA_plot.pdf"))$size, 0) # test that file size > 0
  expect_true(file.exists(file.path(temp_dir, "valid_value_plot.pdf")))
  expect_gt(file.info(file.path(temp_dir, "valid_value_plot.pdf"))$size, 0) # test that file size > 0

})

