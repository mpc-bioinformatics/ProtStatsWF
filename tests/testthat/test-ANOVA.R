test_that("Calculate ANOVA ", {

  # Create a temporary directory so no permanent files are put on a package users directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  data <- prepareTtestData(data_path = test_path("testdata", "test_file_2.xlsx"), intensity_columns = 3:11)

  expect_snapshot(data[["D"]])
  expect_snapshot(data[["ID"]])

  pData <- ANOVA(D = data[["D"]], id = data[["ID"]],
                     group = data[["group"]],  sample = data[["sample"]],
                     paired = FALSE, var.equal = TRUE,
                     log_before_test = TRUE, delog_for_FC = TRUE, log_base = 2,
                     min_obs_per_group = 3, min_perc_per_group = NULL)

  expect_snapshot(pData)

})
