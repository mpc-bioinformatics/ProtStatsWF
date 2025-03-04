test_that("Calculate ttest ", {
  
  # Create temporary directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  data <- prepareTtestData(data_path = system.file("extdata", "test_file_2.xlsx", package = "ProtStatsWF"), intensity_columns = 3:8)
  
  data[["ID"]] = data[["ID"]][,1:2]
  
  expect_snapshot(data[["D"]])
  expect_snapshot(data[["ID"]])
  
  pData <- ttest(D = data[["D"]], id = data[["ID"]], 
                 group = data[["group"]],  sample = data[["sample"]], 
                 paired = FALSE, var.equal = FALSE,
                 log_before_test = TRUE, delog_for_FC = TRUE, log_base = 2,
                 min_obs_per_group = 3, min_obs_per_group_ratio = NULL,
                 filename = paste0(file.path(temp_dir, "result_ttest.xlsx")))
  
  expect_snapshot(pData)
  
})