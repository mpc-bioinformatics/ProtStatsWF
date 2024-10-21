test_that("Calculate ttest ", {
  
  data <- prepareTtestData(data_path = test_path("testdata", "test_file_2.xlsx"), intensity_columns = 3:8)
  
  data[["ID"]] = data[["ID"]][,1:2]
  
  expect_snapshot(data[["D"]])
  expect_snapshot(data[["ID"]])
  
  pData <- ttest(D = data[["D"]], id = data[["ID"]], 
                 group = data[["group"]],  sample = data[["sample"]], 
                 paired = FALSE, var.equal = FALSE,
                 log_before_test = TRUE, delog_for_FC = TRUE, log_base = 2,
                 min_obs_per_group = 3, min_obs_per_group_ratio = NULL,
                 filename = paste0(test_path("testdata", "result_ttest.xlsx")))
  
  expect_snapshot(pData)
  
})