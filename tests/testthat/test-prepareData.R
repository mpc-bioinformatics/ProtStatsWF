
test_that("Data preparation with groups (median normalization) for test_file_1", {
  prepData <- prepareData(data_path = system.file("extdata", "test_file_1.xlsx", package = "ProtStatsWF"), intensity_columns = 3:11, use_groups = TRUE, normalization = "median")
  
  expect_snapshot(prepData[["D"]])
  expect_snapshot(prepData[["ID"]])
  expect_snapshot(prepData[["D_long"]])
})

test_that("Data preparation without groups (loess normalization) for test_file_1", {
  prepData <- prepareData(data_path = system.file("extdata", "test_file_1.xlsx", package = "ProtStatsWF"), intensity_columns = 3:11, use_groups = FALSE, normalization = "loess")
  
  expect_snapshot(prepData[["D"]])
  expect_snapshot(prepData[["ID"]])
  expect_snapshot(prepData[["D_long"]])
})