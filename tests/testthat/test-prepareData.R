test_that("Data preparation with groups (median normalization)", {
  prepData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = TRUE, normalization = "median")
  
  expect_snapshot(prepData[[1]])
  expect_snapshot(prepData[[2]])
  expect_snapshot(prepData[[3]])
})

test_that("Data preparation without groups (loess normalization)", {
  prepData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = FALSE, normalization = "loess")
  
  expect_snapshot(prepData[[1]])
  expect_snapshot(prepData[[2]])
  expect_snapshot(prepData[[3]])
})
