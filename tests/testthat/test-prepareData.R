test_that("Data preparation with groups (median normalization)", {
  prepData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = TRUE, normalization = "median")
  
  expect_snapshot(prepData[["D"]])
  expect_snapshot(prepData[["ID"]])
  expect_snapshot(prepData[["D_long"]])
})

test_that("Data preparation without groups (loess normalization)", {
  prepData <- prepareData(data_path = test_path("testdata", "test01.xlsx"), intensity_columns = 3:11, use_groups = FALSE, normalization = "loess")
  
  expect_snapshot(prepData[["D"]])
  expect_snapshot(prepData[["ID"]])
  expect_snapshot(prepData[["D_long"]])
})
