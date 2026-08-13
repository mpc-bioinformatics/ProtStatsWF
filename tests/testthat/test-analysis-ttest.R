test_that("paired ttest statistical results are consistent", {
  expect_snapshot(ttest_res)
})


test_that("unpaired ttest statistical results are consistent", {


  ttest_res <- ttest(SE = D_hcc$SE, assay = "intensity_norm",
                     groupColumn = "Group", sampleColumn = NULL,
                     logBeforeTest = FALSE, delogForFC = TRUE, logBase = 2,
                     minObs = 3, paired = FALSE,
                     verbose = FALSE)

  expect_snapshot(ttest_res)
})
