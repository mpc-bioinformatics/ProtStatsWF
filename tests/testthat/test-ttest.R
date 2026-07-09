test_that("ttest returns consistent dimensions and column names", {
  expect_snapshot(dim(ttest_res))
  expect_snapshot(colnames(ttest_res))
})

test_that("ttest statistical results are consistent", {
  expect_snapshot(ttest_res[1:5, c("p", "p.fdr", fc_col)])
})
