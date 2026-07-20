test_that("pvalue_foldchange_histogram matches snapshot", {
  hists <- pvalue_foldchange_histogram(RES = ttest_res,
                                       columnNameP = "p", columnNamePadj = "p.fdr",
                                       columnNameFC = fc_col)
  vdiffr::expect_doppelganger("histogram_p_value", hists$histogram_p_value)
  vdiffr::expect_doppelganger("histogram_adjusted_p_value", hists$histogram_adjusted_p_value)
  vdiffr::expect_doppelganger("histogram_fold_change", hists$histogram_fold_change)
})
