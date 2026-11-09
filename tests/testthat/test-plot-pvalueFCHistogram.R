test_that("pvalueFCHistogram matches snapshot", {
  hists <- pvalueFCHistogram(RES = ttest_res,
                                       columnP = "p", columnPadj = "p.fdr",
                                       columnFC = fc_col)
  vdiffr::expect_doppelganger("histogram_p_value", hists$histogram_p_value)
  vdiffr::expect_doppelganger("histogram_adjusted_p_value", hists$histogram_adjusted_p_value)
  vdiffr::expect_doppelganger("histogram_fold_change", hists$histogram_fold_change)
})
