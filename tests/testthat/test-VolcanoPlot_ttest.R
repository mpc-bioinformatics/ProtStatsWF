test_that("VolcanoPlot_ttest matches snapshot", {
  vp <- VolcanoPlot_ttest(RES = ttest_res,
                          columnNameP = "p", columnNamePadj = "p.fdr",
                          columnNameFC = fc_col)
  vdiffr::expect_doppelganger("VolcanoPlot_ttest_default", vp)
})
