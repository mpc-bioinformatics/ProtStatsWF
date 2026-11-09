# Extracted from test-plot-volcano.R:5

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "ProtStatsWF", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
vp <- VolcanoPlot_ttest(RES = ttest_res,
                          columnNameP = "p", columnNamePadj = "p.fdr",
                          columnNameFC = fc_col)
vdiffr::expect_doppelganger("VolcanoPlot_ttest_default", vp)
