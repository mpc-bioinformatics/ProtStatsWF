# Extracted from test-plot-pvalue_FC_histogram.R:5

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "ProtStatsWF", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
hists <- pvalue_foldchange_histogram(RES = ttest_res,
                                       columnNameP = "p", columnNamePadj = "p.fdr",
                                       columnNameFC = fc_col)
vdiffr::expect_doppelganger("histogram_p_value", hists$histogram_p_value)
