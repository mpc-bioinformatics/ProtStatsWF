# Extracted from test-plot-heatmap.R:8

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "ProtStatsWF", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(14)
hm <- Heatmap_with_groups(D = DATA_hcc[candidates, ],
                             id = ID_hcc[candidates, ],
                             groups = group_hcc, verbose = FALSE)
vdiffr::expect_doppelganger("Heatmap_with_groups_default", hm[["heatmap"]])
