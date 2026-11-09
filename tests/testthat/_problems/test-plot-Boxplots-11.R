# Extracted from test-plot-Boxplots.R:11

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "ProtStatsWF", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
boxplots <- Boxplots(D_hcc$D_long,
                       method = "boxplot",
                       groupColumn = "Group",
                       groupColours = NULL,
                       baseSize = 15,
                       lwd = 0.5,
                       outlierSize = 1)
vdiffr::expect_doppelganger("Boxplots_test_file_1", boxplots)
