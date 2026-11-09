# Extracted from test-plot-valid_values.R:9

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "ProtStatsWF", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
vvplot <- ValidValuePlot(D_hcc$D_long,
                           groupColumn = "Group",
                           groupColours = NULL,
                           baseSize = 15)
expect_snapshot(vvplot$table)
vdiffr::expect_doppelganger("ValidValuePlot_test_file_1", vvplot$plot)
