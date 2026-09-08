# Extracted from test-plot-PCA.R:38

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "ProtStatsWF", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
PCA <- PCA_Plot(D_hcc$SE,
                  groupForColour = "Group",
                  colourType = "discrete",
                  groupForShape = "Gender",
                  assay = "intensity_norm",

                  imputeMethod = "mean",
                  propNA = 0,
                  scale. = TRUE,
                  PCx = 1,
                  PCy = 2,

                  groupColours = NULL,
                  alpha = 1,
                  label = FALSE,
                  labelSeed = NA,
                  labelSize = 4,
                  xlim = NULL,
                  ylim = NULL,

                  pointSize = 4,
                  baseSize = 11,
                  NAValueColour = "grey",
                  NAValueShape = 0,
                  verbose = FALSE
  )
os_name <- tolower(Sys.info()[["sysname"]])
expect_snapshot(PCA$D_PCA_plot)
expect_snapshot(PCA$filtered_data)
vdiffr::expect_doppelganger("PCA_test_file_1", PCA$plot)
