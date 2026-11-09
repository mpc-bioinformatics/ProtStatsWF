
test_that("Test standard PCA plot", {

  PCA <- PCA_Plot(D_hcc$SE,
                  groupForColour = "Group",
                  colourType = "discrete",
                  groupForShape = "Gender",
                  assay = "intensity_norm",

                  imputeMethod = "mean",
                  propNA = 0,
                  scale = TRUE,
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
  #expect_snapshot(PCA$pca, variant = os_name)
  expect_snapshot(PCA$filtered_data)
  #expect_snapshot(PCA$loadings)

  vdiffr::expect_doppelganger("PCA_test_file_1", PCA$plot)
})



test_that("Test PCA plot with only shape", {

  PCA <- PCA_Plot(D_hcc$SE,
                  groupForColour = NULL,
                  colourType = "discrete",
                  groupForShape = "Gender",
                  assay = "intensity_norm",

                  imputeMethod = "mean",
                  propNA = 0,
                  scale = TRUE,
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
  #expect_snapshot(PCA$pca, variant = os_name)
  expect_snapshot(PCA$filtered_data)
  #expect_snapshot(PCA$loadings)

  vdiffr::expect_doppelganger("PCA_test_file_2", PCA$plot)
})



test_that("Test PCA plot with only colour and continuous colour", {

  PCA <- PCA_Plot(D_hcc$SE,
                  groupForColour = "Age",
                  colourType = "continuous",
                  groupForShape = NULL,
                  assay = "intensity_norm",

                  imputeMethod = "mean",
                  propNA = 0,
                  scale = TRUE,
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
  #expect_snapshot(PCA$pca, variant = os_name)
  expect_snapshot(PCA$filtered_data)
  #expect_snapshot(PCA$loadings)

  vdiffr::expect_doppelganger("PCA_test_file_3", PCA$plot)
})


## TODO: test different colours and base_size and other settings
