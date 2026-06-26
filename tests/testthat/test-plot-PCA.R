
test_that("Test PCA plot", {
  file_proteins <- system.file("extdata", "proteins_HCC.csv", package = "ProtStatsWF")
  file_clinical <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")
  D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
                     proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
                     sampleNameColumn = "Sample", fileType = "csv", verbose = FALSE)

  PCA <- PCA_Plot(D$SE,
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
  #expect_snapshot(PCA$pca, variant = os_name)
  expect_snapshot(PCA$filtered_data)
  #expect_snapshot(PCA$loadings)

  vdiffr::expect_doppelganger("PCA_test_file_1", PCA$plot)
})

## TODO: test different colours and base_size and other settings
