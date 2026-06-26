
test_that("Test Boxplots", {
  file_proteins <- system.file("extdata", "proteins_HCC.csv", package = "ProtStatsWF")
  file_clinical <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")
  D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
                     proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
                     sampleNameColumn = "Sample", fileType = "csv", verbose = FALSE)

  boxplots <- Boxplots(D$D_long,
                       method = "boxplot",
                       groupColumn = "Group",
                       groupColours = NULL,
                       baseSize = 15,
                       lwd = 0.5,
                       outlierSize = 1)
  vdiffr::expect_doppelganger("Boxplots_test_file_1", boxplots)


  boxplots2 <- Boxplots(D$D_long,
                       method = "boxplot",
                       groupColumn = "Group",
                       groupColours = c("green", "yellow"),
                       baseSize = 25,
                       lwd = 2,
                       outlierSize = 0.2)
  vdiffr::expect_doppelganger("Boxplots_test_file_2", boxplots2)

  violin <- Boxplots(D$D_long,
                       method = "violinplot",
                       groupColumn = "Group",
                       groupColours = NULL,
                       baseSize = 15,
                       lwd = 0.5,
                       outlierSize = 1)
  vdiffr::expect_doppelganger("Boxplots_test_file_3", violin)

})

## TODO: test different colours and base_size and other settings
