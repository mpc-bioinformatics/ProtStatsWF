
test_that("Prepare data from csv file", {
  file_proteins <- system.file("extdata", "proteins_HCC.csv", package = "ProtStatsWF")
  file_clinical <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")
  D <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
                     proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
                     sampleNameColumn = "Sample", fileType = "csv", verbose = FALSE)

  vvplot <- ValidValuePlot(D$D_long,
                           groupColumn = "Group",
                           groupColours = NULL,
                           baseSize = 15)

  expect_snapshot(vvplot$table)
  vdiffr::expect_doppelganger("ValidValuePlot_test_file_1", vvplot$plot)
})


## TODO: test different colours and base_size
