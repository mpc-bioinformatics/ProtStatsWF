
test_that("Regular ValidValue Plot for test_file_2", {
  dataPath <- system.file("extdata", "test_file_2.xlsx", package = "ProtStatsWF")
  sampleInfoPath <- system.file("extdata", "test_file_1_sampleInfo.xlsx", 
                                package = "ProtStatsWF")
  
  prepData <- prepareDataSE(dataPath = dataPath, 
                            intensityColumns = 3:11,
                            proteinNameColumn = "peptide",
                            sampleInfoPath = sampleInfoPath, 
                            verbose = FALSE, 
                            normMethod = "median")
  
  pResult <- ValidValuePlot(D_long = prepData$D_long, groupColumn = "group")
  expect_snapshot(pResult$table)
  vdiffr::expect_doppelganger("ValidValuePlot_test_file_1", pResult$plot)
  
  
  # different order of samples (by replicate, not by group)
  sampleInfoPath2 <- system.file("extdata", "test_file_1_sampleInfo2.xlsx", 
                                package = "ProtStatsWF")
  prepData2 <- prepareDataSE(dataPath = dataPath, 
                            intensityColumns = 3:11,
                            proteinNameColumn = "peptide",
                            sampleInfoPath = sampleInfoPath2, 
                            verbose = FALSE, 
                            normMethod = "median")
  
  pResult2 <- ValidValuePlot(D_long = prepData2$D_long, groupColumn = "group")
  expect_snapshot(pResult2$table)
  vdiffr::expect_doppelganger("ValidValuePlot_test_file_1_ordered", pResult2$plot)
})

