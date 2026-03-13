
test_that("Boxplot with groups for test_file_1", {
  dataPath <- system.file("extdata", "test_file_2.xlsx", package = "ProtStatsWF")
  sampleInfoPath <- system.file("extdata", "test_file_1_sampleInfo.xlsx", package = "ProtStatsWF")
  
  pData <- prepareDataSE(dataPath = dataPath, intensityColumns = 3:11, 
                         proteinNameColumn = "peptide", sampleInfoPath = sampleInfoPath, 
                         verbose = FALSE)
  plot <- Boxplots(D_long = pData$D_long, method = "boxplot", groupColumn = "group")
  vdiffr::expect_doppelganger("Boxplot_test_file_1", plot)
  
  
  ## different order of samples:
  sampleInfoPath2 <- system.file("extdata", "test_file_1_sampleInfo2.xlsx", 
                                 package = "ProtStatsWF")
  pData <- prepareDataSE(dataPath = dataPath, intensityColumns = 3:11, 
                         proteinNameColumn = "peptide", sampleInfoPath = sampleInfoPath2, 
                         verbose = FALSE)
  plot <- Boxplots(D_long = pData$D_long, method = "boxplot", groupColumn = "group")
  vdiffr::expect_doppelganger("Boxplot_test_file_1_ordered", plot)
  
  
})

test_that("Violinplot without groups for test_file_1", {
  dataPath <- system.file("extdata", "test_file_2.xlsx", package = "ProtStatsWF")
  sampleInfoPath <- system.file("extdata", "test_file_1_sampleInfo.xlsx", package = "ProtStatsWF")
  
  pData <- prepareDataSE(dataPath = dataPath, intensityColumns = 3:11, 
                         proteinNameColumn = "peptide", sampleInfoPath = sampleInfoPath, 
                         verbose = FALSE)
  plot <- Boxplots(D_long = pData$D_long, method = "violinplot", groupColumn = NULL)
  vdiffr::expect_doppelganger("Violinplot_test_file_1", plot)
})

