test_that("Calculate candidate boxplots for a ttest ", {
  
  # Skip this test on continuous integration systems like GitHub Actions
  # The function expect_snapshot_file is otherwise too strict 
  # and there is no way to get a few pixel of tolerance
  testthat::skip_on_ci()
  
  # Create temporary directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE)) 
  
  pData <- openxlsx::read.xlsx(xlsxFile = system.file("extdata", "result_ttest.xlsx", package = "ProtStatsWF"), na.strings = c("NA", "NaN", "Filtered","#NV"))
  
  candidates <- c(2:3)
  
  data <- list("D" = pData[,3:8], "ID" = pData[,1:2], 
               "group" = factor(c("state1", "state1", "state1", "state2", "state2", "state2")))

  set.seed(42)
  
  Boxplots_candidates(D = data[["D"]][candidates, ], 
                      protein.names = data[["ID"]][candidates, "protein"],
                      group = data[["group"]],
                      plot_device = "png",
                      output_path = temp_dir )
  
  

  expect_snapshot_file(path = file.path(temp_dir, "boxplots_candidates_P02671.png"), name = "candidate_boxplots_1" )
  expect_snapshot_file(path = file.path(temp_dir, "boxplots_candidates_G3UYD0_G3UYJ6_Q3UHU8.png"), name = "candidate_boxplots_2" )

  })



test_that("Calculate candidate boxplots for an ANOVA ", {
  
  # Skip this test on continuous integration systems like GitHub Actions
  # The function expect_snapshot_file is otherwise too strict 
  # and there is no way to get a few pixel of tolerance
  testthat::skip_on_ci()
  
  # Create temporary directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE)) 
  
  pData <- openxlsx::read.xlsx(xlsxFile = system.file("extdata", "result_ANOVA.xlsx", package = "ProtStatsWF"), na.strings = c("NA", "NaN", "Filtered","#NV"))
  
  candidates <- 11:12

  data <- list("D" = pData[,3:11], "ID" = pData[,1:2], 
               "group" = factor(c("state1", "state1", "state1", "state2", "state2", "state2", "state3", "state3", "state3")))
  
  set.seed(42)
  
  Boxplots_candidates(D = data[["D"]][candidates, ], 
                      protein.names = data[["ID"]][candidates, "protein"],
                      group = data[["group"]],
                      plot_device = "png",
                      output_path = temp_dir )
  
  
  
  expect_snapshot_file(path = file.path(temp_dir, "boxplots_candidates_Q61586.png"), name = "candidate_boxplots_3" )
  expect_snapshot_file(path = file.path(temp_dir, "boxplots_candidates_A0A087WPR7_A0A087WSP0_E9Q9X1.png"), name = "candidate_boxplots_4" )
  
})