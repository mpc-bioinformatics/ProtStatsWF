test_that("Calculate candidate boxplots for a ttest ", {
  
  pData <- openxlsx::read.xlsx(xlsxFile = test_path("testdata", "result_ttest.xlsx"), na.strings = c("NA", "NaN", "Filtered","#NV"))
  
  candidates <- c(2:3)
  
  data <- list("D" = pData[,3:8], "ID" = pData[,1:2], 
               "group" = factor(c("state1", "state1", "state1", "state2", "state2", "state2")))

  set.seed(42)
  
  Boxplots_candidates(D = data[["D"]][candidates, ], 
                      protein.names = data[["ID"]][candidates, "protein"],
                      group = data[["group"]],
                      plot_device = "png",
                      output_path = paste0(test_path("testdata")) )
  
  

  expect_snapshot_file(path = test_path("testdata", "boxplots_candidates_P02671.png"), name = "candidate_boxplots_1" )
  expect_snapshot_file(path = test_path("testdata", "boxplots_candidates_G3UYD0_G3UYJ6_Q3UHU8_Q9ESZ8.png"), name = "candidate_boxplots_2" )

  })



test_that("Calculate candidate boxplots for an ANOVA ", {
  
  pData <- openxlsx::read.xlsx(xlsxFile = test_path("testdata", "result_ANOVA.xlsx"), na.strings = c("NA", "NaN", "Filtered","#NV"))
  
  p_posthoc_columns <- grep("p.posthoc.", colnames(pData))
  fc_columns <- grep("FC_", colnames(pData))
  
  # filter FC columns because FC is calculated "twice" e.g. for state1_divided_state2 and state1_divided_state2
  if(fc_columns[[1]]%%2 == 0){
    fc_columns <- fc_columns[fc_columns %% 2 == 0]
  }else{
    fc_columns <- fc_columns[fc_columns %% 2 != 0]
  }
  
  candidates <- list()
  
  for (i in 1:length(p_posthoc_columns)) {
    significance <- factor(calculate_significance_categories_ANOVA(p_anova = pData[["p.anova"]],
                                                                   p_anova_adj = pData[["p.anova.fdr"]],
                                                                   p_posthoc = pData[[p_posthoc_columns[[i]]]],
                                                                   fc = pData[[fc_columns[[i]]]]))
    candidates[[i]] <- as.character(significance)
    candidates[[i]] <- which(candidates[[i]] == "significant after FDR correction")
  }
  
  candidates <- unique(unlist(candidates))
  candidates <- 11:12

  data <- list("D" = pData[,3:11], "ID" = pData[,1:2], 
               "group" = factor(c("state1", "state1", "state1", "state2", "state2", "state2", "state3", "state3", "state3")))
  
  set.seed(42)
  
  Boxplots_candidates(D = data[["D"]][candidates, ], 
                      protein.names = data[["ID"]][candidates, "protein"],
                      group = data[["group"]],
                      plot_device = "png",
                      output_path = paste0(test_path("testdata")) )
  
  
  
  expect_snapshot_file(path = test_path("testdata", "boxplots_candidates_Q61586.png"), name = "candidate_boxplots_3" )
  expect_snapshot_file(path = test_path("testdata", "boxplots_candidates_A0A087WPR7_A0A087WSP0_E9Q9X1_Q91ZU6_S4R1P5.png"), name = "candidate_boxplots_4" )
  
})