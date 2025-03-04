test_that("Calculate on-off-heatmap for a ttest ", {
  
  pData <- openxlsx::read.xlsx(xlsxFile = system.file("extdata", "result_ttest.xlsx", package = "ProtStatsWF"), na.strings = c("NA", "NaN", "Filtered","#NV"))
  
  data <- list("D" = pData[,3:8], "ID" = pData[,1:2], 
               "group" = factor(c("state1", "state1", "state1", "state2", "state2", "state2")))
  
  min_valid_values_on <- 5 
  max_valid_values_off <- 1
  
  t_on_off_heatmap <- calculate_onoff(D = data[["D"]], 
                                      id = data[["ID"]],
                                      group = data[["group"]],
                                      max_vv_off = max_valid_values_off,
                                      min_vv_on = min_valid_values_on,
                                      protein_id_col = 1)
  
  vdiffr::expect_doppelganger("On-Off-Heatmap_ttest", t_on_off_heatmap)
  
})



test_that("Calculate on-off-heatmap for an ANOVA ", {
  
  pData <- openxlsx::read.xlsx(xlsxFile = system.file("extdata", "result_ANOVA.xlsx", package = "ProtStatsWF"), na.strings = c("NA", "NaN", "Filtered","#NV"))
  
  data <- list("D" = pData[,3:11], "ID" = pData[,1:2], 
               "group" = factor(c("state1", "state1", "state1", "state2", "state2", "state2", "state3", "state3", "state3")))
  
  min_valid_values_on <- 8 
  max_valid_values_off <- 1
  
  t_on_off_heatmap <- calculate_onoff(D = data[["D"]], 
                                      id = data[["ID"]],
                                      group = data[["group"]],
                                      max_vv_off = max_valid_values_off,
                                      min_vv_on = min_valid_values_on,
                                      protein_id_col = 1)
  
  vdiffr::expect_doppelganger("On-Off-Heatmap_ANOVA", t_on_off_heatmap)
  
})