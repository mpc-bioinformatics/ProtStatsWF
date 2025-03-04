test_that("Calculate histograms for a ttest ", {
  
  pData <- openxlsx::read.xlsx(xlsxFile = system.file("extdata", "result_ttest.xlsx", package = "ProtStatsWF"), na.strings = c("NA", "NaN", "Filtered","#NV"))
  
  histograms <- pvalue_foldchange_histogram(RES = pData, 
                                            columnname_p = "p", columnname_padj = "p.fdr", 
                                            columnname_FC = "FC_state1_divided_by_state2")
  
  # Warning: "Removed 5 rows containing non-finite outside the scale range (`stat_bin()`)."
  suppressWarnings(vdiffr::expect_doppelganger("histogram_p_value_ttest_file_2", histograms[1]))
  suppressWarnings(vdiffr::expect_doppelganger("histogram_FC_ttest_file_2", histograms[2]))
})



test_that("Calculate histograms for an ANOVA ", {
  
  pData <- openxlsx::read.xlsx(xlsxFile = system.file("extdata", "result_ANOVA.xlsx", package = "ProtStatsWF"), na.strings = c("NA", "NaN", "Filtered","#NV"))
  
  fc_columns <- grep("FC_", colnames(pData))
  
  # filter FC columns because FC is calculated "twice" e.g. for state1_divided_state2 and state1_divided_state2
  if(fc_columns[[1]]%%2 == 0){
    fc_columns <- fc_columns[fc_columns %% 2 == 0]
  }else{
    fc_columns <- fc_columns[fc_columns %% 2 != 0]
  }
  
  histograms <- pvalue_foldchange_histogram(RES = pData, 
                                            columnname_p = "p.anova", columnname_padj = "p.anova.fdr", 
                                            columnname_FC = colnames(pData)[[fc_columns[[1]]]])
  
  
  # Warning: "Removed 5 rows containing non-finite outside the scale range (`stat_bin()`)."
  suppressWarnings(vdiffr::expect_doppelganger("histogram_p_value_ANOVA_file_2", histograms[1]))
  suppressWarnings(vdiffr::expect_doppelganger("histogram_FC_ANOVA_file_2", histograms[2]))
})