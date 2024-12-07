test_that("Calculate Volcano plot for a ttest ", {
  
  pData <- openxlsx::read.xlsx(xlsxFile = test_path("testdata", "result_ttest.xlsx"), na.strings = c("NA", "NaN", "Filtered","#NV"))
  
  volcano_plot <- VolcanoPlot_ttest(RES = pData, 
                                    columnname_p = "p", columnname_padj = "p.fdr", 
                                    columnname_FC = "FC_state1_divided_by_state2")
  
  # Warning: "Removed 5 rows containing missing values or values outside the scale range (`geom_point()`)."
  # TODO: change volcano plot so missing values are not plotted
  suppressWarnings(vdiffr::expect_doppelganger("volcano_plot_ttest_file_2", volcano_plot))
})



test_that("Calculate Volcano plot for an ANOVA ", {
  
  pData <- openxlsx::read.xlsx(xlsxFile = test_path("testdata", "result_ANOVA.xlsx"), na.strings = c("NA", "NaN", "Filtered","#NV"))
  
  p_posthoc_columns <- grep("p.posthoc.", colnames(pData))
  fc_columns <- grep("FC_", colnames(pData))
  
  # filter FC columns because FC is calculated "twice" e.g. for state1_divided_state2 and state1_divided_state2
  if(fc_columns[[1]]%%2 == 0){
    fc_columns <- fc_columns[fc_columns %% 2 == 0]
  }else{
    fc_columns <- fc_columns[fc_columns %% 2 != 0]
  }
  
  volcano_plots <- VolcanoPlot_ANOVA(RES = pData, 
                                     columnname_p_ANOVA = "p.anova",
                                     columnname_p_ANOVA_adj = "p.anova.fdr",
                                     columns_FC = fc_columns,
                                     columns_p_posthoc = p_posthoc_columns )
  
  # suppressWarnings because rows with NA values are removed, which creates a warning
  # TODO: change volcano plot so missing values are not plotted
  suppressWarnings(vdiffr::expect_doppelganger("volcano_plot_ANOVA_file_2", volcano_plots[[2]]))
})