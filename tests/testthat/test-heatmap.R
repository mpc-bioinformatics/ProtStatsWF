test_that("Calculate heatmap for a ttest ", {

  pData <- openxlsx::read.xlsx(xlsxFile = system.file("extdata", "result_ttest.xlsx", package = "ProtStatsWF"), na.strings = c("NA", "NaN", "Filtered","#NV"))

  significance <- calculate_significance_categories_ttest(p = pData[["p"]],
                                                          p_adj = pData[["p.fdr"]],
                                                          fc = pData[["FC_state1_divided_by_state2"]])
  candidates <- as.character(significance)
  candidates <- which(candidates == "significant after FDR correction")

  data <- list("D" = pData[,3:8], "ID" = pData[,1:2],
               "group" = factor(c("state1", "state1", "state1", "state2", "state2", "state2")))

  set.seed(14)

  t_heatmap <- Heatmap_with_groups(D = data[["D"]][candidates, ],
                                   id = data[["ID"]][candidates, ],
                                   groups = data[["group"]])

  vdiffr::expect_doppelganger("Heatmap_ttest", t_heatmap[["heatmap"]])

})



test_that("Calculate heatmap for an ANOVA ", {

  pData <- openxlsx::read.xlsx(xlsxFile = system.file("extdata", "result_ANOVA.xlsx", package = "ProtStatsWF"), na.strings = c("NA", "NaN", "Filtered","#NV"))

  #candidates <- 11:12

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

  union_candidates <- unique(unlist(candidates))

  data <- list("D" = pData[,3:11], "ID" = pData[,1:2],
               "group" = factor(c("state1", "state1", "state1", "state2", "state2", "state2", "state3", "state3", "state3")))

  set.seed(14)

  t_heatmap <- Heatmap_with_groups(D = data[["D"]][union_candidates, ],
                                   id = data[["ID"]][union_candidates, ],
                                   groups = data[["group"]])

  vdiffr::expect_doppelganger("Heatmap_ANOVA", t_heatmap[["heatmap"]])

})
