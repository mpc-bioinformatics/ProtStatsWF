test_that("Heatmap_with_groups matches snapshot", {

  set.seed(14) # because ComplexHeatmap uses randomly chosen colours
  hm <- Heatmap_with_groups(D = DATA_hcc[candidates, ],
                             id = ID_hcc[candidates, ],
                             groups = data.frame(group = group_hcc), verbose = FALSE)

  vdiffr::expect_doppelganger("Heatmap_with_groups_default", hm[["heatmap"]])
})


