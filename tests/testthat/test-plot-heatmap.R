test_that("default heatmap matches snapshot", {

  set.seed(14) # because ComplexHeatmap uses randomly chosen colours
  hm <- heatmap(D = DATA_hcc[candidates, ],
                             id = ID_hcc[candidates, ],
                             groups = data.frame(group = group_hcc), verbose = FALSE)

  vdiffr::expect_doppelganger("default_heatmap", hm[["heatmap"]])
})


