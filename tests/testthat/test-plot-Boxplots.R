
test_that("Test Boxplots", {

  boxplots <- Boxplots(D_long,
                       valueColumn = "intensity_norm",
                       method = "boxplot",
                       groupColumn = "Group",
                       groupColours = NULL,
                       baseSize = 15,
                       lwd = 0.5,
                       outlierSize = 1)
  vdiffr::expect_doppelganger("Boxplots_test_file_1", boxplots)


  ## test colours
  boxplots2 <- Boxplots(D_long,
                        valueColumn = "intensity_norm",
                       method = "boxplot",
                       groupColumn = "Group",
                       groupColours = c("green", "yellow"),
                       baseSize = 25,
                       lwd = 2,
                       outlierSize = 0.2)
  vdiffr::expect_doppelganger("Boxplots_test_file_2", boxplots2)

  ## test violin plot
  violin <- Boxplots(D_long,
                     valueColumn = "intensity_norm",
                       method = "violinplot",
                       groupColumn = "Group",
                       groupColours = NULL,
                       baseSize = 15,
                       lwd = 0.5,
                       outlierSize = 1)
  vdiffr::expect_doppelganger("Boxplots_test_file_3", violin)

})

test_that("Boxplots works without grouping", {
  boxplots <- Boxplots(D_long,
                       valueColumn = "intensity_norm",
                       method = "boxplot",
                       groupColumn = NULL,
                       groupColours = NULL,
                       baseSize = 15,
                       lwd = 0.5,
                       outlierSize = 1)

  vdiffr::expect_doppelganger("Boxplots_test_file_4", boxplots)
})

