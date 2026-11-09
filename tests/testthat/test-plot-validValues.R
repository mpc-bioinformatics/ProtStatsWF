
test_that("Prepare data from csv file", {
  vvplot <- ValidValuePlot(D_hcc$D_long,
                           groupColumn = "Group",
                           groupColours = NULL,
                           baseSize = 15)

  expect_snapshot(vvplot$table)
  vdiffr::expect_doppelganger("ValidValuePlot_test_file_1", vvplot$plot)
})

test_that("ValidValuePlot works without grouping", {
  vvplot <- ValidValuePlot(D_hcc$D_long,
                           groupColumn = NULL,
                           groupColours = NULL,
                           baseSize = 15)

  expect_s3_class(vvplot$plot, "ggplot")
  expect_true(all(c(".sample", "nrvalid", "meanvalid") %in% colnames(vvplot$table)))
  expect_false("group" %in% colnames(vvplot$table))
})


## TODO: test different colours and base_size
