
test_that("Prepare data from csv file", {
  vvplot <- ValidValuePlot(D_hcc$D_long,
                           groupColumn = "Group",
                           groupColours = NULL,
                           baseSize = 15)

  expect_snapshot(vvplot$table)
  vdiffr::expect_doppelganger("ValidValuePlot_test_file_1", vvplot$plot)
})


## TODO: test different colours and base_size
