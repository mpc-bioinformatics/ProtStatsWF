test_that("calculate_onoff returns consistent results", {
  on_off <- ProtStatsWF:::calculate_onoff(
    D    = DATA_hcc,
    id   = ID_hcc,
    group = group_hcc,
    maxValidValuesOff = 0,
    minValidValuesOn  = min(table(group_hcc)),
    proteinNamesColumn = which(colnames(ID_hcc) == "Protein")
  )

  expect_snapshot(colnames(on_off))
  expect_snapshot(nrow(on_off))
  expect_snapshot(on_off)
})
