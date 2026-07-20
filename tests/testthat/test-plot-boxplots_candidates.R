test_that("Boxplots_candidates creates output PDF", {
  #skip_if(length(candidates) == 0, "No significant candidates in test data")

  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  Boxplots_candidates(D = DATA_hcc[candidates[1:5], ],
                      proteinNames = ID_hcc[candidates[1:5], "Protein"],
                      group = group_hcc,
                      outputPath = temp_dir)

  pdf_path <- file.path(temp_dir, "boxplots_candidates.pdf")
  expect_true(file.exists(pdf_path))
  expect_gt(file.info(pdf_path)$size, 0) # test that file size > 0
})


#### TODO: somehow test if the plots are changing?
