test_that("Boxplots_candidates creates output PDF", {

  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  BoxplotsCandidates(SE = D_hcc$SE[candidates[1:5], ],
                      assay = "intensity_norm",
                      groupColumn = "Group",
                      proteinNameColumn = "Protein",
                      outputPath = temp_dir)


  pdf_path <- file.path(temp_dir, "boxplots_candidates.pdf")
  expect_true(file.exists(pdf_path))
  expect_gt(file.info(pdf_path)$size, 0) # test that file size > 0
})


#### TODO: somehow test if the plots are changing?
