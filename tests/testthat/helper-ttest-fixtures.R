# Shared fixtures for ttest-related tests.
# testthat sources helper-*.R files automatically before running any test file.

file_proteins <- system.file("extdata", "proteins_HCC.csv", package = "ProtStatsWF")
file_clinical  <- system.file("extdata", "clinical_data.csv", package = "ProtStatsWF")

D_hcc <- prepareDataSE(dataPath = file_proteins, intensityColumns = 6:43,
                       proteinNameColumn = "Protein", sampleInfoPath = file_clinical,
                       sampleNameColumn = "Sample", fileType = "csv", verbose = FALSE)

DATA_hcc  <- as.data.frame(SummarizedExperiment::assay(D_hcc$SE, "intensity_norm"))
ID_hcc    <- as.data.frame(SummarizedExperiment::rowData(D_hcc$SE))
group_hcc <- droplevels(factor(SummarizedExperiment::colData(D_hcc$SE)[, "Group"]))

ttest_res <- ttest(D = DATA_hcc, id = ID_hcc, group = group_hcc,
                   sample = SummarizedExperiment::colData(D_hcc$SE)[, "PatientID"],
                   logBeforeTest = FALSE, delogForFC = TRUE, logBase = 2,
                   minObsPerGroup = 3, paired = TRUE)

fc_col <- paste0("FC_", levels(group_hcc)[[1]], "_divided_by_", levels(group_hcc)[[2]])

sig_cats   <- calculate_significance_categories_ttest(
  p = ttest_res$p, pAdj = ttest_res$p.fdr, fc = ttest_res[[fc_col]]
)
candidates <- which(as.character(sig_cats) == "significant after FDR correction")

