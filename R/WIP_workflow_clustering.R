

workflow_clustering <- function(data_path,
                                filetype = "xlsx",
                                sep = ",",
                                dec = ".",
                                header = TRUE,
                                sheet = 1,
                                output_path,
                                na_strings = c("NA", "NaN", "Filtered","#NV"),
                                intensity_columns) {

  prepareData(data_path,
              filetype = filetype,
              sep = sep,
              dec = dec,
              header = header,
              sheet = sheet,
              intensity_columns = intensity_columns,
              na_strings = na_strings,
              zero_to_NA = TRUE,
              do_log_transformation = TRUE, log_base = 2,
              use_groups = FALSE, group_colours = NULL,
              normalization = "loess", lts_quantile = 0.8)


}

