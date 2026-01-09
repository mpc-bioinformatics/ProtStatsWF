


workflow_batchnorm <- function(data_path,
                               intensity_columns,
                               batch,

                               filetype = "xlsx",
                               sep = ",",
                               dec = ".",
                               header = TRUE,
                               sheet = 1,
                               output_path,
                               output_type = "xlsx",
                               na_strings = c("NA", "NaN", "Filtered","#NV"),
                               na_out = "NA",

                               zero_to_NA = TRUE,
                               do_log_transformation = TRUE,
                               log_base = 2,

                               pre_norm = TRUE,
                               pre_norm_method = "loess"
                               ) {



  prepared_data <- prepareData(data_path = data_path, filetype = filetype,
                               sep = sep, dec = dec, header = header,
                               sheet = sheet,
                               intensity_columns = intensity_columns,
                               na_strings = na_strings, zero_to_NA = zero_to_NA,
                               do_log_transformation = do_log_transformation,
                               log_base = log_base, use_groups = FALSE,
                               group_colours = NULL, normalization = "nonorm",
                               lts_quantile = NA)

  D_norm <- regression_normalization(prepared_data$D,
                                     batch,
                                     #log = TRUE,
                                     pre_norm = pre_norm,
                                     pre_norm_method = pre_norm_method,
                                     return_coefs = TRUE)


  # Plots (valid values, MA, boxplots, PCA), if possible coloured by batch
  # new idea: maybe plot also the coefficients (= normalization factors?) as
  #                 a QC step


}
