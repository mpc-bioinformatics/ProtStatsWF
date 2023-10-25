readData <- function (data_path, na_strings, intensity_columns, zero_to_NA, log_data){
  
  D <- openxlsx::read.xlsx(data_path, na.strings = na_strings)
  
  id <- D[, -intensity_columns]
  D <- D[, intensity_columns]
  
  if(zero_to_NA) {
    D[D == 0] <- NA
  }
  
  if(log_data) {
    log_base <- NULL
    D <- log(D, base = log_base)
  }
  
  return (D)
}

makeGroups <- function (D, use_groups, group_colours){
  
  if (use_groups) {
    group <- factor(limma::strsplit2(colnames(D), "_")[,1])
  } else {
    group <- NULL
  }
  
  nr_groups <- length(levels(group))
  
  if (is.null(group_colours) & nr_groups >= 1) group_colours <- scales::hue_pal()(nr_groups)
  
  return(c(group, nr_groups, group_colours))
}