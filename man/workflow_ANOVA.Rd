% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflow_ANOVA.R
\name{workflow_ANOVA}
\alias{workflow_ANOVA}
\title{The workflow for ANOVA of quantitative proteomics data}
\usage{
workflow_ANOVA(
  data_path,
  output_path,
  intensity_columns,
  protein_names_column = NULL,
  paired = FALSE,
  var.equal = TRUE,
  log_before_test = TRUE,
  delog_for_FC = TRUE,
  p_value_zeros_to_min = TRUE,
  significant_after_FDR = TRUE,
  max_valid_values_off = 0,
  min_valid_values_on = NULL,
  suffix = "",
  plot_device = "pdf",
  plot_height = 10,
  plot_width = 10,
  plot_dpi = 300,
  volcano_base_size = 25
)
}
\arguments{
\item{data_path}{\strong{character} \cr
The path to an .xlsx file containing the input data.}

\item{output_path}{\strong{character} \cr
The path to the output folder.}

\item{intensity_columns}{\strong{integer vector} \cr
The numbers of the intensity columns in the table.}

\item{paired}{\strong{logical} \cr
If \code{TRUE}, a paired test will be done, otherwise an unpaired test.}

\item{var.equal}{\strong{logical} \cr
If \code{TRUE}, the variances are assumed to be equal.}

\item{log_before_test}{\strong{logical} \cr
If \code{TRUE}, the data will be log-transformed.}

\item{delog_for_FC}{\strong{logical} \cr
If \code{TRUE}, the fold change will be calculated without the log-transformation.}

\item{p_value_zeros_to_min}{\strong{logical} \cr
If \code{TRUE}, then \code{p_values == 0} will be set to the next smallest value of the p-values.}

\item{significant_after_FDR}{\strong{logical} \cr
If \code{TRUE}, candidates for the boxplots and heatmap need to be significant after FDR correction, otherwise all significant candidates will be used.}

\item{max_valid_values_off}{\strong{integer} \cr
The maximum number of valid values to be an off protein.}

\item{min_valid_values_on}{\strong{integer} \cr
The minimum number of valid values to be an on protein.}

\item{suffix}{\strong{character} \cr
The suffix of the file names should have one.}

\item{plot_device}{\strong{character} \cr
The type of the output file, e.g. "pdf" or "png".}

\item{plot_height}{\strong{numeric} \cr
The plot height in cm.}

\item{plot_width}{\strong{numeric} \cr
The plot width in cm.}

\item{plot_dpi}{\strong{integer} \cr
The "dots per inch" of the plot aka. the plot resolution.}
}
\value{
Message log of the workflow
}
\description{
The workflow for ANOVA of quantitative proteomics data
}
\examples{


# 1. Set the character of your data path, leading to an .xlsx file.
in_path <- "C:/Users/thisuser/Documents/dataFolder/data.xlsx"

# 2. Set the integer vector of the columns, which contain the intensities.
int_col <- 3:17

# 3. Set the character of the output path, leading to a folder for the results.
out_path <- "C:/Users/thisuser/Documents/resultsFolder/"

# 4. Run the ANOVA with the parameters you set.
\dontrun{
result <- workflow_ANOVA(data_path = in_path,
                         output_path = out_path,
                         intensity_columns = int_col) }

}
\seealso{
\code{\link[=workflow_ttest]{workflow_ttest()}} in case of only two groups in the sample.\cr
Functions used in this workflow:
\code{\link[=prepareTtestData]{prepareTtestData()}}, \code{\link[=ttest]{ttest()}}, \code{\link[=VolcanoPlot_ttest]{VolcanoPlot_ttest()}}, \code{\link[=pvalue_foldchange_histogram]{pvalue_foldchange_histogram()}},
\code{\link[=calculate_significance_categories_ttest]{calculate_significance_categories_ttest()}}, \code{\link[=Boxplots_candidates]{Boxplots_candidates()}},
\code{\link[=Heatmap_with_groups]{Heatmap_with_groups()}}, \code{\link[=calculate_onoff]{calculate_onoff()}}
}
