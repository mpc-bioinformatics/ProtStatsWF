% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ValidValuePlot.R
\name{ValidValuePlot}
\alias{ValidValuePlot}
\title{Barplots showing the percentage of valid values for each sample}
\usage{
ValidValuePlot(
  D_long,
  use_groups = NULL,
  groupvar_name = "Group",
  group_colours = NULL,
  base_size = 15
)
}
\arguments{
\item{D_long}{A data.frame of the data set given in long format.}

\item{use_groups}{If \code{TRUE} data will be plotted in groups.}

\item{groupvar_name}{A character containing the name for the group variable.}

\item{group_colours}{A character vector of hex codes for the group colors.}

\item{base_size}{A numeric containing the base size of the font.}
}
\value{
a tibble and a ggplot of the valid values
}
\description{
Barplots showing the percentage of valid values for each sample
}
\examples{
\dontrun{
prepared_data <- prepareData(data_path = path, intensity_columns = intensity_cols)

vvplot <- ValidValuePlot(D_long = prepared_data[["D_long"]])
}

}
