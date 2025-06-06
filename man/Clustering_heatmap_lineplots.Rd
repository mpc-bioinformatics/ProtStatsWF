% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WIP_Clustering_Heatmap_Lineplots.R
\name{Clustering_heatmap_lineplots}
\alias{Clustering_heatmap_lineplots}
\title{Clustering, Heatmap and Lineplots}
\usage{
Clustering_heatmap_lineplots(
  D,
  id,
  output_path,
  suffix = "",
  nr_clusters = NULL,
  row_split = TRUE,
  dist_method = "correlation",
  plot_height_heatmap = 15,
  plot_width_heatmap = 15,
  plot_height_lineplot = 10,
  plot_width_lineplot = 15
)
}
\arguments{
\item{D}{\strong{data.frame} \cr Dataframe with log-transformed protein intensities, e.g. filtered for significant proteins form the ANOVA or t-test results.}

\item{id}{\strong{data.frame} \cr dataframe with id information e.g. protein names, gene names, accessions etc.}

\item{output_path}{\strong{char(1)} \cr Path where results will be saved.}

\item{suffix}{\strong{char(1)} \cr Suffix for the file names, should start with a underscore.}

\item{nr_clusters}{\strong{int(1)} \cr Number of clusters. Default is NULL, meaning that the optimal number of clusters will be determined by \code{\link[dendextend:find_k]{dendextend::find_k()}}.}

\item{row_split}{\strong{logi(1)} \cr If TRUE, there will be space between row clusters in the heatmap.}

\item{dist_method}{\strong{char(1)} \cr distance method for clustering, default is "correlation" (centered Pearson correlation)}
}
\value{
save heatmap and data frame with cluster information, as well as line plots
}
\description{
Clustering, Heatmap and Lineplots
}
\examples{
# TODO
}
