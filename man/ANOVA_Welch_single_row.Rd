% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ANOVA.R
\name{ANOVA_Welch_single_row}
\alias{ANOVA_Welch_single_row}
\title{Welch ANOVA (unequal variances)}
\usage{
ANOVA_Welch_single_row(
  x,
  group,
  min_obs_per_group,
  log_before_test = TRUE,
  log_base = 2,
  delog_for_FC = TRUE,
  min_perc_per_group = NULL
)
}
\arguments{
\item{x}{Numeric vector of protein intensities.}

\item{group}{Factor vector indicating the group membership.}

\item{min_obs_per_group}{Integer indicating the minimum number of observations per group.}

\item{log_before_test}{Logical indicating whether the data should be log-transformed before the test.}

\item{log_base}{Numeric indicating the base of the logarithm.}

\item{delog_for_FC}{Logical indicating whether the fold change should be calculated on the original scale.}

\item{min_perc_per_group}{Numeric indicating the minimum ratio of observations per group.}
}
\value{
Vector with p-values and fold changes.
}
\description{
Welch ANOVA (unequal variances)
}
\examples{
# TODO
}
