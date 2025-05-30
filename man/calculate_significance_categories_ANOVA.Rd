% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers_ttest_ANOVA.R
\name{calculate_significance_categories_ANOVA}
\alias{calculate_significance_categories_ANOVA}
\title{Calculate significance categories for ANOVA.}
\usage{
calculate_significance_categories_ANOVA(
  p_posthoc,
  p_anova_adj,
  p_anova,
  fc,
  thres_fc = 2,
  thres_p = 0.05
)
}
\arguments{
\item{p_posthoc}{\strong{numeric vector} \cr
The posthoc p-values .}

\item{p_anova_adj}{\strong{numeric vector} \cr
The p-values after FDR-correction.}

\item{p_anova}{\strong{numeric vector} \cr
The p-values before FDR-correction}

\item{fc}{\strong{numeric vector} \cr
The values of the fold changes.}

\item{thres_fc}{\strong{numeric} \cr
The threshold for the fold changes.}

\item{thres_p}{\strong{numeric} \cr
The threshold for the p-values.}
}
\value{
A factor containing the significances.
}
\description{
Calculate significance categories for ANOVA.
}
