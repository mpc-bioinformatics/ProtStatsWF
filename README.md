# ProtStatsWF

This package contains workflows for the statistical analysis of proteomics data.

## Installation

You can install the development version of ProtStatsWF from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mpc-bioinformatics/ProtStatsWF")
library(ProtStatsWF)
```

New features will be introduced in github branches, which are merged as soon as the feature has been properly tested.
To install a specific branch of the package please use

``` r
# install.packages("devtools")
devtools::install_github("mpc-bioinformatics/ProtStatsWF", ref = "<branchname>")
library(ProtStatsWF)
```

## Workflows

The following workflows are implemented in ProtStatsWF version 0.1.0:

- `workflow_QC`: Quality control containing normalization, valid values, boxplots, MA-plots, PCA-plots
- `workflow_ttest`: Statistical comparison of two groups using a t-test (paired or unpaired) together with a volcano plot, histograms of p-values and fold changes as well as boxplots and a heatmap of the biomarker candidates
- `workflow_ANOVA`: Different versions of ANOVA together with volcano plots histograms of p-values and fold changes as well as boxplots and a heatmap of the biomarker candidates


The following workflows and features are planned for later versions:

- Adding the calculation of on/off proteins to `workflow_ttest` and `workflow_ANOVA`
- `workflow_clustering`: Hierarchical clustering of proteins based on their intensity patterns, visualized by a heatmap and lineplots of individual clusters
- `workflow_GO`: GO- and Pathway-analysis of a set of biomarker candidates with heatmaps of proteins belonging to a specific term


## Usage

The data sets have to be present as .xlsx files. The names of columns containing protein intensities have to follow the following scheme:
group1_01, group1_02, ...., group2_01, group2_02...

In the following examples the gold standard data set by Uszkoreit et. al, 2022 (PMID: 35845101) is used.

### `workflow_QC`

Normalize the data with the LOESS normalization method:

``` r
library(ProtStatsWF)

file <- system.file("extdata", "proteinGroups_D1.xlsx", package = "ProtStatsWF")

workflow_QC(
  data_path = file,
  filetype = "xlsx",
  sep = ",",
  dec = ".",
  header = TRUE,
  sheet = 1,
  output_path = "results",
  output_type = "xlsx",
  intensity_columns = 4:18,
  normalization_method = "loess",
  plot_device = "png",
  PCA_label = TRUE)


```

### `workflow_ttest`

Use the normalized data to perform a t-test. The data contains 5 groups, we will compare the groups `state1` and `state2` here:

``` r
library(ProtStatsWF)

file <- system.file("extdata", "D1_norm_ID.xlsx", package = "ProtStatsWF")

workflow_ttest(
  data_path = file,
  output_path = "results",
  intensity_columns = 4:9,
  log_before_test = FALSE,
  paired = FALSE,
  column_name_protein = "Gene.names",
  plot_device = "png")

```

### `workflow_ANOVA`

Use the normalized data to perform an ANOVA. Here we can compare all 5 groups directly:

``` r
library(ProtStatsWF)

file <- system.file("extdata", "D1_norm_ID.xlsx", package = "ProtStatsWF")

workflow_ANOVA(
  data_path = file,
  output_path = "results",
  intensity_columns = 4:18,
  paired = FALSE,
  var.equal = FALSE,
  log_before_test = FALSE,
  delog_for_FC = TRUE,
  protein_names_column = "Gene.names",
  plot_device = "png")
```


## Funding

The development and maintanence of ProtStatsWF is funded by CUBiMed.RUB
(<https://www.cubimed.ruhr-uni-bochum.de/index.html.en>). We also also
other cool tools and consulting for statistics, bioninformatics and
machine learning!
