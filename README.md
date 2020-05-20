# MSstatsPTM

The goal of MSstatsPTM is to provide the implementation of general statistical methods for characterization of quantitative changes in global post-translational modification profiling experiments.

## Installation

You can install the development version of MSstatsPTM from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tsunghengtsai/MSstatsPTM")
```

## Functionality

The MSstatsPTM package provides three categories of important functions for quantitative analysis of PTMs: summarization, estimation and comparison. 

### Summarization functions:

The `PTMsummarize` and related functions summarize log2-intensities of spectral features for each PTM site into one value per run.

### Estimation functions:

The `PTMestimate` and related functions take as input the summarized log2-intensities for each PTM site, performs statistical modeling for the abundance of the site, and returns the estimates of model parameters for all sites in all experimental conditions.

### Comparison functions:

The `PTMcompareMeans` and related functions perform significance analysis for detecting changes in PTM mean abundances between conditions.
