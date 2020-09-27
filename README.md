# MSstatsPTM

The goal of MSstatsPTM is to provide the implementation of general statistical methods for characterization of quantitative changes in global post-translational modification profiling experiments.

## Installation

You can install the development version of MSstatsPTM from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tsunghengtsai/MSstatsPTM")
```

## Functionality

Quantitative analyses of PTMs are supported by four main functions of MSstatsPTM.

### Normalization

`PTMnormalize()` normalizes the quantified peak intensities to correct systematic variation across MS runs.

### Summarization

`PTMsummarize()` summarizes log2-intensities of spectral features (i.e., precursor ions in DDA, fragments in DIA, or transitions in SRM) into one value per PTM site per run or one value per protein per run.

### Estimation

`PTMestimate()` takes as input the summarized log2-intensities for each PTM site, performs statistical modeling for the log2-abundance of the site, and returns the estimates of model parameters for all PTM sites in all experimental conditions.

### Comparison

`PTMcompareMeans()` performs statistical testing for detecting changes in PTM mean abundances between conditions.
