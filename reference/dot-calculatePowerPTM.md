# Power calculation for PTM experiment

Power calculation for PTM experiment

## Usage

``` r
.calculatePowerPTM(
  desiredFC,
  FDR,
  delta,
  ptm_median_sigma_error,
  protein_median_sigma_error,
  ptm_median_sigma_subject,
  protein_median_sigma_subject,
  numSample
)
```

## Arguments

- desiredFC:

  the range of a desired fold change which includes the lower and upper
  values of the desired fold change.

- FDR:

  a pre-specified false discovery ratio (FDR) to control the overall
  false positive rate. Default is 0.05

- numSample:

  minimal number of biological replicates per condition. TRUE represents
  you require to calculate the sample size for this category, else you
  should input the exact number of biological replicates.

## Value

`float` of power
