# Get sample size for PTM experiment

Get sample size for PTM experiment

## Usage

``` r
.getNumSamplePTM(
  desiredFC,
  power,
  alpha,
  delta,
  ptm_median_sigma_error,
  protein_median_sigma_error,
  ptm_median_sigma_subject,
  protein_median_sigma_subject
)
```

## Arguments

- desiredFC:

  the range of a desired fold change which includes the lower and upper
  values of the desired fold change.

- power:

  a pre-specified statistical power which defined as the probability of
  detecting a true fold change. TRUE represent you require to calculate
  the power for this category, else you should input the average of
  power you expect. Default is 0.9

## Value

`int` of samples
