# Example annotation file for a label-free Spectronaut experiment.

Must be manually created by the user and input into the
SpectronauttoMSstatsPTMFormat converter. Requires the correct columns
and maps the experimental desing into the MSstats format. Specify unique
bioreplicates for group comparison designs, and the same bioreplicate
for repeated measure designs. The columns and descriptions are below.

## Usage

``` r
spectronaut_annotation
```

## Format

A data.table with 5 columns.

## Details

- Run : Run name that matches exactly with Spectronaut run. Used to join
  evidence and metadata in annotation file.

- Condition : Name of condition that was used for each run.

- BioReplicate : Name of biological replicate. Repeating the same name
  here will tell MSstatsPTM that the experiment is a repeated measure
  design.

- Raw.file : Run name that matches exactly with Spectronaut run. Used to
  join evidence and metadata in annotation file.

## Examples

``` r
head(spectronaut_annotation)
#>                                     Run Fraction TechRepMixture Condition
#> 1 20180815_QE3_nLC3_AH_DIA_Honly_ind_01        1              1   H100_Y0
#> 2 20180815_QE3_nLC3_AH_DIA_Honly_ind_02        1              1   H100_Y0
#> 3 20180815_QE3_nLC3_AH_DIA_Honly_ind_03        1              1   H100_Y0
#> 4 20180815_QE3_nLC3_AH_DIA_Yonly_ind_01        1              1   H0_Y100
#> 5 20180815_QE3_nLC3_AH_DIA_Yonly_ind_02        1              1   H0_Y100
#> 6 20180815_QE3_nLC3_AH_DIA_Yonly_ind_03        1              1   H0_Y100
#>   BioReplicate
#> 1   H100_Y0_04
#> 2   H100_Y0_05
#> 3   H100_Y0_06
#> 4   H0_Y100_01
#> 5   H0_Y100_02
#> 6   H0_Y100_03
```
