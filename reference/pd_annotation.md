# Example annotation file for a label-free Proteome Discoverer experiment.

Must be manually created by the user and input into the
PDtoMSstatsPTMFormat converter. Requires the correct columns and maps
the experimental desing into the MSstats format. Specify unique
bioreplicates for group comparison designs, and the same bioreplicate
for repeated measure designs. The columns and descriptions are below.

## Usage

``` r
pd_annotation
```

## Format

A data.table with 3 columns.

## Details

- Run : Run name that matches exactly with PD run. Used to join evidence
  and metadata in annotation file.

- Condition : Name of condition that was used for each run.

- BioReplicate : Name of biological replicate. Repeating the same name
  here will tell MSstatsPTM that the experiment is a repeated measure
  design.

## Examples

``` r
head(pd_annotation)
#>                                         Run Condition BioReplicate
#> 1 20180810_QE3_nLC3_AH_DDA_Yonly_ind_01.raw   H0_Y100   H0_Y100_01
#> 2 20180810_QE3_nLC3_AH_DDA_Yonly_ind_02.raw   H0_Y100   H0_Y100_02
#> 3 20180810_QE3_nLC3_AH_DDA_Yonly_ind_03.raw   H0_Y100   H0_Y100_03
#> 4 20180810_QE3_nLC3_AH_DDA_Honly_ind_01.raw   H100_Y0   H100_Y0_04
#> 5 20180810_QE3_nLC3_AH_DDA_Honly_ind_02.raw   H100_Y0   H100_Y0_05
#> 6 20180810_QE3_nLC3_AH_DDA_Honly_ind_03.raw   H100_Y0   H100_Y0_06
```
