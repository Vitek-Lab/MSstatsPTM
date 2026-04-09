# Example annotation file for a label-free MaxQuant experiment.

Must be manually created by the user and input into the
MaxQtoMSstatsPTMFormat converter. Requires the correct columns and maps
the experimental desing into the MSstats format. Specify unique
bioreplicates for group comparison designs, and the same bioreplicate
for repeated measure designs. The columns and descriptions are below.

## Usage

``` r
maxq_lf_annotation
```

## Format

A data.table with 5 columns.

## Details

- Run : Run name that matches exactly with MaxQuant run. Used to join
  evidence and metadata in annotation file.

- Condition : Name of condition that was used for each run.

- BioReplicate : Name of biological replicate. Repeating the same name
  here will tell MSstatsPTM that the experiment is a repeated measure
  design.

- Raw.file : Run name that matches exactly with MaxQuant run. Used to
  join evidence and metadata in annotation file.

- IsotopeLabelType: Name of isotope label. May be all `L` or unique
  depending on experimental design.

## Examples

``` r
head(maxq_lf_annotation)
#>                                     Run Condition BioReplicate
#> 1 20180810_QE3_nLC3_AH_DDA_Yonly_ind_01   H0_Y100   H0_Y100_01
#> 2 20180810_QE3_nLC3_AH_DDA_Yonly_ind_02   H0_Y100   H0_Y100_02
#> 3 20180810_QE3_nLC3_AH_DDA_Yonly_ind_03   H0_Y100   H0_Y100_03
#> 4  20180810_QE3_nLC3_AH_DDA_H100_Y25_01  H100_Y25  H100_Y25_07
#> 5  20180810_QE3_nLC3_AH_DDA_H100_Y25_02  H100_Y25  H100_Y25_08
#> 6  20180810_QE3_nLC3_AH_DDA_H100_Y25_03  H100_Y25  H100_Y25_09
#>                                Raw.file IsotopeLabelType
#> 1 20180810_QE3_nLC3_AH_DDA_Yonly_ind_01                L
#> 2 20180810_QE3_nLC3_AH_DDA_Yonly_ind_02                L
#> 3 20180810_QE3_nLC3_AH_DDA_Yonly_ind_03                L
#> 4  20180810_QE3_nLC3_AH_DDA_H100_Y25_01                L
#> 5  20180810_QE3_nLC3_AH_DDA_H100_Y25_02                L
#> 6  20180810_QE3_nLC3_AH_DDA_H100_Y25_03                L
```
