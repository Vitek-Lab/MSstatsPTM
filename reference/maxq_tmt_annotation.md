# Example annotation file for a TMT MaxQuant experiment.

Must be manually created by the user and input into the
MaxQtoMSstatsPTMFormat converter. Requires the correct columns and maps
the experimental desing into the MSstats format. Specify unique
bioreplicates for group comparison designs, and the same bioreplicate
for repeated measure designs. The columns and descriptions are below.

## Usage

``` r
maxq_tmt_annotation
```

## Format

A data.table with 7 columns.

## Details

- Run : Run name that matches exactly with MaxQuant run. Used to join
  evidence and metadata in annotation file.

- Fraction : If multiple fractions were used (i.e. the same mixture
  split into multiple fractions) enter that here. TechRepMixture :
  Multiple runs using the same bioreplicate

- Channel : Mixture channel used

- Condition : Name of condition that was used for each run.

- Mixture : The unique mixture (plex) name

- BioReplicate : Name of biological replicate. Repeating the same name
  here will tell MSstatsPTM that the experiment is a repeated measure
  design.

## Examples

``` r
head(maxq_tmt_annotation)
#>                                             Run Fraction TechRepMixture
#> 1 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1        1              1
#> 2 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1        1              1
#> 3 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1        1              1
#> 4 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1        1              1
#> 5 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1        1              1
#> 6 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1        1              1
#>      Channel Condition  Mixture BioReplicate
#> 1 channel.10 yeast_01x Mixture1  yeast_01x_1
#> 2  channel.1 yeast_04x Mixture1  yeast_04x_1
#> 3  channel.2 yeast_10x Mixture1  yeast_10x_1
#> 4  channel.3 yeast_01x Mixture1  yeast_01x_2
#> 5  channel.4 yeast_04x Mixture1  yeast_04x_2
#> 6  channel.5 yeast_10x Mixture1  yeast_10x_2
```
