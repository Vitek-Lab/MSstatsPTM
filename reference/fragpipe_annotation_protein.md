# Example annotation file for a global profiling run TMT FragPipe experiment.

Automatically created by FragPipe, manually checked by the user and
input into the FragPipetoMSstatsPTMFormat converter. Requires the
correct columns and maps the experimental desing into the MSstats
format. Specify unique bioreplicates for group comparison designs, and
the same bioreplicate for repeated measure designs. The columns and
descriptions are below.

## Usage

``` r
fragpipe_annotation_protein
```

## Format

A data.table with 7 columns.

## Details

- Run : Run name that matches exactly with FragPipe run. Used to join
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
head(fragpipe_annotation_protein)
#>                                      Run Fraction TechRepMixture Mixture
#> 1 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01        1              1  plex16
#> 2 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f02        2              1  plex16
#> 3 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f03        3              1  plex16
#> 4 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01        1              1  plex16
#> 5 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f02        2              1  plex16
#> 6 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f03        3              1  plex16
#>   Channel    BioReplicate Condition
#> 1     126 CPT0088900003_T         T
#> 2     126 CPT0088900003_T         T
#> 3     126 CPT0088900003_T         T
#> 4    127N CPT0079270003_T         T
#> 5    127N CPT0079270003_T         T
#> 6    127N CPT0079270003_T         T
```
