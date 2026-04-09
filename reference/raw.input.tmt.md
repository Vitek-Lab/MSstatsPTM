# Example of input PTM dataset for TMT experiments.

It can be the output of MSstatsPTM converter MaxQtoMSstatsPTMFormat or
other MSstatsTMT converter functions (Please see MSstatsPTM_TMT_Workflow
vignette). The dataset is formatted as a list with two data.tables named
PTM and PROTEIN. In each data.table the variables are as follows:

## Usage

``` r
raw.input.tmt
```

## Format

A list of two data.tables named PTM and PROTEIN with 1716 and 29221 rows
respectively.

## Details

- ProteinName : Name of protein with modification site mapped in with an
  underscore. ie "Protein_4_Y474"

- PeptideSequence

- Charge

- PSM

- Mixture : Mixture of samples labeled with different TMT reagents,
  which can be analyzed in a single mass spectrometry experiment. If the
  channal doesn't have sample, please add
  `Empty' under Condition. \item TechRepMixture : Technical replicate of one mixture. One mixture may have multiple technical replicates. For example, if `TechRepMixture'
  = 1, 2 are the two technical replicates of one mixture, then they
  should match with same
  `Mixture' value. \item Run : MS run ID. \item Channel : Labeling information (126, ... 131). \item Condition : Condition (ex. Healthy, Cancer, Time0) \item BioReplicate : Unique ID for biological subject. If the channal doesn't have sample, please add `Empty'
  under BioReplicate.

- Intensity

## Examples

``` r
head(raw.input.tmt$PTM)
#>       ProteinName PeptideSequence Charge           PSM Mixture TechRepMixture
#> 1 Protein_12_S703      Peptide491      3 Peptide_491_3       1              1
#> 2 Protein_12_S703      Peptide491      3 Peptide_491_3       1              1
#> 3 Protein_12_S703      Peptide491      3 Peptide_491_3       1              1
#> 4 Protein_12_S703      Peptide491      3 Peptide_491_3       1              1
#> 5 Protein_12_S703      Peptide491      3 Peptide_491_3       1              1
#> 6 Protein_12_S703      Peptide491      3 Peptide_491_3       1              1
#>   Run Channel   Condition  BioReplicate Intensity
#> 1 1_1    128N Condition_2 Condition_2_1   48030.0
#> 2 1_1    129C Condition_4 Condition_4_2  100224.4
#> 3 1_1    131C Condition_3 Condition_3_2   66804.6
#> 4 1_1    130N Condition_1 Condition_1_2   46779.8
#> 5 1_1    128C Condition_6 Condition_6_1   77497.9
#> 6 1_1    126C Condition_4 Condition_4_1   81559.7
head(raw.input.tmt$PROTEIN)
#>   ProteinName PeptideSequence Charge             PSM Mixture TechRepMixture Run
#> 1  Protein_12     Peptide9121      3  Peptide_9121_3       1              1 1_1
#> 2  Protein_12    Peptide27963      5 Peptide_27963_5       1              1 1_1
#> 3  Protein_12    Peptide28482      4 Peptide_28482_4       1              1 1_1
#> 4  Protein_12    Peptide10940      2 Peptide_10940_2       2              1 2_1
#> 5  Protein_12     Peptide4900      2  Peptide_4900_2       2              1 2_1
#> 6  Protein_12     Peptide4900      3  Peptide_4900_3       2              1 2_1
#>   Channel   Condition  BioReplicate  Intensity
#> 1    126C Condition_4 Condition_4_1 10996116.9
#> 2    127C Condition_5 Condition_5_1    56965.1
#> 3    131N Condition_2 Condition_2_2   286121.7
#> 4    131N Condition_2 Condition_2_4   534806.0
#> 5    126C Condition_4 Condition_4_3  1134908.7
#> 6    126C Condition_4 Condition_4_3  1605773.2
```
