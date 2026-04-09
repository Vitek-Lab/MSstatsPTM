# Output of FragPipe TMT PTM experiment

This dataset was provided by the FragPipe team at the Nesvilab. It was
processed using Philosopher and targeted Phosphorylation.

## Usage

``` r
fragpipe_input
```

## Format

A data.table with 29 columns and 246 rows.

## Examples

``` r
head(fragpipe_input)
#>                                          Spectrum.Name
#> 1 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.02743.02743.4
#> 2 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.02755.02755.4
#> 3 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.02812.02812.3
#> 4 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.02913.02913.3
#> 5 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.02920.02920.5
#> 6 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.02975.02975.3
#>                                 Spectrum.File Peptide.Sequence
#> 1 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.mzML    RRHSHSHSPMSTR
#> 2 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.mzML     RHSHSHSPMSTR
#> 3 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.mzML      HSHSHSPMSTR
#> 4 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.mzML        HTRDSEAQR
#> 5 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.mzML     QHREPSEQEHRR
#> 6 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01.mzML       RRSTSPDHTR
#>            Modified.Peptide.Sequence Probability Charge Protein.Start
#> 1 n[230]RRHSHS[167]HS[167]PM[147]STR      0.9953      4            92
#> 2  n[230]RHSHS[167]HS[167]PM[147]STR      0.9962      4            93
#> 3   n[230]HSHS[167]HS[167]PM[147]STR      0.9957      3            94
#> 4               n[230]HTRDS[167]EAQR      0.9778      3             6
#> 5            n[230]QHREPS[167]EQEHRR      0.9979      5           123
#> 6         n[230]RRS[167]T[181]SPDHTR      0.9976      3           636
#>   Protein.End    Gene Mapped.Genes               Protein Protein.ID
#> 1         104   TRA2B              sp|P62995|TRA2B_HUMAN     P62995
#> 2         104   TRA2B              sp|P62995|TRA2B_HUMAN     P62995
#> 3         104   TRA2B              sp|P62995|TRA2B_HUMAN     P62995
#> 4          14    STOM               sp|P27105|STOM_HUMAN     P27105
#> 5         134   SNIP1              sp|Q8TAD8|SNIP1_HUMAN     Q8TAD8
#> 6         645 AKAP17A              sp|Q02040|AK17A_HUMAN     Q02040
#>   Mapped.Proteins                Protein.Description Is.Unique Purity Intensity
#> 1                 Transformer-2 protein homolog beta      true   0.00         0
#> 2                 Transformer-2 protein homolog beta      true   1.00  93998240
#> 3                 Transformer-2 protein homolog beta      true   1.00  58713048
#> 4                                           Stomatin      true   0.53         0
#> 5                 Smad nuclear-interacting protein 1      true   0.84  32706532
#> 6                        A-kinase anchor protein 17A      true   1.00  28102432
#>               M.15.9949                                         STY.79.966331
#> 1 RRHSHSHSPM(1.0000)STR RRHS(0.1780)HS(0.8463)HS(0.8392)PMS(0.0709)T(0.0656)R
#> 2  RHSHSHSPM(1.0000)STR  RHS(0.1497)HS(0.8741)HS(0.8698)PMS(0.0523)T(0.0542)R
#> 3   HSHSHSPM(1.0000)STR   HS(0.0582)HS(0.9325)HS(0.9327)PMS(0.0387)T(0.0380)R
#> 4                                                   HT(0.3995)RDS(0.6005)EAQR
#> 5                                                        QHREPS(1.0000)EQEHRR
#> 6                                  RRS(0.7370)T(0.7512)S(0.4257)PDHT(0.0862)R
#>   Channel.126 Channel.127N Channel.127C Channel.128N Channel.128C Channel.129N
#> 1    5578.652     8280.212     7034.635    10747.431     14872.24     17204.29
#> 2   19045.867    25291.979    38326.629    34385.059     42117.77     72897.84
#> 3   18498.551    24321.078    33518.191    31881.815     36766.03     60230.07
#> 4   13825.080    15933.881     8398.121     8001.169     12493.29     22851.57
#> 5   13345.636    24715.256    11790.443    18234.275     34780.09     12546.47
#> 6   15176.378     8430.135    14684.991    19511.988     38792.40     58184.98
#>   Channel.129C Channel.130N Channel.130C Channel.131N
#> 1     15443.60    11442.942     12985.56     11235.75
#> 2     33277.25    50290.703     49428.89     26749.50
#> 3     31366.80    41944.098     47435.18     20533.39
#> 4     11001.78     8394.632     10014.59     14173.92
#> 5     29433.08    18376.287     14489.56     16896.90
#> 6     26905.22    26273.547     24920.71     15653.74
```
