# Output of FragPipe TMT global profiling experiment

This dataset was provided by the FragPipe team at the Nesvilab. It was
processed using Philosopher and targeted Phosphorylation.

## Usage

``` r
fragpipe_input_protein
```

## Format

A data.table with 27 columns and 47 rows.

## Examples

``` r
head(fragpipe_input_protein)
#>                                          Spectrum.Name
#> 1 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.03785.03785.3
#> 2 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.04553.04553.3
#> 3 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.05163.05163.3
#> 4 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.05368.05368.2
#> 5 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.06388.06388.3
#> 6 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.09515.09515.2
#>                                 Spectrum.File Peptide.Sequence
#> 1 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.mzML       ANGMELDGRR
#> 2 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.mzML         EEMDDQDK
#> 3 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.mzML       ANGMELDGRR
#> 4 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.mzML       EYEQDQSSSR
#> 5 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.mzML   DRDTQNLQAQEEER
#> 6 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01.mzML         EEMDDQDK
#>   Modified.Peptide.Sequence Probability Charge Protein.Start Protein.End   Gene
#> 1           ANGM[147]ELDGRR      0.8947      3           180         189  TRA2A
#> 2       n[230]EEM[147]DDQDK      0.8677      3           428         435 THRAP3
#> 3                                0.9888      3           180         189  TRA2A
#> 4          n[230]EYEQDQSSSR      1.0000      2           435         444 ZC3H13
#> 5      n[230]DRDTQNLQAQEEER      1.0000      3           166         179  SNIP1
#> 6            n[230]EEMDDQDK      1.0000      2           428         435 THRAP3
#>   Mapped.Genes               Protein Protein.ID       Mapped.Proteins
#> 1        TRA2B sp|Q13595|TRA2A_HUMAN     Q13595 sp|P62995|TRA2B_HUMAN
#> 2              sp|Q9Y2W1|TR150_HUMAN     Q9Y2W1                      
#> 3        TRA2B sp|Q13595|TRA2A_HUMAN     Q13595 sp|P62995|TRA2B_HUMAN
#> 4              sp|Q5T200|ZC3HD_HUMAN     Q5T200                      
#> 5              sp|Q8TAD8|SNIP1_HUMAN     Q8TAD8                      
#> 6              sp|Q9Y2W1|TR150_HUMAN     Q9Y2W1                      
#>                             Protein.Description Is.Unique Purity Intensity
#> 1           Transformer-2 protein homolog alpha     false   0.66  11404057
#> 2 Thyroid hormone receptor-associated protein 3      true   1.00 214722976
#> 3           Transformer-2 protein homolog alpha     false   1.00         0
#> 4 Zinc finger CCCH domain-containing protein 13      true   1.00         0
#> 5            Smad nuclear-interacting protein 1      true   1.00  23654240
#> 6 Thyroid hormone receptor-associated protein 3      true   1.00 137334144
#>   Channel.126 Channel.127N Channel.127C Channel.128N Channel.128C Channel.129N
#> 1        0.00         0.00         0.00         0.00         0.00         0.00
#> 2   169649.80    128647.68    211484.62    217940.97    285032.38    386665.06
#> 3        0.00         0.00         0.00         0.00         0.00         0.00
#> 4    27110.14     27206.84     18804.96     31722.93     46041.87     56579.66
#> 5    34456.11     45179.79     32160.84     45215.09     75652.51     31730.62
#> 6   540481.56    393964.88    647020.00    672371.25    932212.88   1208951.00
#>   Channel.129C Channel.130N Channel.130C Channel.131N
#> 1         0.00         0.00         0.00         0.00
#> 2    214397.91    235756.58    217907.38    128019.87
#> 3         0.00         0.00         0.00         0.00
#> 4     37022.87     38246.96     35628.89     29679.23
#> 5     46960.10     51399.19     38244.52     36727.02
#> 6    665127.81    719719.12    749849.31    429616.81
```
