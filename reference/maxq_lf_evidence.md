# Example MaxQuant evidence file from the output of a label free experiment

Experiment was performed by the Olsen lab and published on Nat. Commun.
(citation below).

## Usage

``` r
maxq_lf_evidence
```

## Format

a data.table with 63 columns and 511 rows, the output of MaxQuant

## Details

Bekker-Jensen, D.B., Bernhardt, O.M., Hogrebe, A. et al. Rapid and
site-specific deep phosphoproteome profiling by data-independent
acquisition without the need for spectral libraries. Nat Commun 11, 787
(2020). https://doi.org/10.1038/s41467-020-14609-1

The experiment was processed using MaxQuant by the computational
proteomics team at Pfizer (Liang Xue and Pierre Jean).

The experiment did not contain a global profiling run, but we show an
example of extracting the unmodified peptides and using them in place of
the profiling run.

## Examples

``` r
head(maxq_lf_evidence)
#>      Sequence Length Modifications Modified.sequence
#> 1 AAAAAAALQAK     11    Unmodified     _AAAAAAALQAK_
#> 2 AAAAAAALQAK     11    Unmodified     _AAAAAAALQAK_
#> 3 AAAAAAALQAK     11    Unmodified     _AAAAAAALQAK_
#> 4 AAAAAAALQAK     11    Unmodified     _AAAAAAALQAK_
#> 5 AAAAAAALQAK     11    Unmodified     _AAAAAAALQAK_
#> 6 AAAAAAALQAK     11    Unmodified     _AAAAAAALQAK_
#>   Oxidation..M..Probabilities Phospho..STY..Probabilities
#> 1                                                        
#> 2                                                        
#> 3                                                        
#> 4                                                        
#> 5                                                        
#> 6                                                        
#>   Oxidation..M..Score.Diffs Phospho..STY..Score.Diffs Acetyl..Protein.N.term.
#> 1                                                                           0
#> 2                                                                           0
#> 3                                                                           0
#> 4                                                                           0
#> 5                                                                           0
#> 6                                                                           0
#>   Oxidation..M. Phospho..STY. Missed.cleavages Proteins
#> 1             0             0                0   P36578
#> 2             0             0                0   P36578
#> 3             0             0                0   P36578
#> 4             0             0                0   P36578
#> 5             0             0                0   P36578
#> 6             0             0                0   P36578
#>                                                                      Leading.proteins
#> 1 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#> 2 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#> 3 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#> 4 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#> 5 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#> 6 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#>                                                                 Leading.razor.protein
#> 1 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#> 2 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#> 3 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#> 4 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#> 5 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#> 6 sp|P36578|RL4_HUMAN60SribosomalproteinL4OS=Homosapiens(Human)OX=9606GN=RPL4PE=1SV=5
#>         Type                             Raw.file  Experiment MS.MS.m.z Charge
#> 1 MULTI-MSMS 20180810_QE3_nLC3_AH_DDA_H100_Y25_01 H100_Y25_01  478.7820      2
#> 2 MULTI-MSMS 20180810_QE3_nLC3_AH_DDA_H100_Y25_03 H100_Y25_03  478.7810      2
#> 3 MULTI-MSMS 20180810_QE3_nLC3_AH_DDA_H100_Y25_04 H100_Y25_04  478.7803      2
#> 4 MULTI-MSMS 20180810_QE3_nLC3_AH_DDA_H100_Y25_05 H100_Y25_05  478.7815      2
#> 5 MULTI-MSMS 20180810_QE3_nLC3_AH_DDA_H100_Y25_06 H100_Y25_06  478.7802      2
#> 6 MULTI-MSMS 20180810_QE3_nLC3_AH_DDA_H100_Y50_01 H100_Y50_01  478.7803      2
#>        m.z     Mass Uncalibrated...Calibrated.m.z..ppm.
#> 1 478.7798 955.5451                              5.0066
#> 2 478.7798 955.5451                              2.4819
#> 3 478.7798 955.5451                              2.1796
#> 4 478.7798 955.5451                              2.2967
#> 5 478.7798 955.5451                              2.9534
#> 6 478.7798 955.5451                              1.6198
#>   Uncalibrated...Calibrated.m.z..Da. Mass.error..ppm. Mass.error..Da.
#> 1                         0.00239710         -0.46202     -0.00022120
#> 2                         0.00118830          0.76954      0.00036844
#> 3                         0.00104350         -0.27332     -0.00013086
#> 4                         0.00109960          0.50609      0.00024231
#> 5                         0.00141400         -0.77230     -0.00036976
#> 6                         0.00077553         -0.92114     -0.00044102
#>   Uncalibrated.mass.error..ppm. Uncalibrated.mass.error..Da.
#> 1                       4.54460                   0.00217590
#> 2                       3.25140                   0.00155670
#> 3                       1.90620                   0.00091267
#> 4                       2.80280                   0.00134190
#> 5                       2.18110                   0.00104430
#> 6                       0.69867                   0.00033451
#>   Max.intensity.m.z.0 Retention.time Retention.length Calibrated.retention.time
#> 1            478.7797         6.2039         0.079819                    6.2039
#> 2            478.7799         6.2822         0.090469                    6.2822
#> 3            478.7796         6.2327         0.110390                    6.2327
#> 4            478.7799         6.2924         0.090472                    6.2924
#> 5            478.7793         6.9938         0.089908                    6.9938
#> 6            478.7792         6.1802         0.110780                    6.1802
#>   Calibrated.retention.time.start Calibrated.retention.time.finish
#> 1                          6.1683                           6.2481
#> 2                          6.2257                           6.3162
#> 3                          6.1950                           6.3054
#> 4                          6.2470                           6.3375
#> 5                          6.9485                           7.0384
#> 6                          6.1278                           6.2385
#>   Retention.time.calibration Match.time.difference Match.m.z.difference
#> 1                          0                    NA                   NA
#> 2                          0                    NA                   NA
#> 3                          0                    NA                   NA
#> 4                          0                    NA                   NA
#> 5                          0                    NA                   NA
#> 6                          0                    NA                   NA
#>   Match.q.value Match.score Number.of.data.points Number.of.scans
#> 1            NA          NA                    15               7
#> 2            NA          NA                    14               8
#> 3            NA          NA                    19              10
#> 4            NA          NA                    18               8
#> 5            NA          NA                    17               8
#> 6            NA          NA                    17              10
#>   Number.of.isotopic.peaks PIF Fraction.of.total.spectrum Base.peak.fraction
#> 1                        3   0                          0                  0
#> 2                        3   0                          0                  0
#> 3                        2   0                          0                  0
#> 4                        3   0                          0                  0
#> 5                        3   0                          0                  0
#> 6                        2   0                          0                  0
#>          PEP MS.MS.count MS.MS.scan.number   Score Delta.score Combinatorics
#> 1 1.2617e-03           1              4016  78.149      67.201             1
#> 2 8.8997e-05           1              4158 111.170      95.303             1
#> 3 3.4144e-04           1              4087  99.442      84.065             1
#> 4 1.2494e-03           1              4148  76.679      61.302             1
#> 5 6.5027e-05           1              4774 114.890      96.560             1
#> 6 8.7846e-05           1              3994 111.950      96.570             1
#>   Intensity Reverse Potential.contaminant id Protein.group.IDs Peptide.ID
#> 1   7589900                                0              1276          0
#> 2  11810000                                1              1276          0
#> 3  10223000                                2              1276          0
#> 4  10733000                                3              1276          0
#> 5  17840000                                4              1276          0
#> 6   9679200                                5              1276          0
#>   Mod..peptide.ID MS.MS.IDs Best.MS.MS Oxidation..M..site.IDs
#> 1               0         0          0                       
#> 2               0         1          1                       
#> 3               0         2          2                       
#> 4               0         3          3                       
#> 5               0         4          4                       
#> 6               0         5          5                       
#>   Phospho..STY..site.IDs Taxonomy.IDs
#> 1                                  NA
#> 2                                  NA
#> 3                                  NA
#> 4                                  NA
#> 5                                  NA
#> 6                                  NA
```
