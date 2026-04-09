# Example MaxQuant evidence file from the output of a TMT experiment

Experiment was performed by the Olsen lab and published on Nat. Commun.
(citation below).

## Usage

``` r
maxq_tmt_evidence
```

## Format

a data.table with 96 columns and 199 rows, the output of MaxQuant

## Details

Hogrebe, A., von Stechow, L., Bekker-Jensen, D.B. et al. Benchmarking
common quantification strategies for large-scale phosphoproteomics. Nat
Commun 9, 1045 (2018). https://doi.org/10.1038/s41467-018-03309-6

The experiment was processed using MaxQuant by the computational
proteomics team at Pfizer (Liang Xue and Pierre Jean).

The experiment did not contain a global profiling run, but we show an
example of extracting the unmodified peptides and using them in place of
the profiling run.

## Examples

``` r
head(maxq_tmt_evidence)
#>               Sequence Length Modifications
#> 1 AALLAQYADVTDEEDEADEK     20 Phospho (STY)
#> 2 AALLAQYADVTDEEDEADEK     20 Phospho (STY)
#> 3 AALLAQYADVTDEEDEADEK     20 Phospho (STY)
#> 4 AALLAQYADVTDEEDEADEK     20 Phospho (STY)
#> 5 AALLAQYADVTDEEDEADEK     20 Phospho (STY)
#> 6 AALLAQYADVTDEEDEADEK     20 Phospho (STY)
#>                       Modified.sequence Oxidation..M..Probabilities
#> 1 _AALLAQYADVT(Phospho (STY))DEEDEADEK_                            
#> 2 _AALLAQYADVT(Phospho (STY))DEEDEADEK_                            
#> 3 _AALLAQYADVT(Phospho (STY))DEEDEADEK_                            
#> 4 _AALLAQYADVT(Phospho (STY))DEEDEADEK_                            
#> 5 _AALLAQYADVT(Phospho (STY))DEEDEADEK_                            
#> 6 _AALLAQYADVT(Phospho (STY))DEEDEADEK_                            
#>          Phospho..STY..Probabilities Oxidation..M..Score.Diffs
#> 1            AALLAQYADVT(1)DEEDEADEK                          
#> 2 AALLAQY(0.316)ADVT(0.684)DEEDEADEK                          
#> 3            AALLAQYADVT(1)DEEDEADEK                          
#> 4 AALLAQY(0.001)ADVT(0.999)DEEDEADEK                          
#> 5 AALLAQY(0.001)ADVT(0.999)DEEDEADEK                          
#> 6            AALLAQYADVT(1)DEEDEADEK                          
#>         Phospho..STY..Score.Diffs Acetyl..N.term. Oxidation..M. Phospho..STY.
#> 1   AALLAQY(-74)ADVT(74)DEEDEADEK               0             0             1
#> 2 AALLAQY(-3.3)ADVT(3.3)DEEDEADEK               0             0             1
#> 3   AALLAQY(-75)ADVT(75)DEEDEADEK               0             0             1
#> 4   AALLAQY(-29)ADVT(29)DEEDEADEK               0             0             1
#> 5   AALLAQY(-29)ADVT(29)DEEDEADEK               0             0             1
#> 6   AALLAQY(-39)ADVT(39)DEEDEADEK               0             0             1
#>   Missed.cleavages Proteins Leading.proteins Leading.razor.protein Gene.names
#> 1                0   Q96MW1           Q96MW1                Q96MW1     CCDC43
#> 2                0   Q96MW1           Q96MW1                Q96MW1     CCDC43
#> 3                0   Q96MW1           Q96MW1                Q96MW1     CCDC43
#> 4                0   Q96MW1           Q96MW1                Q96MW1     CCDC43
#> 5                0   Q96MW1           Q96MW1                Q96MW1     CCDC43
#> 6                0   Q96MW1           Q96MW1                Q96MW1     CCDC43
#>                              Protein.names       Type
#> 1 Coiled-coil domain-containing protein 43 MULTI-MSMS
#> 2 Coiled-coil domain-containing protein 43       MSMS
#> 3 Coiled-coil domain-containing protein 43 MULTI-MSMS
#> 4 Coiled-coil domain-containing protein 43 MULTI-MSMS
#> 5 Coiled-coil domain-containing protein 43       MSMS
#> 6 Coiled-coil domain-containing protein 43       MSMS
#>                                        Raw.file Fraction Experiment MS.MS.m.z
#> 1 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1        1        MS2  912.0978
#> 2 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1        1        MS2  911.7875
#> 3 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_2        2        MS2  912.4321
#> 4 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_2        2        MS2  912.4330
#> 5 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_2        2        MS2  912.0959
#> 6 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_2        2        MS2  912.1791
#>   Charge      m.z     Mass Uncalibrated...Calibrated.m.z..ppm.
#> 1      3 759.3212 2274.942                              2.5755
#> 2      3 759.3212 2274.942                                 NaN
#> 3      3 759.3212 2274.942                              2.0210
#> 4      3 759.3212 2274.942                              2.2832
#> 5      3 759.3212 2274.942                                 NaN
#> 6      3 759.3212 2274.942                                 NaN
#>   Uncalibrated...Calibrated.m.z..Da. Mass.error..ppm. Mass.error..Da.
#> 1                          0.0019556           201200          152.78
#> 2                                NaN              NaN             NaN
#> 3                          0.0015346           201200          152.78
#> 4                          0.0017337           201200          152.78
#> 5                                NaN              NaN             NaN
#> 6                                NaN              NaN             NaN
#>   Uncalibrated.mass.error..ppm. Uncalibrated.mass.error..Da.
#> 1                        201200                       152.78
#> 2                           NaN                          NaN
#> 3                        201200                       152.78
#> 4                        201200                       152.78
#> 5                           NaN                          NaN
#> 6                           NaN                          NaN
#>   Max.intensity.m.z.0 Retention.time Retention.length Calibrated.retention.time
#> 1            912.4306         115.23          2.47460                    115.23
#> 2                 NaN         118.50          1.00000                    118.50
#> 3            912.4308         114.72          2.35110                    114.72
#> 4            912.4312         116.72          0.53974                    116.72
#> 5                 NaN         116.50          1.00000                    116.50
#> 6                 NaN         117.68          1.00000                    117.68
#>   Calibrated.retention.time.start Calibrated.retention.time.finish
#> 1                          114.78                           117.25
#> 2                          118.00                           119.00
#> 3                          114.22                           116.57
#> 4                          116.57                           117.11
#> 5                          116.00                           117.00
#> 6                          117.18                           118.18
#>   Retention.time.calibration Match.time.difference Match.m.z.difference
#> 1                          0                    NA                   NA
#> 2                          0                    NA                   NA
#> 3                          0                    NA                   NA
#> 4                          0                    NA                   NA
#> 5                          0                    NA                   NA
#> 6                          0                    NA                   NA
#>   Match.q.value Match.score Number.of.data.points Number.of.scans
#> 1            NA          NA                   382             108
#> 2            NA          NA                    NA              NA
#> 3            NA          NA                   356             106
#> 4            NA          NA                    42              23
#> 5            NA          NA                    NA              NA
#> 6            NA          NA                    NA              NA
#>   Number.of.isotopic.peaks       PIF Fraction.of.total.spectrum
#> 1                        6 0.9277212                0.004005980
#> 2                       NA       NaN                        NaN
#> 3                        6 0.9841691                0.006051585
#> 4                        3 0.7940387                0.001093331
#> 5                       NA       NaN                        NaN
#> 6                       NA       NaN                        NaN
#>   Base.peak.fraction        PEP MS.MS.count MS.MS.scan.number   Score
#> 1         0.09994941 1.8662e-46           2             38414 189.110
#> 2                NaN 1.6689e-02           1             39380  47.726
#> 3         0.08806564 5.8955e-39           2             38486 176.840
#> 4         0.02217602 5.4630e-08           1             38940 104.880
#> 5                NaN 8.6608e-27           1             38878 127.050
#> 6                NaN 4.8577e-26           1             39319 122.300
#>   Delta.score Combinatorics Intensity Reporter.intensity.corrected.1
#> 1     180.010             2 240530000                         125940
#> 2      30.847             2        NA                              0
#> 3     167.750             2 102250000                         155390
#> 4      89.116             2  11613000                          69254
#> 5     118.860             2        NA                              0
#> 6     114.110             2        NA                              0
#>   Reporter.intensity.corrected.2 Reporter.intensity.corrected.3
#> 1                         133100                         113620
#> 2                              0                              0
#> 3                         164220                         143320
#> 4                          80370                          60918
#> 5                              0                              0
#> 6                              0                              0
#>   Reporter.intensity.corrected.4 Reporter.intensity.corrected.5
#> 1                         119420                         108710
#> 2                              0                              0
#> 3                         160600                         163330
#> 4                          75915                          67276
#> 5                              0                              0
#> 6                              0                              0
#>   Reporter.intensity.corrected.6 Reporter.intensity.corrected.7
#> 1                         155830                         138320
#> 2                              0                              0
#> 3                         205080                         184670
#> 4                          86921                          80510
#> 5                              0                              0
#> 6                              0                              0
#>   Reporter.intensity.corrected.8 Reporter.intensity.corrected.9
#> 1                         158440                         174430
#> 2                              0                              0
#> 3                         178690                         203130
#> 4                          80892                          87856
#> 5                              0                              0
#> 6                              0                              0
#>   Reporter.intensity.corrected.10 Reporter.intensity.1 Reporter.intensity.2
#> 1                          130450               125940               133100
#> 2                               0                    0                    0
#> 3                          155170               155390               164220
#> 4                           82300                69254                80370
#> 5                               0                    0                    0
#> 6                               0                    0                    0
#>   Reporter.intensity.3 Reporter.intensity.4 Reporter.intensity.5
#> 1               113620               119420               108710
#> 2                    0                    0                    0
#> 3               143320               160600               163330
#> 4                60918                75915                67276
#> 5                    0                    0                    0
#> 6                    0                    0                    0
#>   Reporter.intensity.6 Reporter.intensity.7 Reporter.intensity.8
#> 1               155830               138320               158440
#> 2                    0                    0                    0
#> 3               205080               184670               178690
#> 4                86921                80510                80892
#> 5                    0                    0                    0
#> 6                    0                    0                    0
#>   Reporter.intensity.9 Reporter.intensity.10 Reporter.intensity.count.1
#> 1               174430                130450                          2
#> 2                    0                     0                          0
#> 3               203130                155170                          2
#> 4                87856                 82300                          1
#> 5                    0                     0                          0
#> 6                    0                     0                          0
#>   Reporter.intensity.count.2 Reporter.intensity.count.3
#> 1                          2                          2
#> 2                          0                          0
#> 3                          2                          2
#> 4                          1                          1
#> 5                          0                          0
#> 6                          0                          0
#>   Reporter.intensity.count.4 Reporter.intensity.count.5
#> 1                          2                          2
#> 2                          0                          0
#> 3                          2                          2
#> 4                          1                          1
#> 5                          0                          0
#> 6                          0                          0
#>   Reporter.intensity.count.6 Reporter.intensity.count.7
#> 1                          2                          2
#> 2                          0                          0
#> 3                          2                          2
#> 4                          1                          1
#> 5                          0                          0
#> 6                          0                          0
#>   Reporter.intensity.count.8 Reporter.intensity.count.9
#> 1                          2                          2
#> 2                          0                          0
#> 3                          2                          2
#> 4                          1                          1
#> 5                          0                          0
#> 6                          0                          0
#>   Reporter.intensity.count.10 Reverse Potential.contaminant  id
#> 1                           2                               125
#> 2                           0                               126
#> 3                           2                               127
#> 4                           1                               128
#> 5                           0                               129
#> 6                           0                               130
#>   Protein.group.IDs Peptide.ID Mod..peptide.ID MS.MS.IDs Best.MS.MS
#> 1              2449         42              50   153;154        154
#> 2              2449         42              50       155        155
#> 3              2449         42              50   156;157        157
#> 4              2449         42              50       158        158
#> 5              2449         42              50       159        159
#> 6              2449         42              50       160        160
#>   Oxidation..M..site.IDs Phospho..STY..site.IDs Taxonomy.IDs
#> 1                                          6884           NA
#> 2                                          6884           NA
#> 3                                          6884           NA
#> 4                                          6884           NA
#> 5                                          6884           NA
#> 6                                          6884           NA
```
