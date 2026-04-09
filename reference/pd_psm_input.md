# Example Proteome Discoverer evidence file from the output of a label free experiment

Experiment was performed by the Olsen lab and published on Nat. Commun.
(citation below).

## Usage

``` r
pd_psm_input
```

## Format

a data.table with 60 columns and 1657 rows, the output of PD

## Details

Bekker-Jensen, D.B., Bernhardt, O.M., Hogrebe, A. et al. Rapid and
site-specific deep phosphoproteome profiling by data-independent
acquisition without the need for spectral libraries. Nat Commun 11, 787
(2020). https://doi.org/10.1038/s41467-020-14609-1

The experiment was processed using Proteome Discoverer by the
computational proteomics team at Pfizer (Liang Xue and Pierre Jean).

The experiment did not contain a global profiling run, but we show an
example of extracting the unmodified peptides and using them in place of
the profiling run.

## Examples

``` r
head(pd_psm_input)
#>   Checked Tags Confidence Identifying.Node.Type Identifying.Node Search.ID
#> 1   False   NA       High            Sequest HT  Sequest HT (B4)         B
#> 2   False   NA       High            Sequest HT  Sequest HT (B4)         B
#> 3   False   NA       High            Sequest HT  Sequest HT (B4)         B
#> 4   False   NA       High            Sequest HT  Sequest HT (B4)         B
#> 5   False   NA       High            Sequest HT  Sequest HT (B4)         B
#> 6   False   NA       High            Sequest HT  Sequest HT (B4)         B
#>   Identifying.Node.No PSM.Ambiguity                     Sequence
#> 1                   4   Unambiguous AAEAGETGAATSATEGDNNNNTAAGDKK
#> 2                   4   Unambiguous                    LEAGLSDSK
#> 3                   4   Unambiguous                  LQETNPEEVPK
#> 4                   4   Unambiguous            LREENFSSNTSELGNKK
#> 5                   4   Unambiguous                   LLDNTNTDVK
#> 6                   4      Selected   HSDSYSENETNHTNVPISSTGGTNNK
#>                     Annotated.Sequence Modifications X..Proteins
#> 1 [K].AAEAGETGAATSATEGDNNNNTAAGDKK.[G]                         1
#> 2                    [K].LEAGLsDSK.[Q]   S6(Phospho)           1
#> 3                  [K].LQETNPEEVPK.[F]                         1
#> 4            [K].LREENFSSNtSELGNKK.[H]  T10(Phospho)           1
#> 5                   [K].LLDNTNTDVK.[I]                         1
#> 6   [R].HSDSYsENETNHTNVPISSTGGTNNK.[T]   S6(Phospho)          20
#>   Master.Protein.Accessions
#> 1                    Q07478
#> 2                    P17536
#> 3                    P35691
#> 4                    P33419
#> 5                    Q08421
#> 6                    Q99231
#>                                                                                                              Master.Protein.Descriptions
#> 1                   ATP-dependent RNA helicase SUB2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=SUB2 PE=1 SV=1
#> 2                                     Tropomyosin-1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TPM1 PE=1 SV=1
#> 3 Translationally-controlled tumor protein homolog OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TMA19 PE=1 SV=1
#> 4                        Spindle pole component 29 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=SPC29 PE=1 SV=1
#> 5             Enhancer of translation termination 1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=ETT1 PE=1 SV=1
#> 6        Transposon Ty1-DR3 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-DR3 PE=1 SV=2
#>                                                                                                                                               Protein.Accessions
#> 1                                                                                                                                                         Q07478
#> 2                                                                                                                                                         P17536
#> 3                                                                                                                                                         P35691
#> 4                                                                                                                                                         P33419
#> 5                                                                                                                                                         Q08421
#> 6 P0C2I9; P47100; Q12141; P47098; P0C2J0; P0C2I3; Q12269; Q03855; P0C2I2; Q99231; Q12414; Q07793; Q04711; Q04670; Q12273; Q03612; Q92393; P0C2I5; Q03434; O13527
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Protein.Descriptions
#> 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               ATP-dependent RNA helicase SUB2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=SUB2 PE=1 SV=1
#> 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Tropomyosin-1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TPM1 PE=1 SV=1
#> 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Translationally-controlled tumor protein homolog OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TMA19 PE=1 SV=1
#> 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    Spindle pole component 29 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=SPC29 PE=1 SV=1
#> 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         Enhancer of translation termination 1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=ETT1 PE=1 SV=1
#> 6 Transposon Ty1-PR1 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-PR1 PE=3 SV=1; Transposon Ty1-JR2 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-JR2 PE=3 SV=3; Transposon Ty1-GR1 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-GR1 PE=3 SV=1; Transposon Ty1-JR1 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-JR1 PE=3 SV=3; Transposon Ty1-PR2 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-PR2 PE=3 SV=1; Transposon Ty1-DR6 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-DR6 PE=3 SV=1; Transposon Ty1-GR2 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-GR2 PE=3 SV=3; Transposon Ty1-DR1 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-DR1 PE=3 SV=2; Transposon Ty1-DR5 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-DR5 PE=3 SV=1; Transposon Ty1-DR3 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-DR3 PE=1 SV=2; Transposon Ty1-PL Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-PL PE=3 SV=1; Transposon Ty1-DR4 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-DR4 PE=3 SV=2; Transposon Ty1-ML1 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-ML1 PE=3 SV=2; Transposon Ty1-MR2 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-MR2 PE=3 SV=2; Transposon Ty1-OL Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-OL PE=3 SV=1; Transposon Ty1-ER1 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-ER1 PE=3 SV=1; Transposon Ty1-OR Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-OR PE=3 SV=1; Transposon Ty1-LR2 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-LR2 PE=3 SV=1; Transposon Ty1-ML2 Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-ML2 PE=3 SV=2; Truncated transposon Ty1-A Gag-Pol polyprotein OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=TY1B-A PE=3 SV=3
#>   X..Missed.Cleavages Charge Original.Precursor.Charge DeltaScore DeltaCn Rank
#> 1                   1      3                         3     0.7340       0    1
#> 2                   0      2                         2     0.1038       0    1
#> 3                   0      2                         2     0.4962       0    1
#> 4                   2      3                         3     0.1122       0    1
#> 5                   0      2                         2     0.5480       0    1
#> 6                   0      3                         3     0.0091       0    1
#>   Search.Engine.Rank Concatenated.Rank m.z..Da. Contaminant Aligned.m.z..Da.
#> 1                  1                 1 879.3940       False         879.3940
#> 2                  1                 1 500.2237       False         500.2237
#> 3                  1                 1 642.3285       False         642.3285
#> 4                  1                 1 678.3174       False         678.3174
#> 5                  1                 1 566.7962       False         566.7962
#> 6                  1                 1 957.4013       False         957.4013
#>    MH...Da. Theo..MH...Da. DeltaM..ppm. Deltam.z..Da. Ions.Matched Matched.Ions
#> 1 2636.1673      2636.1667         0.25       0.00022          0/0            0
#> 2  999.4401       999.4394         0.64       0.00032          0/0            0
#> 3 1283.6497      1283.6478         1.47       0.00094          0/0            0
#> 4 2032.9376      2032.9335         2.00       0.00135          0/0            0
#> 5 1132.5851      1132.5844         0.57       0.00032          0/0            0
#> 6 2870.1893      2870.1861         1.14       0.00109          0/0            0
#>   Total.Ions Intensity Activation.Type NCE.... MS.Order
#> 1          0   7769382             HCD      28      MS2
#> 2          0   2978492             HCD      28      MS2
#> 3          0   2776155             HCD      28      MS2
#> 4          0   2369740             HCD      28      MS2
#> 5          0   2483459             HCD      28      MS2
#> 6          0   4993433             HCD      28      MS2
#>   Isolation.Interference.... Ion.Inject.Time..ms. RT..min. First.Scan Last.Scan
#> 1                   1.321094                   22   4.9337       1851      1851
#> 2                  25.545490                   22   7.5681       5317      5317
#> 3                  10.363530                   22   7.5687       5318      5318
#> 4                   6.098627                   22   7.5720       5320      5320
#> 5                  24.051600                   22   7.5726       5321      5321
#> 6                   0.000000                   22   7.5738       5323      5323
#>   Master.Scan.s.                             Spectrum.File File.ID Quan.Info
#> 1           1848 20180810_QE3_nLC3_AH_DDA_Yonly_ind_02.raw     F35          
#> 2           5306 20180810_QE3_nLC3_AH_DDA_Yonly_ind_02.raw     F35          
#> 3           5306 20180810_QE3_nLC3_AH_DDA_Yonly_ind_02.raw     F35          
#> 4           5319 20180810_QE3_nLC3_AH_DDA_Yonly_ind_02.raw     F35          
#> 5           5319 20180810_QE3_nLC3_AH_DDA_Yonly_ind_02.raw     F35          
#> 6           5319 20180810_QE3_nLC3_AH_DDA_Yonly_ind_02.raw     F35          
#>   Peptides.Matched XCorr X..Protein.Groups Percolator.q.Value Percolator.PEP
#> 1             3606  8.57                 1       0.0001142857   6.305117e-16
#> 2              679  2.60                 1       0.0001143000   3.711000e-04
#> 3              511  2.62                 1       0.0001143000   1.498000e-05
#> 4             1894  5.26                 1       0.0001143000   3.015000e-07
#> 5              505  2.50                 1       0.0001143000   7.418000e-05
#> 6             2549  7.66                 1       0.0001143000   6.116000e-14
#>   Percolator.SVMScore Precursor.Area Apex.RT..min.
#> 1            9.869511        5859362          4.94
#> 2            1.135000        5397380          7.58
#> 3            1.719000       66213752          7.61
#> 4            2.428000             NA            NA
#> 5            1.428000        1715639          7.56
#> 6            5.227000       22296560          7.59
#>   ptmRS..Binomial.Peptide.Score ptmRS..Isoform.Confidence.Probability
#> 1                            NA                                    NA
#> 2                        308.92                                 1.000
#> 3                            NA                                    NA
#> 4                        439.90                                 1.000
#> 5                            NA                                    NA
#> 6                        611.48                                 0.006
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ptmRS..Isoform.Group.Report
#> 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
#> 2                                                                                                                                                                                                                                                                                                          y1+-H(3) O(4) P(8): 788.38; y1+-H(3) O(4) P(7): 659.34; y1+-H(3) O(4) P(6): 588.30; y1+-H(3) O(4) P(5): 531.28; y1+-H(3) O(4) P(4): 418.19; y2+-H(3) O(4) P(6): 294.65; y1+(3): 349.17; y1+(2): 234.14
#> 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
#> 4                                                                                                                                                                                                                                           b1+-H(3) O(4) P(10): 1160.53; b1+-H(3) O(4) P(12): 1376.61; y1+-H(3) O(4) P(10): 1059.54; y1+-H(3) O(4) P(9): 972.51; y1+-H(3) O(4) P(8): 858.47; y2+-H(3) O(4) P(10): 530.28; y2+-H(3) O(4) P(9): 486.76; y2+-H(3) O(4) P(8): 429.74; y1+(7): 775.43; y2+(7): 388.22
#> 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
#> 6 b1+-H(3) O(4) P(4): 409.15; b1+-H(3) O(4) P(5): 572.21; b1+-H(3) O(4) P(6): 659.24; b1+-H(3) O(4) P(9): 1031.37; b1+-H(3) O(4) P(10): 1132.42; b1+-H(3) O(4) P(11): 1246.46; b1+-H(3) O(4) P(12): 1383.52; b1+-H(3) O(4) P(14): 1598.61; b1+-H(3) O(4) P(15): 1697.68; b2+-H(3) O(4) P(9): 516.19; b2+-H(3) O(4) P(12): 692.26; b2+-H(3) O(4) P(13): 742.79; b2+-H(3) O(4) P(14): 799.81; b2+-H(3) O(4) P(15): 849.34; y2+-H(3) O(4) P(25): 1318.08; y2+-H(3) O(4) P(24): 1274.56; y2+-H(3) O(4) P(23): 1217.05
#>   ptmRS..Best.Site.Probabilities
#> 1                               
#> 2               S6(Phospho): 100
#> 3                               
#> 4              T10(Phospho): 100
#> 5                               
#> 6              S4(Phospho): 99.4
#>                                                                                     ptmRS..Phospho.Site.Probabilities
#> 1                                                                                                                    
#> 2                                                                                              S(6): 100.0; S(8): 0.0
#> 3                                                                                                                    
#> 4                                                                      S(7): 0.0; S(8): 0.0; T(10): 100.0; S(11): 0.0
#> 5                                                                                                                    
#> 6 S(2): 0.0; S(4): 99.4; Y(5): 0.0; S(6): 0.6; T(10): 0.0; T(13): 0.0; S(18): 0.0; S(19): 0.0; T(20): 0.0; T(23): 0.0
```
