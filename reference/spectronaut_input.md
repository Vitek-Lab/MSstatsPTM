# Example Spectronaut evidence file from the output of a label free experiment

Experiment was performed by the Olsen lab and published on Nat. Commun.
(citation below).

## Usage

``` r
spectronaut_input
```

## Format

a data.table with 23 columns and 2683 rows, the output of Spectronaut

## Details

Bekker-Jensen, D.B., Bernhardt, O.M., Hogrebe, A. et al. Rapid and
site-specific deep phosphoproteome profiling by data-independent
acquisition without the need for spectral libraries. Nat Commun 11, 787
(2020). https://doi.org/10.1038/s41467-020-14609-1

The experiment was processed using Spectronaut by the computational
proteomics team at Pfizer (Liang Xue and Pierre Jean).

The experiment did not contain a global profiling run, but we show an
example of extracting the unmodified peptides and using them in place of
the profiling run.

## Examples

``` r
head(spectronaut_input)
#>   R.Condition                            R.FileName R.Replicate PG.Genes
#> 1       Honly 20180815_QE3_nLC3_AH_DIA_Honly_ind_01           1    DNM1L
#> 2       Honly 20180815_QE3_nLC3_AH_DIA_Honly_ind_01           1     BIN1
#> 3       Honly 20180815_QE3_nLC3_AH_DIA_Honly_ind_01           1     BIN1
#> 4       Honly 20180815_QE3_nLC3_AH_DIA_Honly_ind_01           1    KMT2D
#> 5       Honly 20180815_QE3_nLC3_AH_DIA_Honly_ind_01           1 PPP1R12A
#> 6       Honly 20180815_QE3_nLC3_AH_DIA_Honly_ind_01           1   SEC16A
#>                         PG.ProteinDescriptions PG.ProteinGroups PG.ProteinNames
#> 1                       Dynamin-1-like protein           O00429     DNM1L_HUMAN
#> 2      Myc box-dependent-interacting protein 1           O00499      BIN1_HUMAN
#> 3      Myc box-dependent-interacting protein 1           O00499      BIN1_HUMAN
#> 4        Histone-lysine N-methyltransferase 2D           O14686     KMT2D_HUMAN
#> 5 Protein phosphatase 1 regulatory subunit 12A           O14974     MYPT1_HUMAN
#> 6             Protein transport protein Sec16A           O15027     SC16A_HUMAN
#>   PEP.PeptidePosition EG.IsDecoy
#> 1                 607      FALSE
#> 2                 293      FALSE
#> 3                 293      FALSE
#> 4                4736      FALSE
#> 5                 443      FALSE
#> 6                 879      FALSE
#>                                             EG.PrecursorId
#> 1          _SKPIPIM[Oxidation (M)]PAS[Phospho (STY)]PQK_.2
#> 2                   _GNKSPS[Phospho (STY)]PPDGSPAATPEIR_.3
#> 3                   _GNKS[Phospho (STY)]PSPPDGSPAATPEIR_.3
#> 4                           _ALS[Phospho (STY)]PVIPLIPR_.2
#> 5                        _TGS[Phospho (STY)]YGALAEITASK_.2
#> 6 _AQQELVPPQQQ[Deamidation (NQ)]AS[Phospho (STY)]PPQLPK_.3
#>   EG.PTMAssayCandidateScore EG.PTMAssayProbability EG.PTMLocalizationConfidence
#> 1                 29.064455              0.9999999                    0.9999999
#> 2                  6.009665              0.4966855                    0.4966855
#> 3                  6.009665              0.4966855                    0.4966855
#> 4                       NaN                    NaN                    1.0000000
#> 5                 24.307762              0.5848936                    0.5848936
#> 6                 15.286304              0.3310838                    0.6655419
#>                                                                                                                                                             EG.PTMLocalizationProbabilities
#> 1                                                                                                              _S[Phospho (STY): 0%]KPIPIM[Oxidation (M): 100%]PAS[Phospho (STY): 100%]PQK_
#> 2                                                                               _GNKS[Phospho (STY): 49.7%]PS[Phospho (STY): 49.7%]PPDGS[Phospho (STY): 0.3%]PAAT[Phospho (STY): 0.3%]PEIR_
#> 3                                                                               _GNKS[Phospho (STY): 49.7%]PS[Phospho (STY): 49.7%]PPDGS[Phospho (STY): 0.3%]PAAT[Phospho (STY): 0.3%]PEIR_
#> 4                                                                                                                                                        _ALS[Phospho (STY): 100%]PVIPLIPR_
#> 5                                                                     _T[Phospho (STY): 41.5%]GS[Phospho (STY): 58.5%]Y[Phospho (STY): 0%]GALAEIT[Phospho (STY): 0%]AS[Phospho (STY): 0%]K_
#> 6 _AQ[Deamidation (NQ): 0%]Q[Deamidation (NQ): 0%]ELVPPQ[Deamidation (NQ): 33.1%]Q[Deamidation (NQ): 33.1%]Q[Deamidation (NQ): 33.1%]AS[Phospho (STY): 100%]PPQ[Deamidation (NQ): 0.7%]LPK_
#>   EG.NormalizationFactor EG.TotalQuantity..Settings. FG.Charge F.Charge
#> 1              1558039.8                 24796966912         2        2
#> 2              1766846.4                  9006234624         3        3
#> 3              1766846.4                  9006234624         3        3
#> 4               887979.7                 26086424576         2        2
#> 5              1128734.0                 98855960576         2        2
#> 6              1189332.0                356193000000         3        3
#>                                      EG.ModifiedSequence F.FrgIon F.FrgLossType
#> 1          _SKPIPIM[Oxidation (M)]PAS[Phospho (STY)]PQK_       NA        noloss
#> 2                   _GNKSPS[Phospho (STY)]PPDGSPAATPEIR_       NA        noloss
#> 3                   _GNKS[Phospho (STY)]PSPPDGSPAATPEIR_       NA        noloss
#> 4                           _ALS[Phospho (STY)]PVIPLIPR_       NA        noloss
#> 5                        _TGS[Phospho (STY)]YGALAEITASK_       NA        noloss
#> 6 _AQQELVPPQQQ[Deamidation (NQ)]AS[Phospho (STY)]PPQLPK_       NA        noloss
#>   F.ExcludedFromQuantification   F.PeakArea
#> 1                        FALSE  24796966912
#> 2                        FALSE   9006234624
#> 3                        FALSE   9006234624
#> 4                        FALSE  26086424576
#> 5                        FALSE  98855960576
#> 6                        FALSE 356193000000
```
