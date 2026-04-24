# Convert Spectronaut output into MSstatsPTM format

Converters label-free Spectronaut data into MSstatsPTM format. Requires
PSM output from Spectronaut and a custom made annotation file, mapping
the run name to the condition and bioreplicate. Can optionally take a
seperate PSM file for a global profiling run. If no global profiling run
provided, the function can extract the unmodified peptides from the PTM
PSM file and use them as a global profiling run (not recommended).

## Usage

``` r
SpectronauttoMSstatsPTMFormat(
  input,
  annotation = NULL,
  fasta_path = NULL,
  protein_input = NULL,
  annotation_protein = NULL,
  use_unmod_peptides = FALSE,
  intensity = "PeakArea",
  mod_id = "\\[Phospho \\(STY\\)\\]",
  fasta_protein_name = "uniprot_iso",
  remove_other_mods = TRUE,
  filter_with_Qvalue = TRUE,
  qvalue_cutoff = 0.01,
  useUniquePeptide = TRUE,
  removeFewMeasurements = TRUE,
  removeProtein_with1Feature = FALSE,
  summaryforMultipleRows = max,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = NULL
)
```

## Arguments

- input:

  name of Spectronaut PTM output, which is long-format. ProteinName,
  PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge,
  IsotopeLabelType, Condition, BioReplicate, Run, Intensity,
  F.ExcludedFromQuantification are required. Rows with
  F.ExcludedFromQuantification=True will be removed.

- annotation:

  name of 'annotation.txt' data which includes Condition, BioReplicate,
  Run. If annotation is already complete in Spectronaut, use
  annotation=NULL (default). It will use the annotation information from
  input.

- fasta_path:

  string containing path to the corresponding fasta file for the
  modified peptide dataset.

- protein_input:

  name of Spectronaut global protein output, which is as in the same
  format as `input` parameter.

- annotation_protein:

  name of annotation file for global protein data, in the same format as
  above.

- use_unmod_peptides:

  If `protein_input` is not provided, unmodified peptides can be
  extracted from `input` to be used in place of a global profiling run.
  Default is `FALSE`.

- intensity:

  'PeakArea'(default) uses not normalized peak area.
  'NormalizedPeakArea' uses peak area normalized by Spectronaut. Default
  is NULL

- mod_id:

  Character that indicates the modification of interest. Default is
  `\\(Phospho\\)`. Note `\\` must be included before special characters.

- fasta_protein_name:

  Name of fasta column that matches with protein name in evidence file.
  Default is `uniprot_iso`.

- remove_other_mods:

  Remove peptides which include modfications other than the one listed
  in `mod_id`. Default is `TRUE`. For example, in an experiment
  targeting Phosphorylation, setting this parameter to `TRUE` would
  remove peptides like (Acetyl (Protein N-term))AAAAPDSRVS(Phospho
  (STY))EEENLK. Set this parameter to `FALSE` to keep peptides with
  extraneous modifications.

- filter_with_Qvalue:

  TRUE(default) will filter out the intensities that have greater than
  qvalue_cutoff in EG.Qvalue column. Those intensities will be replaced
  with zero and will be considered as censored missing values for
  imputation purpose.

- qvalue_cutoff:

  Cutoff for EG.Qvalue. Default is 0.01.

- useUniquePeptide:

  TRUE (default) removes peptides that are assigned for more than one
  proteins. We assume to use unique peptide for each protein.

- removeFewMeasurements:

  TRUE (default) will remove the features that have 1 or 2 measurements
  across runs.

- removeProtein_with1Feature:

  TRUE will remove the proteins which have only 1 feature, which is the
  combination of peptide, precursor charge, fragment and charge. FALSE
  is default.

- summaryforMultipleRows:

  max(default) or sum - when there are multiple measurements for certain
  feature and certain run, use highest or sum of multiple intensities.

- use_log_file:

  logical. If TRUE, information about data processing will be saved to a
  file.

- append:

  logical. If TRUE, information about data processing will be added to
  an existing log file.

- verbose:

  logical. If TRUE, information about data processing wil be printed to
  the console.

- log_file_path:

  character. Path to a file to which information about data processing
  will be saved. If not provided, such a file will be created
  automatically. If 'append = TRUE', has to be a valid path to a file.

## Value

a list of two data.tables named 'PTM' and 'PROTEIN' in the format
required by MSstatsPTM.

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
head(spectronaut_annotation)
#>                                     Run Fraction TechRepMixture Condition
#> 1 20180815_QE3_nLC3_AH_DIA_Honly_ind_01        1              1   H100_Y0
#> 2 20180815_QE3_nLC3_AH_DIA_Honly_ind_02        1              1   H100_Y0
#> 3 20180815_QE3_nLC3_AH_DIA_Honly_ind_03        1              1   H100_Y0
#> 4 20180815_QE3_nLC3_AH_DIA_Yonly_ind_01        1              1   H0_Y100
#> 5 20180815_QE3_nLC3_AH_DIA_Yonly_ind_02        1              1   H0_Y100
#> 6 20180815_QE3_nLC3_AH_DIA_Yonly_ind_03        1              1   H0_Y100
#>   BioReplicate
#> 1   H100_Y0_04
#> 2   H100_Y0_05
#> 3   H100_Y0_06
#> 4   H0_Y100_01
#> 5   H0_Y100_02
#> 6   H0_Y100_03

msstats_input = SpectronauttoMSstatsPTMFormat(spectronaut_input, 
                  annotation=spectronaut_annotation, 
                  fasta_path=system.file("extdata", "spectronaut_fasta.fasta", package="MSstatsPTM"),
                  use_unmod_peptides=TRUE,
                  mod_id = "\\[Phospho \\(STY\\)\\]",
                  fasta_protein_name = "uniprot_iso"
                  )
#> INFO  [2026-04-24 19:22:57] ** Raw data from Spectronaut imported successfully.
#> INFO  [2026-04-24 19:22:57] ** Raw data from Spectronaut cleaned successfully.
#> INFO  [2026-04-24 19:22:57] ** Using provided annotation.
#> INFO  [2026-04-24 19:22:57] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-24 19:22:57] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be removed.
#> INFO  [2026-04-24 19:22:57] ** Intensities with values of FExcludedFromQuantification equal to TRUE are replaced with NA
#> WARN  [2026-04-24 19:22:57] ** PGQvalue not found in input columns.
#> WARN  [2026-04-24 19:22:57] ** EGQvalue not found in input columns.
#> INFO  [2026-04-24 19:22:57] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-24 19:22:57] ** Shared peptides are removed.
#> INFO  [2026-04-24 19:22:57] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: max
#> INFO  [2026-04-24 19:22:57] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-04-24 19:22:57] ** Run annotation merged with quantification data.
#> INFO  [2026-04-24 19:22:57] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-04-24 19:22:57] ** Fractionation handled.
#> INFO  [2026-04-24 19:22:57] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-24 19:22:57] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.

head(msstats_input$PTM)
#>    ProteinName                PeptideSequence PrecursorCharge FragmentIon
#> 37  P09938_S22 AAADALSDLEIKDS[Phospho (STY)]K               2        <NA>
#> 38  P09938_S22 AAADALSDLEIKDS[Phospho (STY)]K               2        <NA>
#> 39  P09938_S22 AAADALSDLEIKDS[Phospho (STY)]K               2        <NA>
#> 40  P09938_S22 AAADALSDLEIKDS[Phospho (STY)]K               2        <NA>
#> 41  P09938_S22 AAADALSDLEIKDS[Phospho (STY)]K               2        <NA>
#> 42  P09938_S22 AAADALSDLEIKDS[Phospho (STY)]K               2        <NA>
#>    ProductCharge IsotopeLabelType Condition BioReplicate
#> 37             2                L   H100_Y0   H100_Y0_04
#> 38             2                L   H100_Y0   H100_Y0_05
#> 39             2                L   H100_Y0   H100_Y0_06
#> 40             2                L   H0_Y100   H0_Y100_01
#> 41             2                L   H0_Y100   H0_Y100_02
#> 42             2                L   H0_Y100   H0_Y100_03
#>                                      Run Fraction Intensity
#> 37 20180815_QE3_nLC3_AH_DIA_Honly_ind_01        1        NA
#> 38 20180815_QE3_nLC3_AH_DIA_Honly_ind_02        1        NA
#> 39 20180815_QE3_nLC3_AH_DIA_Honly_ind_03        1        NA
#> 40 20180815_QE3_nLC3_AH_DIA_Yonly_ind_01        1 201390.72
#> 41 20180815_QE3_nLC3_AH_DIA_Yonly_ind_02        1  75962.33
#> 42 20180815_QE3_nLC3_AH_DIA_Yonly_ind_03        1 281808.72
head(msstats_input$PROTEIN)
#>   ProteinName PeptideSequence PrecursorCharge FragmentIon ProductCharge
#> 1      P36578     AAAAAAALQAK               2        <NA>             2
#> 2      P36578     AAAAAAALQAK               2        <NA>             2
#> 3      P36578     AAAAAAALQAK               2        <NA>             2
#> 4      P36578     AAAAAAALQAK               2        <NA>             2
#> 5      P36578     AAAAAAALQAK               2        <NA>             2
#> 6      P36578     AAAAAAALQAK               2        <NA>             2
#>   IsotopeLabelType Condition BioReplicate                                   Run
#> 1                L   H100_Y0   H100_Y0_04 20180815_QE3_nLC3_AH_DIA_Honly_ind_01
#> 2                L   H100_Y0   H100_Y0_05 20180815_QE3_nLC3_AH_DIA_Honly_ind_02
#> 3                L   H100_Y0   H100_Y0_06 20180815_QE3_nLC3_AH_DIA_Honly_ind_03
#> 4                L   H0_Y100   H0_Y100_01 20180815_QE3_nLC3_AH_DIA_Yonly_ind_01
#> 5                L   H0_Y100   H0_Y100_02 20180815_QE3_nLC3_AH_DIA_Yonly_ind_02
#> 6                L   H0_Y100   H0_Y100_03 20180815_QE3_nLC3_AH_DIA_Yonly_ind_03
#>   Fraction Intensity
#> 1        1        NA
#> 2        1        NA
#> 3        1        NA
#> 4        1        NA
#> 5        1        NA
#> 6        1        NA
```
