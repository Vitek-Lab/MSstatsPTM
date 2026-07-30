# Convert output of label-free or TMT MaxQuant experiments into MSstatsPTM format

Takes as input LF/TMT experiments from MaxQ and converts the data into
the format needed for MSstatsPTM. Requires modified evidence.txt file
from MaxQ and an annotation file for PTM data. To adjust modified
peptides for changes in global protein level, unmodified TMT
experimental data must also be returned. Optionally can use
`Phospho(STY)Sites.txt` (or other PTM specific files) from MaxQuant, but
this is not recommended. If PTM specific file provided, the raw
intensities must be provided, not a ratio.

## Usage

``` r
MaxQtoMSstatsPTMFormat(
  evidence = NULL,
  annotation = NULL,
  fasta_path,
  fasta_protein_name = "uniprot_ac",
  mod_id = "\\(Phospho \\(STY\\)\\)",
  sites_data = NULL,
  evidence_prot = NULL,
  proteinGroups = NULL,
  annotation_protein = NULL,
  use_unmod_peptides = FALSE,
  labeling_type = "LF",
  mod_num = "Single",
  TMT_keyword = "TMT",
  ptm_keyword = "phos",
  which_proteinid_ptm = "Proteins",
  which_proteinid_protein = "Proteins",
  remove_other_mods = TRUE,
  removeMpeptides = FALSE,
  removeOxidationMpeptides = FALSE,
  removeProtein_with1Peptide = FALSE,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = NULL
)
```

## Arguments

- evidence:

  name of 'evidence.txt' data, which includes feature-level data for
  enriched (PTM) data.

- annotation:

  data frame annotation file for the ptm level data. Contains column
  Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate,
  Condition.

- fasta_path:

  A string of path to a FASTA file, used to match PTM peptides.

- fasta_protein_name:

  Name of fasta column that matches with protein name in evidence file.
  Default is `uniprot_ac`.

- mod_id:

  Character that indicates the modification of interest. Default is
  `\\(Phospho\\)`. Note `\\` must be included before special characters.

- sites_data:

  (Not recommended. Only used if evidence file not provided. Only works
  for TMT labeled data) Modified peptide output from MaxQuant. For
  example, a phosphorylation experiment would require the
  Phospho(STY)Sites.txt file

- evidence_prot:

  name of 'evidence.txt' data, which includes feature-level data for
  global profiling (unmodified) data.

- proteinGroups:

  name of 'proteinGroups.txt' data. It needs to matching protein group
  ID in `evidence_prot`.

- annotation_protein:

  data frame annotation file for the protein level data. Contains column
  Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate,
  Condition.

- use_unmod_peptides:

  Boolean if the unmodified peptides in the input file should be used to
  construct the unmodified protein output. Only used if `input_protein`
  is not provided. Default is `FALSE`.

- labeling_type:

  Either `TMT` or `LF` (Label-Free) depending on experimental design.
  Default is `LF`.

- mod_num:

  (Only if `sites.data` is used) For modified peptide dataset. The
  number modifications per peptide to be used. If "Single", only
  peptides with one modification will be used. Otherwise "Total" can be
  selected which does not cap the number of modifications per peptide.
  "Single" is the default. Selecting "Total" may confound the effect of
  different modifications.

- TMT_keyword:

  (Only if `sites.data` is used) the sub-name of columns in sites.data
  file. Default is `TMT`. This corresponds to the columns in the format
  `Reporter.intensity.corrected.1.TMT1phos___1`. Specifically, this
  parameter indicates the first section of the string `TMT1phos` (Before
  the mixture number). If `TMT` is present in the string, set this value
  to `TMT`. Else if `TMT` is not there (ie string is in the format
  `1phos`) leave this parameter as an empty string (”).

- ptm_keyword:

  (Only if `sites.data` is used) the sub-name of columns in the
  sites.data file. Default is `phos`. This corresponds to the columns in
  the format `Reporter.intensity.corrected.1.TMT1phos___1`.
  Specifically, this parameter indicates the second section of the
  string `TMT1phos` (After the mixture number). If the string is
  present, set this parameter. Else if this part of the string is empty
  (ie string is in the format `TMT1`) leave this parameter as an empty
  string (”).

- which_proteinid_ptm:

  For PTM dataset, which column to use for protein name. Use
  'Proteins'(default) column for protein name. 'Leading.proteins' or
  'Leading.razor.protein' or 'Gene.names' can be used instead to get the
  protein ID with single protein. However, those can potentially have
  the shared peptides.

- which_proteinid_protein:

  For Protein dataset, which column to use for protein name. Same
  options as above.

- remove_other_mods:

  Remove peptides which include modfications other than the one listed
  in `mod_id`. Default is `TRUE`. For example, in an experiment
  targeting Phosphorylation, setting this parameter to `TRUE` would
  remove peptides like (Acetyl (Protein N-term))AAAAPDSRVS(Phospho
  (STY))EEENLK. Set this parameter to `FALSE` to keep peptides with
  extraneous modifications.

- removeMpeptides:

  If Oxidation (M) modifications should be removed. Default is TRUE.

- removeOxidationMpeptides:

  TRUE will remove the peptides including 'oxidation (M)' in
  modification. FALSE is default.

- removeProtein_with1Peptide:

  TRUE will remove the proteins which have only 1 peptide and charge.
  FALSE is default.

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
# TMT experiment
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

msstats_format_tmt = MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                        annotation=maxq_tmt_annotation,
                        fasta=system.file("extdata", "maxq_tmt_fasta.fasta", package="MSstatsPTM"),
                        fasta_protein_name="uniprot_ac",
                        mod_id="\\(Phospho \\(STY\\)\\)",
                        use_unmod_peptides=TRUE,
                        labeling_type = "TMT",
                        which_proteinid_ptm = "Proteins")
#> INFO  [2026-07-30 21:35:47] ** Raw data from MaxQuant imported successfully.
#> INFO  [2026-07-30 21:35:47] ** Rows with values of Potentialcontaminant equal to + are removed 
#> INFO  [2026-07-30 21:35:47] ** Rows with values of Reverse equal to + are removed 
#> INFO  [2026-07-30 21:35:47] ** Features with all missing measurements across channels within each run are removed.
#> INFO  [2026-07-30 21:35:47] ** Using provided annotation.
#> INFO  [2026-07-30 21:35:47] ** Run and Channel labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-07-30 21:35:47] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements within each run will be kept.
#> INFO  [2026-07-30 21:35:47] ** Features with all missing measurements across channels within each run are removed.
#> INFO  [2026-07-30 21:35:47] ** Shared peptides are removed.
#> INFO  [2026-07-30 21:35:47] ** Features with all missing measurements across channels within each run are removed.
#> INFO  [2026-07-30 21:35:47] ** PSMs have been aggregated to peptide ions.
#> INFO  [2026-07-30 21:35:47] ** Run annotation merged with quantification data.
#> INFO  [2026-07-30 21:35:47] ** Features with one or two measurements across channels within each run are removed.
#> INFO  [2026-07-30 21:35:47] ** Fractionation handled.
#> INFO  [2026-07-30 21:35:47] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-07-30 21:35:47] ** Finished preprocessing. The dataset is ready to be processed by the proteinSummarization function.

head(msstats_format_tmt$PTM)
#>   ProteinName                     PeptideSequence Charge
#> 1 Q96MW1_T139 AALLAQYADVT(Phospho (STY))DEEDEADEK      3
#> 2 Q96MW1_T139 AALLAQYADVT(Phospho (STY))DEEDEADEK      3
#> 3 Q96MW1_T139 AALLAQYADVT(Phospho (STY))DEEDEADEK      3
#> 4 Q96MW1_T139 AALLAQYADVT(Phospho (STY))DEEDEADEK      3
#> 5 Q96MW1_T139 AALLAQYADVT(Phospho (STY))DEEDEADEK      3
#> 6 Q96MW1_T139 AALLAQYADVT(Phospho (STY))DEEDEADEK      3
#>                                     PSM  Mixture TechRepMixture
#> 1 AALLAQYADVT(Phospho (STY))DEEDEADEK_3 Mixture1              1
#> 2 AALLAQYADVT(Phospho (STY))DEEDEADEK_3 Mixture1              1
#> 3 AALLAQYADVT(Phospho (STY))DEEDEADEK_3 Mixture1              1
#> 4 AALLAQYADVT(Phospho (STY))DEEDEADEK_3 Mixture1              1
#> 5 AALLAQYADVT(Phospho (STY))DEEDEADEK_3 Mixture1              1
#> 6 AALLAQYADVT(Phospho (STY))DEEDEADEK_3 Mixture1              1
#>                                             Run   Channel BioReplicate
#> 1 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1  channel1  yeast_04x_1
#> 2 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1 channel10  yeast_01x_1
#> 3 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1  channel2  yeast_10x_1
#> 4 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1  channel3  yeast_01x_2
#> 5 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1  channel4  yeast_04x_2
#> 6 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1  channel5  yeast_10x_2
#>   Condition Intensity
#> 1 yeast_04x    125940
#> 2 yeast_01x    130450
#> 3 yeast_10x    133100
#> 4 yeast_01x    113620
#> 5 yeast_04x    119420
#> 6 yeast_10x    108710
head(msstats_format_tmt$PROTEIN)
#>     ProteinName   PeptideSequence Charge                 PSM  Mixture
#> 351      P29966 EAGEGGEAEAPAAEGGK      2 EAGEGGEAEAPAAEGGK_2 Mixture1
#> 352      P29966 EAGEGGEAEAPAAEGGK      2 EAGEGGEAEAPAAEGGK_2 Mixture1
#> 353      P29966 EAGEGGEAEAPAAEGGK      2 EAGEGGEAEAPAAEGGK_2 Mixture1
#> 354      P29966 EAGEGGEAEAPAAEGGK      2 EAGEGGEAEAPAAEGGK_2 Mixture1
#> 355      P29966 EAGEGGEAEAPAAEGGK      2 EAGEGGEAEAPAAEGGK_2 Mixture1
#> 356      P29966 EAGEGGEAEAPAAEGGK      2 EAGEGGEAEAPAAEGGK_2 Mixture1
#>     TechRepMixture                                           Run   Channel
#> 351              1 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1  channel1
#> 352              1 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1 channel10
#> 353              1 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1  channel2
#> 354              1 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1  channel3
#> 355              1 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1  channel4
#> 356              1 20171106_LUMOS1_nLC13_AH_TechBench2_TMTMS2L_1  channel5
#>     BioReplicate Condition Intensity
#> 351  yeast_04x_1 yeast_04x     23384
#> 352  yeast_01x_1 yeast_01x     29612
#> 353  yeast_10x_1 yeast_10x     26450
#> 354  yeast_01x_2 yeast_01x     33341
#> 355  yeast_04x_2 yeast_04x     22335
#> 356  yeast_10x_2 yeast_10x     27212

# LF experiment
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
head(maxq_lf_annotation)
#>                                     Run Condition BioReplicate
#> 1 20180810_QE3_nLC3_AH_DDA_Yonly_ind_01   H0_Y100   H0_Y100_01
#> 2 20180810_QE3_nLC3_AH_DDA_Yonly_ind_02   H0_Y100   H0_Y100_02
#> 3 20180810_QE3_nLC3_AH_DDA_Yonly_ind_03   H0_Y100   H0_Y100_03
#> 4  20180810_QE3_nLC3_AH_DDA_H100_Y25_01  H100_Y25  H100_Y25_07
#> 5  20180810_QE3_nLC3_AH_DDA_H100_Y25_02  H100_Y25  H100_Y25_08
#> 6  20180810_QE3_nLC3_AH_DDA_H100_Y25_03  H100_Y25  H100_Y25_09
#>                                Raw.file IsotopeLabelType
#> 1 20180810_QE3_nLC3_AH_DDA_Yonly_ind_01                L
#> 2 20180810_QE3_nLC3_AH_DDA_Yonly_ind_02                L
#> 3 20180810_QE3_nLC3_AH_DDA_Yonly_ind_03                L
#> 4  20180810_QE3_nLC3_AH_DDA_H100_Y25_01                L
#> 5  20180810_QE3_nLC3_AH_DDA_H100_Y25_02                L
#> 6  20180810_QE3_nLC3_AH_DDA_H100_Y25_03                L

msstats_format_lf = MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                        annotation=maxq_lf_annotation,
                        fasta=system.file("extdata", "maxq_lf_fasta.fasta", package="MSstatsPTM"),
                        fasta_protein_name="uniprot_ac",
                        mod_id="\\(Phospho \\(STY\\)\\)",
                        use_unmod_peptides=TRUE,
                        labeling_type = "LF",
                        which_proteinid_ptm = "Proteins")
#> INFO  [2026-07-30 21:35:47] ** Raw data from MaxQuant imported successfully.
#> INFO  [2026-07-30 21:35:47] ** Rows with values of Potentialcontaminant equal to + are removed 
#> INFO  [2026-07-30 21:35:47] ** Rows with values of Reverse equal to + are removed 
#> INFO  [2026-07-30 21:35:47] ** Using provided annotation.
#> INFO  [2026-07-30 21:35:47] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-07-30 21:35:47] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be removed.
#> INFO  [2026-07-30 21:35:47] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-07-30 21:35:47] ** Shared peptides are removed.
#> INFO  [2026-07-30 21:35:47] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: max
#> INFO  [2026-07-30 21:35:47] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-07-30 21:35:47] ** Run annotation merged with quantification data.
#> INFO  [2026-07-30 21:35:47] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-07-30 21:35:47] ** Fractionation handled.
#> INFO  [2026-07-30 21:35:47] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-07-30 21:35:47] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.
head(msstats_format_lf$PTM)
#>       ProteinName                               PeptideSequence PrecursorCharge
#> 34 P09938_S15_S22 AAADALS(Phospho (STY))DLEIKDS(Phospho (STY))K               2
#> 35 P09938_S15_S22 AAADALS(Phospho (STY))DLEIKDS(Phospho (STY))K               2
#> 36 P09938_S15_S22 AAADALS(Phospho (STY))DLEIKDS(Phospho (STY))K               2
#> 37 P09938_S15_S22 AAADALS(Phospho (STY))DLEIKDS(Phospho (STY))K               2
#> 38 P09938_S15_S22 AAADALS(Phospho (STY))DLEIKDS(Phospho (STY))K               2
#> 39 P09938_S15_S22 AAADALS(Phospho (STY))DLEIKDS(Phospho (STY))K               2
#>    FragmentIon ProductCharge IsotopeLabelType Condition BioReplicate
#> 34          NA            NA                L H100_Y100 H100_Y100_19
#> 35          NA            NA                L H100_Y100 H100_Y100_20
#> 36          NA            NA                L H100_Y100 H100_Y100_21
#> 37          NA            NA                L H100_Y100 H100_Y100_22
#> 38          NA            NA                L H100_Y100 H100_Y100_23
#> 39          NA            NA                L H100_Y100 H100_Y100_24
#>                                      Run Fraction Intensity
#> 34 20180810_QE3_nLC3_AH_DDA_H100_Y100_01        1        NA
#> 35 20180810_QE3_nLC3_AH_DDA_H100_Y100_02        1        NA
#> 36 20180810_QE3_nLC3_AH_DDA_H100_Y100_03        1        NA
#> 37 20180810_QE3_nLC3_AH_DDA_H100_Y100_04        1        NA
#> 38 20180810_QE3_nLC3_AH_DDA_H100_Y100_05        1        NA
#> 39 20180810_QE3_nLC3_AH_DDA_H100_Y100_06        1        NA
head(msstats_format_lf$PROTEIN)
#>   ProteinName PeptideSequence PrecursorCharge FragmentIon ProductCharge
#> 1      P36578     AAAAAAALQAK               2          NA            NA
#> 2      P36578     AAAAAAALQAK               2          NA            NA
#> 3      P36578     AAAAAAALQAK               2          NA            NA
#> 4      P36578     AAAAAAALQAK               2          NA            NA
#> 5      P36578     AAAAAAALQAK               2          NA            NA
#> 6      P36578     AAAAAAALQAK               2          NA            NA
#>   IsotopeLabelType Condition BioReplicate                                   Run
#> 1                L H100_Y100 H100_Y100_19 20180810_QE3_nLC3_AH_DDA_H100_Y100_01
#> 2                L H100_Y100 H100_Y100_20 20180810_QE3_nLC3_AH_DDA_H100_Y100_02
#> 3                L H100_Y100 H100_Y100_21 20180810_QE3_nLC3_AH_DDA_H100_Y100_03
#> 4                L H100_Y100 H100_Y100_22 20180810_QE3_nLC3_AH_DDA_H100_Y100_04
#> 5                L H100_Y100 H100_Y100_23 20180810_QE3_nLC3_AH_DDA_H100_Y100_05
#> 6                L H100_Y100 H100_Y100_24 20180810_QE3_nLC3_AH_DDA_H100_Y100_06
#>   Fraction Intensity
#> 1        1  13697000
#> 2        1   8738100
#> 3        1  10827000
#> 4        1   9628900
#> 5        1   9485600
#> 6        1        NA
```
