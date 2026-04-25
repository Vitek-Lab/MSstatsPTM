# Convert output of TMT labeled Fragpipe data into MSstatsPTM format.

Takes as input TMT experiments which are the output of Fragpipe and
converts into MSstatsPTM format. Requires `msstats.csv` file and an
annotation file. Optionally an additional `msstats.csv` file can be
uploaded if a corresponding global profiling run was performed. Site
localization is performed and only high probability localizations are
kept.

## Usage

``` r
FragPipetoMSstatsPTMFormat(
  input,
  annotation = NULL,
  input_protein = NULL,
  annotation_protein = NULL,
  use_unmod_peptides = FALSE,
  label_type = "TMT",
  protein_id_col = "Protein",
  peptide_id_col = "Peptide.Sequence",
  mod_id_col = "STY",
  localization_cutoff = 0.75,
  remove_unlocalized_peptides = TRUE,
  Purity_cutoff = 0.6,
  PeptideProphet_prob_cutoff = 0.7,
  useUniquePeptide = TRUE,
  rmPSM_withfewMea_withinRun = FALSE,
  rmPeptide_OxidationM = TRUE,
  rmProtein_with1Feature = FALSE,
  summaryforMultipleRows = sum,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = NULL
)
```

## Arguments

- input:

  data.frame of `msstats.csv` file produced by Philosopher

- annotation:

  annotation with Run, Fraction, TechRepMixture, Mixture, Channel,
  BioReplicate, Condition columns or a path to file. Refer to the
  example 'annotation' for the meaning of each column. Channel column
  should be consistent with the channel columns (Ignore the prefix
  "Channel ") in msstats.csv file. Run column should be consistent with
  the Spectrum.File columns in msstats.csv file.

- input_protein:

  same as `input` for global profiling run. Default is NULL.

- annotation_protein:

  same as `annotation` for global profiling run. Default is NULL.

- use_unmod_peptides:

  Boolean if the unmodified peptides in the input file should be used to
  construct the unmodified protein output. Only used if `input_protein`
  is not provided. Default is `FALSE`.

- label_type:

  Type of labeling used for experiment. Must be one of "LF" or "TMT".
  Default is "TMT".

- protein_id_col:

  Use 'Protein'(default) column for TMT. This needs to be changed to
  "ProteinName" for label free. For TMT, 'Master.Protein.Accessions' can
  be used instead to get the protein ID with single protein.

- peptide_id_col:

  Use 'Peptide.Sequence'(default) column for TMT. Must be changed to
  "PeptideSequence" for label free. "Modified.Peptide.Sequence" can be
  used instead to get the modified peptide sequence.

- mod_id_col:

  Column containing the modified Amino Acids. For example, a
  Phosphorylation experiment may pass `STY`. The corresponding column
  with `STY` combined with the mass (e.x. `STY.79.9663`) will be
  selected. Default is `STY`.

- localization_cutoff:

  Minimum localization score required to keep modification. Default is
  .75.

- remove_unlocalized_peptides:

  Boolean indicating if peptides without all sites localized should be
  kept. Default is TRUE (non-localized sites will be removed).

- Purity_cutoff:

  Cutoff for purity. Default is 0.6. Purity refers to how much of the
  detected ion signal within a specific inclusion window belongs to the
  target molecule or its closely related forms, compared to any other
  unwanted signals or noise. Higher values indicate greater purity.

- PeptideProphet_prob_cutoff:

  Cutoff for the peptide identification probability. Default is 0.7. The
  probability is confidence score determined by PeptideProphet and
  higher values indicate greater confidence.

- useUniquePeptide:

  logical, if TRUE (default) removes peptides that are assigned for more
  than one proteins. We assume to use unique peptide for each protein.

- rmPSM_withfewMea_withinRun:

  TRUE will remove the features that have 1 or 2 measurements within
  each Run. Default is FALSE.

- rmPeptide_OxidationM:

  TRUE (default) will remove the peptides including oxidation (M)
  sequence.

- rmProtein_with1Feature:

  TRUE will remove the proteins which have only 1 peptide and charge.
  Default is FALSE.

- summaryforMultipleRows:

  sum (default) or max - when there are multiple measurements for
  certain feature in certain run, select the feature with the largest
  summation or maximal value.

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

`list` of one or two `data.frame` of class `MSstatsTMT`, named `PTM` and
`PROTEIN`

## Examples

``` r
# TMT Example (with global profiling run)
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
head(fragpipe_annotation)
#>                                      Run Fraction TechRepMixture Mixture
#> 1 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01        1              1  plex16
#> 2 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f02        2              1  plex16
#> 3 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f03        3              1  plex16
#> 4 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01        1              1  plex16
#> 5 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f02        2              1  plex16
#> 6 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f03        3              1  plex16
#>   Channel    BioReplicate Condition
#> 1     126 CPT0088900003_T         T
#> 2     126 CPT0088900003_T         T
#> 3     126 CPT0088900003_T         T
#> 4    127N CPT0079270003_T         T
#> 5    127N CPT0079270003_T         T
#> 6    127N CPT0079270003_T         T
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

msstats_data = FragPipetoMSstatsPTMFormat(fragpipe_input,
                                          fragpipe_annotation,
                                          fragpipe_input_protein, 
                                          fragpipe_annotation_protein,
                                          label_type="TMT",
                                          mod_id_col = "STY",
                                          localization_cutoff=.75,
                                          remove_unlocalized_peptides=TRUE)
#> INFO  [2026-04-25 15:36:35] ** Raw data from Philosopher imported successfully.
#> INFO  [2026-04-25 15:36:35] ** Using provided annotation.
#> INFO  [2026-04-25 15:36:35] ** Run and Channel labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-25 15:36:35] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements within each run will be kept.
#> INFO  [2026-04-25 15:36:35] ** Rows with values not greater than 0.6 in Purity are removed 
#> WARN  [2026-04-25 15:36:35] ** PeptideProphetProbability not found in input columns.
#> INFO  [2026-04-25 15:36:35] ** Sequences containing Oxidation are removed.
#> INFO  [2026-04-25 15:36:35] ** Features with all missing measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:35] ** Shared peptides are removed.
#> INFO  [2026-04-25 15:36:35] ** Features with all missing measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:35] ** PSMs have been aggregated to peptide ions.
#> INFO  [2026-04-25 15:36:35] ** Run annotation merged with quantification data.
#> INFO  [2026-04-25 15:36:35] ** Features with one or two measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:35] ** Fractionation handled.
#> INFO  [2026-04-25 15:36:35] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-25 15:36:35] ** Finished preprocessing. The dataset is ready to be processed by the proteinSummarization function.
#> INFO  [2026-04-25 15:36:35] ** Raw data from Philosopher imported successfully.
#> INFO  [2026-04-25 15:36:35] ** Using provided annotation.
#> INFO  [2026-04-25 15:36:35] ** Run and Channel labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-25 15:36:35] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements within each run will be kept.
#> INFO  [2026-04-25 15:36:35] ** Rows with values not greater than 0.6 in Purity are removed 
#> WARN  [2026-04-25 15:36:35] ** PeptideProphetProbability not found in input columns.
#> INFO  [2026-04-25 15:36:35] ** Sequences containing Oxidation are removed.
#> INFO  [2026-04-25 15:36:35] ** Features with all missing measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:35] ** Shared peptides are removed.
#> INFO  [2026-04-25 15:36:35] ** Features with all missing measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:35] ** PSMs have been aggregated to peptide ions.
#> INFO  [2026-04-25 15:36:35] ** Run annotation merged with quantification data.
#> INFO  [2026-04-25 15:36:35] ** Features with one or two measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:35] ** Fractionation handled.
#> INFO  [2026-04-25 15:36:35] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-25 15:36:35] ** Finished preprocessing. The dataset is ready to be processed by the proteinSummarization function.
head(msstats_data$PTM)
#>                  ProteinName       PeptideSequence Charge
#> 1 sp|Q9Y2W1|TR150_HUMAN_S805 AEEYTEETEEREES*TTGFDK      3
#> 2 sp|Q9Y2W1|TR150_HUMAN_S805 AEEYTEETEEREES*TTGFDK      3
#> 3 sp|Q9Y2W1|TR150_HUMAN_S805 AEEYTEETEEREES*TTGFDK      3
#> 4 sp|Q9Y2W1|TR150_HUMAN_S805 AEEYTEETEEREES*TTGFDK      3
#> 5 sp|Q9Y2W1|TR150_HUMAN_S805 AEEYTEETEEREES*TTGFDK      3
#> 6 sp|Q9Y2W1|TR150_HUMAN_S805 AEEYTEETEEREES*TTGFDK      3
#>                       PSM Mixture TechRepMixture
#> 1 AEEYTEETEEREES*TTGFDK_3  plex16              1
#> 2 AEEYTEETEEREES*TTGFDK_3  plex16              1
#> 3 AEEYTEETEEREES*TTGFDK_3  plex16              1
#> 4 AEEYTEETEEREES*TTGFDK_3  plex16              1
#> 5 AEEYTEETEEREES*TTGFDK_3  plex16              1
#> 6 AEEYTEETEEREES*TTGFDK_3  plex16              1
#>                                      Run Channel    BioReplicate Condition
#> 1 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01     126 CPT0088900003_T         T
#> 2 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01    127C CPT0088920001_N         N
#> 3 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01    127N CPT0079270003_T         T
#> 4 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01    128C CPT0088550004_T         T
#> 5 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01    128N CPT0079300001_N         N
#> 6 16CPTAC_CCRCC_P_JHU_20180326_LUMOS_f01    129C CPT0014450004_T         T
#>   Intensity
#> 1  47545.47
#> 2  45316.41
#> 3  80388.11
#> 4  66856.88
#> 5 118057.66
#> 6 192263.72
head(msstats_data$PROTEIN)
#>             ProteinName      PeptideSequence Charge                    PSM
#> 1 sp|Q9Y2W1|TR150_HUMAN AEEYTEETEEREESTTGFDK      3 AEEYTEETEEREESTTGFDK_3
#> 2 sp|Q9Y2W1|TR150_HUMAN AEEYTEETEEREESTTGFDK      3 AEEYTEETEEREESTTGFDK_3
#> 3 sp|Q9Y2W1|TR150_HUMAN AEEYTEETEEREESTTGFDK      3 AEEYTEETEEREESTTGFDK_3
#> 4 sp|Q9Y2W1|TR150_HUMAN AEEYTEETEEREESTTGFDK      3 AEEYTEETEEREESTTGFDK_3
#> 5 sp|Q9Y2W1|TR150_HUMAN AEEYTEETEEREESTTGFDK      3 AEEYTEETEEREESTTGFDK_3
#> 6 sp|Q9Y2W1|TR150_HUMAN AEEYTEETEEREESTTGFDK      3 AEEYTEETEEREESTTGFDK_3
#>   Mixture TechRepMixture                                    Run Channel
#> 1  plex16              1 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01     126
#> 2  plex16              1 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01    127C
#> 3  plex16              1 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01    127N
#> 4  plex16              1 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01    128C
#> 5  plex16              1 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01    128N
#> 6  plex16              1 16CPTAC_CCRCC_W_JHU_20180322_LUMOS_f01    129C
#>      BioReplicate Condition Intensity
#> 1 CPT0088900003_T         T  97711.92
#> 2 CPT0088920001_N         N 104208.69
#> 3 CPT0079270003_T         T 107628.15
#> 4 CPT0088550004_T         T 194282.36
#> 5 CPT0079300001_N         N 152762.88
#> 6 CPT0014450004_T         T 143260.84

# LFQ Example
input = system.file("tinytest/raw_data/Fragpipe/MSstats.csv", 
                                        package = "MSstatsPTM")
input = data.table::fread(input)
annot = system.file("tinytest/raw_data/Fragpipe/experiment_annotation.tsv", 
                                        package = "MSstatsPTM")
annot = data.table::fread(annot)       
input_protein = system.file("tinytest/raw_data/Fragpipe/msstats_proteome_lf.csv",
                                        package = "MSstatsPTM")                                  
input_protein = data.table::fread(input_protein)

msstats_data = FragPipetoMSstatsPTMFormat(input,
                                          annot,
                                          input_protein = input_protein,
                                          label_type="LF",
                                          mod_id_col = "STY",
                                          localization_cutoff=.75,
                                          protein_id_col = "ProteinName",
                                          peptide_id_col = "PeptideSequence")
#> INFO  [2026-04-25 15:36:35] ** Raw data from FragPipe imported successfully.
#> INFO  [2026-04-25 15:36:35] ** Using annotation extracted from quantification data.
#> INFO  [2026-04-25 15:36:35] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-25 15:36:35] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be kept.
#> INFO  [2026-04-25 15:36:35] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-25 15:36:35] ** Shared peptides are removed.
#> INFO  [2026-04-25 15:36:35] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: sum
#> INFO  [2026-04-25 15:36:35] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-25 15:36:35] ** Run annotation merged with quantification data.
#> INFO  [2026-04-25 15:36:35] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-25 15:36:35] ** Fractionation handled.
#> INFO  [2026-04-25 15:36:35] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-25 15:36:35] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.
#> INFO  [2026-04-25 15:36:35] ** Raw data from FragPipe imported successfully.
#> INFO  [2026-04-25 15:36:35] ** Using annotation extracted from quantification data.
#> INFO  [2026-04-25 15:36:35] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-25 15:36:35] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be kept.
#> INFO  [2026-04-25 15:36:35] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-25 15:36:35] ** Shared peptides are removed.
#> INFO  [2026-04-25 15:36:35] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: sum
#> INFO  [2026-04-25 15:36:35] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-25 15:36:35] ** Run annotation merged with quantification data.
#> INFO  [2026-04-25 15:36:35] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-25 15:36:35] ** Fractionation handled.
#> INFO  [2026-04-25 15:36:35] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-25 15:36:35] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.
                                          
# If no global profiling run is available, omit input_protein and set:
# msstats_data = FragPipetoMSstatsPTMFormat(input, annot,
#   label_type = "LF", mod_id_col = "STY",
#   localization_cutoff = .75, protein_id_col = "ProteinName",
#   peptide_id_col = "PeptideSequence", use_unmod_peptides = FALSE)

head(msstats_data$PTM)
#>                 ProteinName                                   PeptideSequence
#> 1 sp|P02400|RLA4_YEAST_S100 FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEES*DDDMGFGLFD
#> 2 sp|P02400|RLA4_YEAST_S100 FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEES*DDDMGFGLFD
#> 3 sp|P02400|RLA4_YEAST_S100 FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEES*DDDMGFGLFD
#> 4 sp|P02400|RLA4_YEAST_S100 FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEES*DDDMGFGLFD
#> 5 sp|P02400|RLA4_YEAST_S100 FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEES*DDDMGFGLFD
#> 6 sp|P02400|RLA4_YEAST_S100 FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEES*DDDMGFGLFD
#>   PrecursorCharge FragmentIon ProductCharge IsotopeLabelType Condition
#> 1               3        <NA>            NA                L        WT
#> 2               3        <NA>            NA                L       MUT
#> 3               3        <NA>            NA                L        WT
#> 4               3        <NA>            NA                L       MUT
#> 5               4        <NA>            NA                L        WT
#> 6               4        <NA>            NA                L       MUT
#>   BioReplicate                             Run Fraction Intensity
#> 1            3  JCI-TiO2-DDAphosMethod-WT_1338        1 143712.34
#> 2            2 JCI-TiO2-DDAphosMethod-rho_1339        1        NA
#> 3            4         JCI-TiO2-DDAstd-WT_1335        1 128482.80
#> 4            3        JCI-TiO2-DDAstd-rho_1336        1        NA
#> 5            3  JCI-TiO2-DDAphosMethod-WT_1338        1  64697.76
#> 6            2 JCI-TiO2-DDAphosMethod-rho_1339        1  27288.38
head(msstats_data$PROTEIN)
#>            ProteinName                                  PeptideSequence
#> 1 sp|P02400|RLA4_YEAST FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEESDDDMGFGLFD
#> 2 sp|P02400|RLA4_YEAST FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEESDDDMGFGLFD
#> 3 sp|P02400|RLA4_YEAST FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEESDDDMGFGLFD
#> 4 sp|P02400|RLA4_YEAST FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEESDDDMGFGLFD
#> 5 sp|P02400|RLA4_YEAST FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEESDDDMGFGLFD
#> 6 sp|P02400|RLA4_YEAST FATVPTGGASSAAAGAAGAAAGGDAAEEEKEEEAKEESDDDMGFGLFD
#>   PrecursorCharge FragmentIon ProductCharge IsotopeLabelType Condition
#> 1               3        <NA>            NA                L        WT
#> 2               3        <NA>            NA                L        WT
#> 3               3        <NA>            NA                L       MUT
#> 4               3        <NA>            NA                L        WT
#> 5               3        <NA>            NA                L       MUT
#> 6               3        <NA>            NA                L        WT
#>   BioReplicate                                        Run Fraction Intensity
#> 1            1 3-9-2023_JCI_230308_DDA_TiO2-10L_S1-A5_199        1       0.0
#> 2            2                            CR_Tio2-WT_1715        1       0.0
#> 3            1                       CR_Tio2-phos85D_1716        1       0.0
#> 4            3             JCI-TiO2-DDAphosMethod-WT_1338        1  213712.3
#> 5            2            JCI-TiO2-DDAphosMethod-rho_1339        1  141983.5
#> 6            4                    JCI-TiO2-DDAstd-WT_1335        1  235076.7
```
