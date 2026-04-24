# Convert the output of DIA-NN PSM file into MSstatsPTM format

Takes as input the `report.tsv` file from DIA-NN and converts it into
MSstatsPTM format. Requires PSM and an annotation file. Optionally an
additional `report.tsv` file for a corresponding global profiling run
can be included.

## Usage

``` r
DIANNtoMSstatsPTMFormat(
  input,
  annotation,
  input_protein = NULL,
  annotation_protein = NULL,
  fasta_path = NULL,
  use_unmod_peptides = FALSE,
  protein_id_col = "Protein.Group",
  fasta_protein_name = "uniprot_ac",
  global_qvalue_cutoff = 0.01,
  qvalue_cutoff = 0.01,
  pg_qvalue_cutoff = 0.01,
  useUniquePeptide = TRUE,
  removeFewMeasurements = TRUE,
  removeOxidationMpeptides = TRUE,
  removeProtein_with1Feature = FALSE,
  MBR = TRUE,
  quantificationColumn = "FragmentQuantCorrected",
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = NULL
)
```

## Arguments

- input:

  name of MSstats input report from Diann, which includes fragment-level
  data. Output fragment data with –export-quant flag in DIA-NN 2.0

- annotation:

  name of 'annotation.txt' data which includes Condition, BioReplicate,
  Run.

- input_protein:

  same as `input` for global profiling run. Default is NULL.

- annotation_protein:

  same as `annotation` for global profiling run. Default is NULL.

- fasta_path:

  A string of path to a FASTA file, used to match PTM peptides.

- use_unmod_peptides:

  Boolean if the unmodified peptides in the input file should be used to
  construct the unmodified protein output. Only used if `input_protein`
  is not provided. Default is `FALSE`.

- protein_id_col:

  Use 'Protein.Groups'(default) column for protein name.

- fasta_protein_name:

  Name of column that matches with the protein names in
  `protein_id_col`. The protein names in these two columns must match in
  order to join the FASTA file with the DIA-NN output. Default is
  "uniprot_ac" for uniprot ID. For uniprot mnemonic ID, use "entry_name"

- global_qvalue_cutoff:

  The qvalue cutoff for the Q.Value column, i.e. the run-specific
  precursor q-value. Default is 0.01.

- qvalue_cutoff:

  If MBR is false, the qvalue cutoff for the Global.Q.Value column, i.e.
  global precursor q-value. If MBR is true, the qvalue cutoff for the
  Lib.Q.Value column, i.e. the q-value for the library created after the
  first MBR pass. Default is 0.01.

- pg_qvalue_cutoff:

  If MBR is false, the qvalue cutoff for the Global.PG.Q.Value column,
  i.e. the global q-value for the protein group. If MBR is true, the
  qvalue cutoff for the Lib.PG.Q.Value column, i.e. the protein group
  q-value for the library created after the first MBR pass. Default is
  0.01.

- useUniquePeptide:

  should unique peptides be removed

- removeFewMeasurements:

  should proteins with few measurements be removed

- removeOxidationMpeptides:

  should peptides with oxidation be removed

- removeProtein_with1Feature:

  should proteins with a single feature be removed

- MBR:

  True if analysis was done with match between runs

- quantificationColumn:

  Use 'FragmentQuantCorrected'(default) column for quantified
  intensities for DIANN 1.8.x. Use 'FragmentQuantRaw' for quantified
  intensities for DIANN 1.9.x. Use 'auto' for quantified intensities for
  DIANN 2.x where each fragment intensity is a separate column, e.g.
  Fr0Quantity.

- use_log_file:

  logical. If TRUE, information about data processing will be saved to a
  file.

- append:

  logical. If TRUE, information about data processing will be added to
  an existing log file.

- verbose:

  logical. If TRUE, information about data processing will be printed to
  the console.

- log_file_path:

  character. Path to a file to which information about data processing
  will be saved. If not provided, such a file will be created
  automatically. If `append = TRUE`, has to be a valid path to a file.

## Value

`list` of one or two `data.frame` of class `MSstatsTMT`, named `PTM` and
`PROTEIN`

## Examples

``` r
# Example from PRIDE ID PXD053502
input = system.file("tinytest/raw_data/DIANN/report.tsv", 
                                        package = "MSstatsPTM")
input = data.table::fread(input)
annot = system.file("tinytest/raw_data/DIANN/annot.csv", 
                                        package = "MSstatsPTM")
annot = data.table::fread(annot)
fasta_path = system.file("extdata", "diann.fasta", 
                       package="MSstatsPTM")

msstatsptm_format = DIANNtoMSstatsPTMFormat(
    input, 
    annot, 
    protein_id_col = "Protein.Names", 
    fasta_path = fasta_path, 
    fasta_protein_name = "entry_name", 
    use_log_file = FALSE
)
#> INFO  [2026-04-24 21:03:26] ** Raw data from DIANN imported successfully.
#> INFO  [2026-04-24 21:03:26] ** Raw data from DIANN cleaned successfully.
#> INFO  [2026-04-24 21:03:26] ** Using provided annotation.
#> INFO  [2026-04-24 21:03:26] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-24 21:03:26] ** Filtering on Global Q Value < 0.01
#> INFO  [2026-04-24 21:03:26] ** MBR was used to analyze the data. Now setting names and filtering
#> INFO  [2026-04-24 21:03:26] -- LibPGQValue < 0.01
#> INFO  [2026-04-24 21:03:26] -- LibQValue < 0.01
#> INFO  [2026-04-24 21:03:26] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be removed.
#> INFO  [2026-04-24 21:03:26] ** Sequences containing DECOY, Decoys are removed.
#> INFO  [2026-04-24 21:03:26] ** Sequences containing \(UniMod\:35\) are removed.
#> INFO  [2026-04-24 21:03:26] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-24 21:03:26] ** Shared peptides are removed.
#> INFO  [2026-04-24 21:03:26] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: max
#> INFO  [2026-04-24 21:03:26] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-04-24 21:03:26] ** Run annotation merged with quantification data.
#> WARN  [2026-04-24 21:03:26] The following features have missing values in at least one run. ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK_3_Frag3_1,
#>  ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK_3_Frag5_1,
#>  ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK_3_Frag10_1,
#>  ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK_3_Frag12_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_3_Frag1_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_3_Frag2_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_3_Frag3_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_3_Frag4_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_3_Frag5_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_3_Frag6_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_3_Frag9_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag1_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag2_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag3_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag4_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag5_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag6_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag7_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag8_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag9_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag10_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag11_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag12_1,
#>  DSNPEEIEIDFETLK(UniMod:121)PSTLR_3_Frag1_1,
#>  DSNPEEIEIDFETLK(UniMod:121)PSTLR_3_Frag8_1,
#>  DSNPEEIEIDFETLK(UniMod:121)PSTLR_3_Frag12_1,
#>  HPMDLSTVK(UniMod:121)R_3_Frag1_1,
#>  HPMDLSTVK(UniMod:121)R_3_Frag2_1,
#>  HPMDLSTVK(UniMod:121)R_3_Frag3_1,
#>  HPMDLSTVK(UniMod:121)R_3_Frag4_1,
#>  HPMDLSTVK(UniMod:121)R_3_Frag5_1,
#>  HPMDLSTVK(UniMod:121)R_3_Frag6_1,
#>  HPMDLSTVK(UniMod:121)R_3_Frag7_1,
#>  HPMDLSTVK(UniMod:121)R_3_Frag8_1,
#>  HPMDLSTVK(UniMod:121)R_3_Frag10_1,
#>  HPMDLSTVK(UniMod:121)R_3_Frag12_1,
#>  K(UniMod:121)LQDVFEFR_2_Frag1_1,
#>  K(UniMod:121)LQDVFEFR_2_Frag2_1,
#>  K(UniMod:121)LQDVFEFR_2_Frag3_1,
#>  K(UniMod:121)LQDVFEFR_2_Frag4_1,
#>  K(UniMod:121)LQDVFEFR_2_Frag5_1,
#>  K(UniMod:121)LQDVFEFR_2_Frag6_1,
#>  K(UniMod:121)LQDVFEFR_2_Frag7_1,
#>  K(UniMod:121)LQDVFEFR_2_Frag8_1,
#>  K(UniMod:121)LQDVFEFR_3_Frag1_1,
#>  K(UniMod:121)LQDVFEFR_3_Frag2_1,
#>  K(UniMod:121)LQDVFEFR_3_Frag3_1,
#>  K(UniMod:121)LQDVFEFR_3_Frag4_1,
#>  K(UniMod:121)LQDVFEFR_3_Frag5_1,
#>  K(UniMod:121)LQDVFEFR_3_Frag6_1,
#>  K(UniMod:121)LQDVFEFR_3_Frag12_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHDVVAMAR_4_Frag1_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHDVVAMAR_4_Frag4_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHDVVAMAR_4_Frag5_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHDVVAMAR_4_Frag6_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHDVVAMAR_4_Frag7_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHDVVAMAR_4_Frag10_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHDVVAMAR_4_Frag12_1,
#>  QPMDMGTIK(UniMod:121)R_3_Frag1_1,
#>  QPMDMGTIK(UniMod:121)R_3_Frag2_1,
#>  QPMDMGTIK(UniMod:121)R_3_Frag3_1,
#>  QPMDMGTIK(UniMod:121)R_3_Frag4_1,
#>  QPMDMGTIK(UniMod:121)R_3_Frag6_1,
#>  QPMDMGTIK(UniMod:121)R_3_Frag7_1,
#>  QPMDMGTIK(UniMod:121)R_3_Frag8_1,
#>  QPMDMGTIK(UniMod:121)R_3_Frag10_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_3_Frag1_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_3_Frag2_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_3_Frag3_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_3_Frag4_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_3_Frag5_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_3_Frag6_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_3_Frag11_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_4_Frag1_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_4_Frag2_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_4_Frag3_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_4_Frag4_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_4_Frag5_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_4_Frag6_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_4_Frag7_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_4_Frag8_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_4_Frag9_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_4_Frag10_1,
#>  QPVDAVK(UniMod:121)LGLPDYHK_4_Frag12_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag1_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag2_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag3_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag4_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag5_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag6_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag7_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag8_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag9_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag10_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag11_1,
#>  SLHSAGPPLLAVTAAPPAQPLAK(UniMod:121)K_4_Frag12_1,
#>  TK(UniMod:121)EELALEK_2_Frag1_1,
#>  TK(UniMod:121)EELALEK_2_Frag2_1,
#>  TK(UniMod:121)EELALEK_2_Frag3_1,
#>  TK(UniMod:121)EELALEK_2_Frag4_1,
#>  TK(UniMod:121)EELALEK_2_Frag6_1,
#>  TK(UniMod:121)EELALEK_2_Frag7_1,
#>  TK(UniMod:121)EELALEK_2_Frag8_1,
#>  TK(UniMod:121)EELALEK_2_Frag9_1,
#>  TK(UniMod:121)EELALEK_2_Frag10_1,
#>  TK(UniMod:121)EELALEK_2_Frag12_1,
#>  TK(UniMod:121)EELALEK_3_Frag1_1,
#>  TK(UniMod:121)EELALEK_3_Frag3_1,
#>  TK(UniMod:121)EELALEK_3_Frag4_1,
#>  TK(UniMod:121)EELALEK_3_Frag5_1,
#>  TK(UniMod:121)EELALEK_3_Frag6_1,
#>  TK(UniMod:121)EELALEK_3_Frag8_1,
#>  AVHEQLAALSQAPVNK(UniMod:121)PK_4_Frag1_1,
#>  AVHEQLAALSQAPVNK(UniMod:121)PK_4_Frag2_1,
#>  AVHEQLAALSQAPVNK(UniMod:121)PK_4_Frag4_1,
#>  AVHEQLAALSQAPVNK(UniMod:121)PK_4_Frag6_1,
#>  AVHEQLAALSQAPVNK(UniMod:121)PK_4_Frag7_1,
#>  AVHEQLAALSQAPVNK(UniMod:121)PK_4_Frag9_1,
#>  AVHEQLAALSQAPVNK(UniMod:121)PK_4_Frag11_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_3_Frag1_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_3_Frag2_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_3_Frag3_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_3_Frag4_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_3_Frag5_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_3_Frag6_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_3_Frag7_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_3_Frag9_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_3_Frag11_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_4_Frag1_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_4_Frag2_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_4_Frag3_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_4_Frag4_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_4_Frag5_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_4_Frag6_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_4_Frag7_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_4_Frag9_1,
#>  DLEDGEVPQHAGK(UniMod:121)K_4_Frag10_1,
#>  FAK(UniMod:121)MPDEPVEAPALPAPAAPMVSK_3_Frag1_1,
#>  FAK(UniMod:121)MPDEPVEAPALPAPAAPMVSK_3_Frag2_1,
#>  FAK(UniMod:121)MPDEPVEAPALPAPAAPMVSK_3_Frag3_1,
#>  FAK(UniMod:121)MPDEPVEAPALPAPAAPMVSK_3_Frag4_1,
#>  FAK(UniMod:121)MPDEPVEAPALPAPAAPMVSK_3_Frag5_1,
#>  FAK(UniMod:121)MPDEPVEAPALPAPAAPMVSK_3_Frag6_1,
#>  FAK(UniMod:121)MPDEPVEAPALPAPAAPMVSK_3_Frag7_1,
#>  FAK(UniMod:121)MPDEPVEAPALPAPAAPMVSK_3_Frag8_1,
#>  FAK(UniMod:121)MPDEPVEAPALPAPAAPMVSK_3_Frag9_1,
#>  FAK(UniMod:121)MPDEPVEAPALPAPAAPMVSK_3_Frag10_1,
#>  FAK(UniMod:121)MPDEPVEAPALPAPAAPMVSK_3_Frag11_1,
#>  LNLPDYHK(UniMod:121)IIK_3_Frag1_1,
#>  LNLPDYHK(UniMod:121)IIK_3_Frag2_1,
#>  LNLPDYHK(UniMod:121)IIK_3_Frag3_1,
#>  LNLPDYHK(UniMod:121)IIK_3_Frag4_1,
#>  LNLPDYHK(UniMod:121)IIK_3_Frag6_1,
#>  LNLPDYHK(UniMod:121)IIK_3_Frag7_1,
#>  LNLPDYHK(UniMod:121)IIK_3_Frag9_1,
#>  LNLPDYHK(UniMod:121)IIK_3_Frag11_1,
#>  LNLPDYHK(UniMod:121)IIK_3_Frag12_1,
#>  LQDVSGQLSSSK(UniMod:121)K_3_Frag1_1,
#>  LQDVSGQLSSSK(UniMod:121)K_3_Frag2_1,
#>  LQDVSGQLSSSK(UniMod:121)K_3_Frag3_1,
#>  LQDVSGQLSSSK(UniMod:121)K_3_Frag4_1,
#>  LQDVSGQLSSSK(UniMod:121)K_3_Frag5_1,
#>  LQDVSGQLSSSK(UniMod:121)K_3_Frag6_1,
#>  LQDVSGQLSSSK(UniMod:121)K_3_Frag8_1,
#>  LQDVSGQLSSSK(UniMod:121)K_3_Frag10_1,
#>  LQDVSGQLSSSK(UniMod:121)K_3_Frag11_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_3_Frag1_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_3_Frag2_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_3_Frag3_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_3_Frag4_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_3_Frag5_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_3_Frag6_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_3_Frag12_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_4_Frag1_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_4_Frag2_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_4_Frag3_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_4_Frag4_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_4_Frag5_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_4_Frag6_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_4_Frag8_1,
#>  QLSLDINRLPGEK(UniMod:121)LGR_4_Frag11_1,
#>  VAQMPQEEVELLPPAPK(UniMod:121)GK_3_Frag1_1,
#>  VAQMPQEEVELLPPAPK(UniMod:121)GK_3_Frag3_1,
#>  VAQMPQEEVELLPPAPK(UniMod:121)GK_3_Frag4_1,
#>  VAQMPQEEVELLPPAPK(UniMod:121)GK_3_Frag5_1,
#>  VAQMPQEEVELLPPAPK(UniMod:121)GK_3_Frag6_1,
#>  VAQMPQEEVELLPPAPK(UniMod:121)GK_3_Frag8_1,
#>  VAQMPQEEVELLPPAPK(UniMod:121)GK_3_Frag12_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_2_Frag1_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_2_Frag2_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_2_Frag3_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_2_Frag4_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_2_Frag5_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_2_Frag7_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_2_Frag8_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_2_Frag9_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_2_Frag10_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_2_Frag11_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag1_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag2_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag3_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag4_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag5_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag6_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag7_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag8_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag9_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag10_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag11_1,
#>  C(UniMod:4)C(UniMod:4)SGILK(UniMod:121)EMFAK_3_Frag12_1,
#>  DVPDSQQHPAPEK(UniMod:121)SSK_3_Frag1_1,
#>  DVPDSQQHPAPEK(UniMod:121)SSK_3_Frag4_1,
#>  DVPDSQQHPAPEK(UniMod:121)SSK_3_Frag5_1,
#>  DVPDSQQHPAPEK(UniMod:121)SSK_3_Frag9_1,
#>  DVPDSQQHPAPEK(UniMod:121)SSK_3_Frag11_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag1_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag2_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag3_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag4_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag5_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag6_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag7_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag8_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag9_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag10_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag11_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_3_Frag12_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_4_Frag1_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_4_Frag2_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_4_Frag3_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_4_Frag4_1,
#>  FAK(UniMod:121)MPDEPEEPVVAVSSPAVPPPTK_4_Frag6_1,
#>  HPMDMSTIK(UniMod:121)SK_3_Frag1_1,
#>  HPMDMSTIK(UniMod:121)SK_3_Frag2_1,
#>  HPMDMSTIK(UniMod:121)SK_3_Frag3_1,
#>  HPMDMSTIK(UniMod:121)SK_3_Frag4_1,
#>  HPMDMSTIK(UniMod:121)SK_3_Frag5_1,
#>  HPMDMSTIK(UniMod:121)SK_3_Frag6_1,
#>  INELPTEETEIMIVQAK(UniMod:121)GR_3_Frag1_1,
#>  INELPTEETEIMIVQAK(UniMod:121)GR_3_Frag2_1,
#>  INELPTEETEIMIVQAK(UniMod:121)GR_3_Frag3_1,
#>  INELPTEETEIMIVQAK(UniMod:121)GR_3_Frag4_1,
#>  INELPTEETEIMIVQAK(UniMod:121)GR_3_Frag5_1,
#>  INELPTEETEIMIVQAK(UniMod:121)GR_3_Frag6_1,
#>  INELPTEETEIMIVQAK(UniMod:121)GR_3_Frag7_1,
#>  INELPTEETEIMIVQAK(UniMod:121)GR_3_Frag8_1,
#>  INELPTEETEIMIVQAK(UniMod:121)GR_3_Frag9_1,
#>  INELPTEETEIMIVQAK(UniMod:121)GR_3_Frag12_1,
#>  K(UniMod:121)LQDVFEMR_2_Frag1_1,
#>  K(UniMod:121)LQDVFEMR_2_Frag2_1,
#>  K(UniMod:121)LQDVFEMR_2_Frag3_1,
#>  K(UniMod:121)LQDVFEMR_2_Frag4_1,
#>  K(UniMod:121)LQDVFEMR_2_Frag5_1,
#>  K(UniMod:121)LQDVFEMR_2_Frag6_1,
#>  K(UniMod:121)LQDVFEMR_2_Frag7_1,
#>  K(UniMod:121)LQDVFEMR_2_Frag8_1,
#>  K(UniMod:121)LQDVFEMR_2_Frag9_1,
#>  K(UniMod:121)LQDVFEMR_2_Frag10_1,
#>  K(UniMod:121)LQDVFEMR_3_Frag1_1,
#>  K(UniMod:121)LQDVFEMR_3_Frag2_1,
#>  K(UniMod:121)LQDVFEMR_3_Frag3_1,
#>  K(UniMod:121)LQDVFEMR_3_Frag4_1,
#>  K(UniMod:121)LQDVFEMR_3_Frag5_1,
#>  K(UniMod:121)LQDVFEMR_3_Frag6_1,
#>  K(UniMod:121)LQDVFEMR_3_Frag7_1,
#>  K(UniMod:121)LQDVFEMR_3_Frag9_1,
#>  K(UniMod:121)LQDVFEMR_3_Frag10_1,
#>  K(UniMod:121)LQDVFEMR_3_Frag12_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_3_Frag1_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_3_Frag2_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_3_Frag3_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_3_Frag6_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_3_Frag9_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_3_Frag10_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_3_Frag12_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag1_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag2_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag3_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag4_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag5_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag6_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag7_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag8_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag9_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag10_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag11_1,
#>  LMFSNC(UniMod:4)YK(UniMod:121)YNPPDHEVVAMAR_4_Frag12_1,
#>  LNLPDYYK(UniMod:121)IIK_2_Frag1_1,
#>  LNLPDYYK(UniMod:121)IIK_2_Frag3_1,
#>  LNLPDYYK(UniMod:121)IIK_2_Frag4_1,
#>  LNLPDYYK(UniMod:121)IIK_2_Frag6_1,
#>  LNLPDYYK(UniMod:121)IIK_2_Frag7_1,
#>  LNLPDYYK(UniMod:121)IIK_2_Frag8_1,
#>  LNLPDYYK(UniMod:121)IIK_2_Frag9_1,
#>  LNLPDYYK(UniMod:121)IIK_2_Frag10_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag1_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag2_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag3_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag4_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag5_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag6_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag7_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag8_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag9_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag10_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag11_1,
#>  LNLPDYYK(UniMod:121)IIK_3_Frag12_1,
#>  LPGEK(UniMod:121)LGR_2_Frag1_1,
#>  LPGEK(UniMod:121)LGR_2_Frag2_1,
#>  LPGEK(UniMod:121)LGR_2_Frag3_1,
#>  LPGEK(UniMod:121)LGR_2_Frag4_1,
#>  LPGEK(UniMod:121)LGR_2_Frag5_1,
#>  LPGEK(UniMod:121)LGR_2_Frag6_1,
#>  LPGEK(UniMod:121)LGR_2_Frag7_1,
#>  LPGEK(UniMod:121)LGR_2_Frag8_1,
#>  LPGEK(UniMod:121)LGR_2_Frag9_1,
#>  LPGEK(UniMod:121)LGR_2_Frag10_1,
#>  LPGEK(UniMod:121)LGR_2_Frag11_1,
#>  NSNPDEIEIDFETLK(UniMod:121)PSTLR_3_Frag1_1,
#>  NSNPDEIEIDFETLK(UniMod:121)PSTLR_3_Frag2_1,
#>  NSNPDEIEIDFETLK(UniMod:121)PSTLR_3_Frag3_1,
#>  NSNPDEIEIDFETLK(UniMod:121)PSTLR_3_Frag4_1,
#>  NSNPDEIEIDFETLK(UniMod:121)PSTLR_3_Frag6_1,
#>  NSNPDEIEIDFETLK(UniMod:121)PSTLR_3_Frag8_1,
#>  NSNPDEIEIDFETLK(UniMod:121)PSTLR_3_Frag10_1,
#>  QLSLDINK(UniMod:121)LPGEK_3_Frag1_1,
#>  QLSLDINK(UniMod:121)LPGEK_3_Frag2_1,
#>  QLSLDINK(UniMod:121)LPGEK_3_Frag3_1,
#>  QLSLDINK(UniMod:121)LPGEK_3_Frag4_1,
#>  QLSLDINK(UniMod:121)LPGEK_3_Frag5_1,
#>  QLSLDINK(UniMod:121)LPGEK_3_Frag6_1,
#>  QLSLDINK(UniMod:121)LPGEK_3_Frag8_1,
#>  QLSLDINK(UniMod:121)LPGEK_3_Frag9_1,
#>  QLSLDINK(UniMod:121)LPGEK_3_Frag10_1,
#>  QLSLDINK(UniMod:121)LPGEK_3_Frag12_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_3_Frag1_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_3_Frag2_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_3_Frag3_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_3_Frag4_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_3_Frag5_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_3_Frag6_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_3_Frag7_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_3_Frag8_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_3_Frag9_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_3_Frag10_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_3_Frag12_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag1_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag2_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag3_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag4_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag5_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag6_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag7_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag8_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag9_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag10_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag11_1,
#>  QLSLDINKLPGEK(UniMod:121)LGR_4_Frag12_1,
#>  SSK(UniMod:121)VSEQLK_2_Frag1_1,
#>  SSK(UniMod:121)VSEQLK_2_Frag2_1,
#>  SSK(UniMod:121)VSEQLK_2_Frag3_1,
#>  SSK(UniMod:121)VSEQLK_2_Frag8_1,
#>  TPMDMGTIK(UniMod:121)K_2_Frag1_1,
#>  TPMDMGTIK(UniMod:121)K_2_Frag2_1,
#>  TPMDMGTIK(UniMod:121)K_2_Frag3_1,
#>  TPMDMGTIK(UniMod:121)K_2_Frag4_1,
#>  TPMDMGTIK(UniMod:121)K_2_Frag5_1,
#>  TPMDMGTIK(UniMod:121)K_2_Frag6_1,
#>  TPMDMGTIK(UniMod:121)K_2_Frag7_1,
#>  TPMDMGTIK(UniMod:121)K_2_Frag9_1,
#>  TPMDMGTIK(UniMod:121)K_3_Frag1_1,
#>  TPMDMGTIK(UniMod:121)K_3_Frag2_1,
#>  TPMDMGTIK(UniMod:121)K_3_Frag3_1,
#>  TPMDMGTIK(UniMod:121)K_3_Frag4_1,
#>  TPMDMGTIK(UniMod:121)K_3_Frag5_1,
#>  TPMDMGTIK(UniMod:121)K_3_Frag6_1,
#>  TPMDMGTIK(UniMod:121)K_3_Frag7_1,
#>  TPMDMGTIK(UniMod:121)K_3_Frag8_1,
#>  TPMDMGTIK(UniMod:121)K_3_Frag10_1,
#>  TPMDMGTIK(UniMod:121)K_3_Frag11_1,
#>  VDVIAGSSK(UniMod:121)MK_2_Frag1_1,
#>  VDVIAGSSK(UniMod:121)MK_2_Frag2_1,
#>  VDVIAGSSK(UniMod:121)MK_2_Frag3_1,
#>  VDVIAGSSK(UniMod:121)MK_2_Frag4_1,
#>  VDVIAGSSK(UniMod:121)MK_2_Frag9_1,
#>  VSEQLK(UniMod:121)C(UniMod:4)C(UniMod:4)SGILK_3_Frag1_1,
#>  VSEQLK(UniMod:121)C(UniMod:4)C(UniMod:4)SGILK_3_Frag2_1,
#>  VSEQLK(UniMod:121)C(UniMod:4)C(UniMod:4)SGILK_3_Frag3_1,
#>  VSEQLK(UniMod:121)C(UniMod:4)C(UniMod:4)SGILK_3_Frag4_1,
#>  VSEQLK(UniMod:121)C(UniMod:4)C(UniMod:4)SGILK_3_Frag5_1,
#>  VSEQLK(UniMod:121)C(UniMod:4)C(UniMod:4)SGILK_3_Frag6_1,
#>  VSEQLK(UniMod:121)C(UniMod:4)C(UniMod:4)SGILK_3_Frag7_1,
#>  VSEQLK(UniMod:121)C(UniMod:4)C(UniMod:4)SGILK_3_Frag8_1,
#>  VSEQLK(UniMod:121)C(UniMod:4)C(UniMod:4)SGILK_3_Frag9_1,
#>  VSEQLK(UniMod:121)C(UniMod:4)C(UniMod:4)SGILK_3_Frag10_1,
#>  VSEQLK(UniMod:121)C(UniMod:4)C(UniMod:4)SGILK_3_Frag12_1,
#>  VVLK(UniMod:121)TLWK_2_Frag1_1,
#>  VVLK(UniMod:121)TLWK_2_Frag2_1,
#>  VVLK(UniMod:121)TLWK_2_Frag3_1,
#>  VVLK(UniMod:121)TLWK_2_Frag4_1,
#>  VVLK(UniMod:121)TLWK_2_Frag5_1,
#>  VVLK(UniMod:121)TLWK_2_Frag7_1
#> INFO  [2026-04-24 21:03:26] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-24 21:03:26] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.

head(msstatsptm_format$PTM)
#>       ProteinName                          PeptideSequence PrecursorCharge
#> 1 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 2 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 3 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 4 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 5 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 6 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#>   FragmentIon ProductCharge IsotopeLabelType Condition BioReplicate
#> 1       Frag3             1            Light      DMSO            2
#> 2       Frag5             1            Light      DMSO            2
#> 3      Frag10             1            Light      DMSO            2
#> 4      Frag12             1            Light      DMSO            2
#> 5       Frag3             1            Light      MZ-1            4
#> 6       Frag5             1            Light      MZ-1            4
#>                   Run Fraction Intensity
#> 1 144-2024-GS-DMSO-R2        1   973.400
#> 2 144-2024-GS-DMSO-R2        1   408.661
#> 3 144-2024-GS-DMSO-R2        1  1259.200
#> 4 144-2024-GS-DMSO-R2        1   383.783
#> 5  144-2024-GS-MZ1-R1        1  1735.130
#> 6  144-2024-GS-MZ1-R1        1   686.085

# Example DIANN 2.0
input = system.file("tinytest/raw_data/DIANN/diann_2_ptm.parquet", 
                                        package = "MSstatsPTM")
input = arrow::read_parquet(input)
annot = system.file("tinytest/raw_data/DIANN/annotation_diann_2.0_ptm.csv", 
                                        package = "MSstatsPTM")
annot = data.table::fread(annot)
fasta_path = system.file("extdata", "diann.fasta", 
                       package="MSstatsPTM")

msstatsptm_format = DIANNtoMSstatsPTMFormat(
    input, 
    annot, 
    protein_id_col = "Protein.Names", 
    fasta_path = fasta_path, 
    fasta_protein_name = "entry_name", 
    use_log_file = FALSE,
    quantificationColumn = "auto"
)
#> INFO  [2026-04-24 21:03:26] ** Raw data from DIANN imported successfully.
#> INFO  [2026-04-24 21:03:26] ** Raw data from DIANN cleaned successfully.
#> INFO  [2026-04-24 21:03:26] ** Using provided annotation.
#> INFO  [2026-04-24 21:03:26] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-24 21:03:26] ** Filtering on Global Q Value < 0.01
#> INFO  [2026-04-24 21:03:26] ** MBR was used to analyze the data. Now setting names and filtering
#> INFO  [2026-04-24 21:03:26] -- LibPGQValue < 0.01
#> INFO  [2026-04-24 21:03:26] -- LibQValue < 0.01
#> INFO  [2026-04-24 21:03:26] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be removed.
#> INFO  [2026-04-24 21:03:26] ** Sequences containing DECOY, Decoys are removed.
#> INFO  [2026-04-24 21:03:26] ** Sequences containing \(UniMod\:35\) are removed.
#> INFO  [2026-04-24 21:03:26] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-24 21:03:26] ** Shared peptides are removed.
#> INFO  [2026-04-24 21:03:26] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: max
#> INFO  [2026-04-24 21:03:26] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-04-24 21:03:26] ** Run annotation merged with quantification data.
#> WARN  [2026-04-24 21:03:26] The following features have missing values in at least one run. AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag1_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag2_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag3_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag4_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag5_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag6_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag7_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag8_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag9_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag10_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag11_1,
#>  AVHEQLAALSQGPISK(UniMod:121)PK_4_Frag12_1
#> INFO  [2026-04-24 21:03:26] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-24 21:03:26] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.

head(msstatsptm_format$PTM)
#>       ProteinName                          PeptideSequence PrecursorCharge
#> 1 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 2 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 3 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 4 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 5 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 6 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#>   FragmentIon ProductCharge IsotopeLabelType Condition BioReplicate  Run
#> 1       Frag1             1            Light   Control            1 Run1
#> 2       Frag2             1            Light   Control            1 Run1
#> 3       Frag3             1            Light   Control            1 Run1
#> 4       Frag4             1            Light   Control            1 Run1
#> 5       Frag5             1            Light   Control            1 Run1
#> 6       Frag6             1            Light   Control            1 Run1
#>   Fraction  Intensity
#> 1        1 102640.406
#> 2        1  20122.303
#> 3        1   4221.658
#> 4        1  31232.246
#> 5        1   5321.414
#> 6        1   7399.935
```
