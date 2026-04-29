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
#> INFO  [2026-04-29 18:52:22] ** Raw data from DIANN imported successfully.
#> INFO  [2026-04-29 18:52:22] ** Filtering on Q.Value < 0.01
#> INFO  [2026-04-29 18:52:22] ** MBR was used to analyze the data. Now setting names and filtering
#> INFO  [2026-04-29 18:52:22] -- LibPGQValue < 0.01
#> INFO  [2026-04-29 18:52:22] -- LibQValue < 0.01
#> INFO  [2026-04-29 18:52:22] ** Raw data from DIANN cleaned successfully.
#> INFO  [2026-04-29 18:52:22] ** Using provided annotation.
#> INFO  [2026-04-29 18:52:22] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-29 18:52:22] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be removed.
#> INFO  [2026-04-29 18:52:22] ** Sequences containing DECOY, Decoys are removed.
#> INFO  [2026-04-29 18:52:22] ** Sequences containing \(UniMod\:35\) are removed.
#> INFO  [2026-04-29 18:52:22] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-29 18:52:22] ** Shared peptides are removed.
#> INFO  [2026-04-29 18:52:22] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: max
#> INFO  [2026-04-29 18:52:22] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-04-29 18:52:22] ** Run annotation merged with quantification data.
#> INFO  [2026-04-29 18:52:22] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-29 18:52:22] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.

head(msstatsptm_format$PTM)
#>       ProteinName                          PeptideSequence PrecursorCharge
#> 1 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 2 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 3 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 4 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 5 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 6 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#>   FragmentIon ProductCharge IsotopeLabelType Condition BioReplicate
#> 1      Frag10             1            Light      DMSO            1
#> 2      Frag10             1            Light      DMSO            2
#> 3      Frag10             1            Light      DMSO            3
#> 4      Frag10             1            Light      MZ-1            4
#> 5      Frag10             1            Light      MZ-1            5
#> 6      Frag10             1            Light      MZ-1            6
#>                   Run Fraction Intensity
#> 1 144-2024-GS-DMSO-R1        1        NA
#> 2 144-2024-GS-DMSO-R2        1  1259.200
#> 3 144-2024-GS-DMSO-R3        1        NA
#> 4  144-2024-GS-MZ1-R1        1  1776.290
#> 5  144-2024-GS-MZ1-R2        1        NA
#> 6  144-2024-GS-MZ1-R3        1   733.464

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
#> INFO  [2026-04-29 18:52:23] ** Raw data from DIANN imported successfully.
#> INFO  [2026-04-29 18:52:23] ** Filtering on Q.Value < 0.01
#> INFO  [2026-04-29 18:52:23] ** MBR was used to analyze the data. Now setting names and filtering
#> INFO  [2026-04-29 18:52:23] -- LibPGQValue < 0.01
#> INFO  [2026-04-29 18:52:23] -- LibQValue < 0.01
#> INFO  [2026-04-29 18:52:23] ** Raw data from DIANN cleaned successfully.
#> INFO  [2026-04-29 18:52:23] ** Using provided annotation.
#> INFO  [2026-04-29 18:52:23] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-29 18:52:23] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be removed.
#> INFO  [2026-04-29 18:52:23] ** Sequences containing DECOY, Decoys are removed.
#> INFO  [2026-04-29 18:52:23] ** Sequences containing \(UniMod\:35\) are removed.
#> INFO  [2026-04-29 18:52:23] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-29 18:52:23] ** Shared peptides are removed.
#> INFO  [2026-04-29 18:52:23] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: max
#> INFO  [2026-04-29 18:52:23] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-04-29 18:52:23] ** Run annotation merged with quantification data.
#> INFO  [2026-04-29 18:52:23] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-29 18:52:23] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.

head(msstatsptm_format$PTM)
#>       ProteinName                          PeptideSequence PrecursorCharge
#> 1 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 2 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 3 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 4 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 5 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#> 6 BRD2_HUMAN_K614 ASGSGGGSAALGPSGFGPSGGSGTK(UniMod:121)LPK               3
#>   FragmentIon ProductCharge IsotopeLabelType Condition BioReplicate  Run
#> 1      Frag10             1            Light   Control            1 Run1
#> 2      Frag10             1            Light   Control            2 Run2
#> 3      Frag10             1            Light   Control            3 Run3
#> 4      Frag10             1            Light   Control            4 Run4
#> 5      Frag10             1            Light Treatment            5 Run5
#> 6      Frag10             1            Light Treatment            6 Run6
#>   Fraction Intensity
#> 1        1  5521.414
#> 2        1  5652.361
#> 3        1  5783.308
#> 4        1  5914.254
#> 5        1  6045.201
#> 6        1  6176.148
```
