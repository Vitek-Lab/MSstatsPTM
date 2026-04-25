# Import Metamorpheus files into PTM format

Import Metamorpheus files into PTM format

## Usage

``` r
MetamorpheusToMSstatsPTMFormat(
  input,
  annotation,
  fasta_path,
  input_protein = NULL,
  annotation_protein = NULL,
  use_unmod_peptides = FALSE,
  mod_ids = c("\\[Common Biological:Phosphorylation on S\\]"),
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

  name of Metamorpheus output file, which is tabular format. Use the
  AllQuantifiedPeaks.tsv file from the Metamorpheus output.

- annotation:

  name of 'annotation.txt' data which includes Condition, BioReplicate.

- fasta_path:

  string containing path to the corresponding fasta file for the
  modified peptide dataset.

- input_protein:

  same as `input` for global profiling run. Default is NULL.

- annotation_protein:

  same as `annotation` for global profiling run. Default is NULL.

- use_unmod_peptides:

  If `protein_input` is not provided, unmodified peptides can be
  extracted from `input` to be used in place of a global profiling run.
  Default is `FALSE`.

- mod_ids:

  List of modifications of interest. Default is a list with only
  `Common Biological:Phosphorylation on S`. Please note that the
  'mod_ids' parameter currently supports lists of size 1 only. Future
  updates aim to extend its functionality to accommodate lists of
  greater sizes. Note `\\` must be included before special characters.

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

## Author

Anthony Wu

## Examples

``` r
input = system.file("tinytest/raw_data/Metamorpheus/AllQuantifiedPeaks.tsv", 
                                package = "MSstatsPTM")
input = data.table::fread(input)
annot = system.file("tinytest/raw_data/Metamorpheus/ExperimentalDesign.tsv", 
                                package = "MSstatsPTM")
annot = data.table::fread(annot)
input_protein = system.file("tinytest/raw_data/Metamorpheus/AllQuantifiedPeaksGlobalProteome.tsv",
                                package = "MSstatsPTM")
input_protein = data.table::fread(input_protein)
annot_protein = system.file("tinytest/raw_data/Metamorpheus/ExperimentalDesignGlobalProteome.tsv", 
                                package = "MSstatsPTM")
annot_protein = data.table::fread(annot_protein)
fasta_path=system.file("extdata", "metamorpheus_fasta.fasta", 
                                package="MSstatsPTM")
metamorpheus_imported = MetamorpheusToMSstatsPTMFormat(
    input, 
    annot, 
    fasta_path=fasta_path,
    input_protein=input_protein,
    annotation_protein=annot_protein,
    use_unmod_peptides=FALSE,
    mod_ids = c("\\[Common Fixed:Carbamidomethyl on C\\]")
)
#> [1] "FASTA file missing 3 Proteins. These will be removed. This may be due to non-unique identifications."
#> INFO  [2026-04-25 15:36:36] ** Raw data from Metamorpheus imported successfully.
#> INFO  [2026-04-25 15:36:36] ** Raw data from Metamorpheus cleaned successfully.
#> INFO  [2026-04-25 15:36:36] ** Using provided annotation.
#> INFO  [2026-04-25 15:36:36] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-25 15:36:36] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be removed.
#> INFO  [2026-04-25 15:36:36] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-25 15:36:36] ** Shared peptides are removed.
#> INFO  [2026-04-25 15:36:36] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: max
#> INFO  [2026-04-25 15:36:36] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-04-25 15:36:36] ** Run annotation merged with quantification data.
#> INFO  [2026-04-25 15:36:36] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-04-25 15:36:36] ** Fractionation handled.
#> INFO  [2026-04-25 15:36:36] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-25 15:36:36] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.
#> INFO  [2026-04-25 15:36:36] ** Raw data from Metamorpheus imported successfully.
#> INFO  [2026-04-25 15:36:36] ** Raw data from Metamorpheus cleaned successfully.
#> INFO  [2026-04-25 15:36:36] ** Using provided annotation.
#> INFO  [2026-04-25 15:36:36] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-25 15:36:36] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be removed.
#> INFO  [2026-04-25 15:36:36] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-04-25 15:36:36] ** Shared peptides are removed.
#> INFO  [2026-04-25 15:36:36] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: max
#> INFO  [2026-04-25 15:36:36] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-04-25 15:36:36] ** Run annotation merged with quantification data.
#> INFO  [2026-04-25 15:36:36] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-04-25 15:36:36] ** Fractionation handled.
#> INFO  [2026-04-25 15:36:36] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-25 15:36:36] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.
head(metamorpheus_imported$PTM)
#>   ProteinName
#> 1 P06748_C104
#> 2 P06748_C104
#> 3 P06748_C104
#> 4 P06748_C104
#> 5 P06748_C104
#> 6 P06748_C104
#>                                                                                  PeptideSequence
#> 1 C[Common Fixed:Carbamidomethyl on C]GSGPVHISGQHLVAVEEDAES[UniProt:Phosphoserine on S]EDEEEEDVK
#> 2 C[Common Fixed:Carbamidomethyl on C]GSGPVHISGQHLVAVEEDAES[UniProt:Phosphoserine on S]EDEEEEDVK
#> 3 C[Common Fixed:Carbamidomethyl on C]GSGPVHISGQHLVAVEEDAES[UniProt:Phosphoserine on S]EDEEEEDVK
#> 4 C[Common Fixed:Carbamidomethyl on C]GSGPVHISGQHLVAVEEDAES[UniProt:Phosphoserine on S]EDEEEEDVK
#> 5 C[Common Fixed:Carbamidomethyl on C]GSGPVHISGQHLVAVEEDAES[UniProt:Phosphoserine on S]EDEEEEDVK
#> 6 C[Common Fixed:Carbamidomethyl on C]GSGPVHISGQHLVAVEEDAES[UniProt:Phosphoserine on S]EDEEEEDVK
#>   PrecursorCharge FragmentIon ProductCharge IsotopeLabelType Condition
#> 1               3          NA            NA                L   control
#> 2               3          NA            NA                L   control
#> 3               3          NA            NA                L     1_min
#> 4               3          NA            NA                L     1_min
#> 5               3          NA            NA                L    10_min
#> 6               3          NA            NA                L    10_min
#>   BioReplicate                                                   Run Fraction
#> 1            1    180323acs_001_HUVEC_P5_no_VEGF_control_1_001-calib        1
#> 2            2    180323acs_002_HUVEC_P5_no_VEGF_control_1_002-calib        1
#> 3            3  180323acs_003_HUVEC_P5_50ng_mL_VEGF_1min_1_001-calib        1
#> 4            4  180323acs_004_HUVEC_P5_50ng_mL_VEGF_1min_1_002-calib        1
#> 5            5 180323acs_007_HUVEC_P5_50ng_mL_VEGF_10min_1_001-calib        1
#> 6            6 180323acs_008_HUVEC_P5_50ng_mL_VEGF_10min_1_002-calib        1
#>   Intensity
#> 1  224148.4
#> 2        NA
#> 3        NA
#> 4  251253.7
#> 5  188734.7
#> 6        NA
head(metamorpheus_imported$PROTEIN)
#>        ProteinName                                PeptideSequence
#> 1 Q8N3V7;UNDEFINED AAS[UniProt:Phosphoserine on S]PAKPSSLDLVPNLPK
#> 2 Q8N3V7;UNDEFINED AAS[UniProt:Phosphoserine on S]PAKPSSLDLVPNLPK
#> 3 Q8N3V7;UNDEFINED AAS[UniProt:Phosphoserine on S]PAKPSSLDLVPNLPK
#> 4 Q8N3V7;UNDEFINED AAS[UniProt:Phosphoserine on S]PAKPSSLDLVPNLPK
#> 5 Q8N3V7;UNDEFINED AAS[UniProt:Phosphoserine on S]PAKPSSLDLVPNLPK
#> 6 Q8N3V7;UNDEFINED AAS[UniProt:Phosphoserine on S]PAKPSSLDLVPNLPK
#>   PrecursorCharge FragmentIon ProductCharge IsotopeLabelType Condition
#> 1               3          NA            NA                L   control
#> 2               3          NA            NA                L   control
#> 3               3          NA            NA                L     1_min
#> 4               3          NA            NA                L     1_min
#> 5               3          NA            NA                L    10_min
#> 6               3          NA            NA                L    10_min
#>   BioReplicate                                                   Run Fraction
#> 1            7    180323acs_001_HUVEC_P5_no_VEGF_control_1_001-calib        1
#> 2            8    180323acs_002_HUVEC_P5_no_VEGF_control_1_002-calib        1
#> 3            9  180323acs_003_HUVEC_P5_50ng_mL_VEGF_1min_1_001-calib        1
#> 4           10  180323acs_004_HUVEC_P5_50ng_mL_VEGF_1min_1_002-calib        1
#> 5           11 180323acs_007_HUVEC_P5_50ng_mL_VEGF_10min_1_001-calib        1
#> 6           12 180323acs_008_HUVEC_P5_50ng_mL_VEGF_10min_1_002-calib        1
#>   Intensity
#> 1  602295.6
#> 2  572816.1
#> 3  815920.8
#> 4  856970.5
#> 5  776769.5
#> 6  853108.7
```
