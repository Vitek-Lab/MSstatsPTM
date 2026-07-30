# Converts non-TMT Progenesis output into the format needed for MSstatsPTM

Converts non-TMT Progenesis output into the format needed for MSstatsPTM

## Usage

``` r
ProgenesistoMSstatsPTMFormat(
  ptm_input,
  annotation,
  global_protein_input = FALSE,
  fasta_path = FALSE,
  useUniquePeptide = TRUE,
  summaryforMultipleRows = max,
  removeFewMeasurements = TRUE,
  removeOxidationMpeptides = FALSE,
  removeProtein_with1Peptide = FALSE,
  mod.num = "Single"
)
```

## Arguments

- ptm_input:

  name of Progenesis output with modified peptides, which is
  wide-format. 'Accession', Sequence', 'Modification', 'Charge' and one
  column for each run are required

- annotation:

  name of 'annotation.txt' or 'annotation.csv' data which includes
  Condition, BioReplicate, Run, and Type (PTM or Protein) information.
  It will be matched with the column name of input for MS runs. Please
  note PTM and global Protein run names are often different, which is
  why an additional Type column indicating Protein or PTM is required.

- global_protein_input:

  name of Progenesis output with unmodified peptides, which is
  wide-format. 'Accession', Sequence', 'Modification', 'Charge' and one
  column for each run are required

- fasta_path:

  string containing path to the corresponding fasta file for the
  modified peptide dataset.

- useUniquePeptide:

  TRUE(default) removes peptides that are assigned for more than one
  proteins. We assume to use unique peptide for each protein.

- summaryforMultipleRows:

  max(default) or sum - when there are multiple measurements for certain
  feature and certain run, use highest or sum of multiple intensities.

- removeFewMeasurements:

  TRUE (default) will remove the features that have 1 or 2 measurements
  across runs.

- removeOxidationMpeptides:

  TRUE will remove the modified peptides including 'Oxidation (M)'
  sequence. FALSE is default.

- removeProtein_with1Peptide:

  TRUE will remove the proteins which have only 1 peptide and charge.
  FALSE is default.

- mod.num:

  For modified peptide dataset, must be one of `Single` or `Total`. The
  default is `Single`. The number modifications per peptide to be used.
  If "Single", only peptides with one modification will be used.
  Otherwise "Total" includes peptides with more than one modification.
  Selecting "Total" may confound the effect of different modifications.

## Value

a list of two data.tables named 'PTM' and 'PROTEIN' in the format
required by MSstatsPTM.

## Examples

``` r

input = system.file("tinytest/raw_data/Progenesis/progenesis_peptide.csv", 
                                package = "MSstatsPTM")
input = data.table::fread(input)
colnames(input) = unlist(input[1,])
input = input[-1,]
annot = system.file("tinytest/raw_data/Progenesis/phospho_annotation.csv", 
                                package = "MSstatsPTM")
annot = data.table::fread(annot)
prog_imported = ProgenesistoMSstatsPTMFormat(
    input, 
    annot
)
#> INFO  [2026-07-30 21:35:49] ** Raw data from Progenesis imported successfully.
#> INFO  [2026-07-30 21:35:49] ** Raw data from Progenesis cleaned successfully.
#> INFO  [2026-07-30 21:35:49] ** Using provided annotation.
#> INFO  [2026-07-30 21:35:49] ** Run labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-07-30 21:35:49] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements across runs will be removed.
#> INFO  [2026-07-30 21:35:49] ** Features with all missing measurements across runs are removed.
#> INFO  [2026-07-30 21:35:49] ** Shared peptides are removed.
#> INFO  [2026-07-30 21:35:49] ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows: max
#> INFO  [2026-07-30 21:35:49] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-07-30 21:35:49] ** Run annotation merged with quantification data.
#> INFO  [2026-07-30 21:35:49] ** Features with one or two measurements across runs are removed.
#> INFO  [2026-07-30 21:35:49] ** Fractionation handled.
#> INFO  [2026-07-30 21:35:49] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-07-30 21:35:49] ** Finished preprocessing. The dataset is ready to be processed by the dataProcess function.
head(prog_imported$PTM)
#>                                                               ProteinName
#> 1 sp|A2ASS6|TITIN_MOUSE_AVTSPPRVKSPEPR_[4] Phospho (ST)|[10] Phospho (ST)
#> 2 sp|A2ASS6|TITIN_MOUSE_AVTSPPRVKSPEPR_[4] Phospho (ST)|[10] Phospho (ST)
#> 3 sp|A2ASS6|TITIN_MOUSE_AVTSPPRVKSPEPR_[4] Phospho (ST)|[10] Phospho (ST)
#> 4 sp|A2ASS6|TITIN_MOUSE_AVTSPPRVKSPEPR_[4] Phospho (ST)|[10] Phospho (ST)
#> 5 sp|A2ASS6|TITIN_MOUSE_AVTSPPRVKSPEPR_[4] Phospho (ST)|[10] Phospho (ST)
#> 6 sp|A2ASS6|TITIN_MOUSE_AVTSPPRVKSPEPR_[4] Phospho (ST)|[10] Phospho (ST)
#>                            PeptideModifiedSequence PrecursorCharge FragmentIon
#> 1 AVTSPPRVKSPEPR[4] Phospho (ST)|[10] Phospho (ST)               2          NA
#> 2 AVTSPPRVKSPEPR[4] Phospho (ST)|[10] Phospho (ST)               2          NA
#> 3 AVTSPPRVKSPEPR[4] Phospho (ST)|[10] Phospho (ST)               2          NA
#> 4 AVTSPPRVKSPEPR[4] Phospho (ST)|[10] Phospho (ST)               2          NA
#> 5 AVTSPPRVKSPEPR[4] Phospho (ST)|[10] Phospho (ST)               2          NA
#> 6 AVTSPPRVKSPEPR[4] Phospho (ST)|[10] Phospho (ST)               2          NA
#>   ProductCharge IsotopeLabelType Condition BioReplicate           Run Fraction
#> 1            NA                L         A           S1 001_PH003_G23        1
#> 2            NA                L         B           S2 002_PH003_G23        1
#> 3            NA                L         A           S3 003_PH003_G23        1
#> 4            NA                L         A           S4 004_PH003_G23        1
#> 5            NA                L         B           S5 005_PH003_G23        1
#> 6            NA                L         B           S6 006_PH003_G23        1
#>   Intensity
#> 1  31694748
#> 2  71177637
#> 3  63130816
#> 4  56450448
#> 5  43183206
#> 6  41647175
```
