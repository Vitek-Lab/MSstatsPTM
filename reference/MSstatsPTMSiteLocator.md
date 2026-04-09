# Locate modification site number and amino acid

Locate modification site number and amino acid

## Usage

``` r
MSstatsPTMSiteLocator(
  data,
  protein_name_col = "ProteinName",
  unmod_pep_col = "PeptideSequence",
  mod_pep_col = "PeptideModifiedSequence",
  clean_mod = FALSE,
  fasta_file = NULL,
  fasta_protein_name = "header",
  mod_id = "\\*",
  localization_scores = FALSE,
  localization_cutoff = 0.75,
  remove_unlocalized_peptides = TRUE,
  terminus_included = FALSE,
  terminus_id = "\\.",
  mod_id_is_numeric = FALSE,
  remove_underscores = FALSE,
  remove_other_mods = FALSE,
  bracket = FALSE,
  replace_text = FALSE
)
```

## Arguments

- data:

  `data.table` of enriched experimental run. Must include `ProteinName`,
  `PeptideSequence`, `PeptideModifiedSequence`, and (optionally) `Start`
  columns.

- protein_name_col:

  Name of column indicating protein. Default is `ProteinName`.

- unmod_pep_col:

  Name of column indicating unmodified peptide sequence. Default is
  `PeptideSequence`.

- mod_pep_col:

  Name of column indicating modified peptide sequence. Default is
  `PeptideModifiedSequence`.

- clean_mod:

  Remove special characters and numbers around modification name.
  Default is `FALSE`

- fasta_file:

  File path to FASTA file that matches with proteins in `data`. Can be
  either string or `data.table` processed with
  [`tidyFasta()`](https://vitek-lab.github.io/MSstatsPTM/reference/tidyFasta.md)
  function. Default to NULL if peptide number included in `data`.

- fasta_protein_name:

  Name of fasta file column that matches with `protein_name_col`.
  Default is `header`.

- mod_id:

  String that indicates what amino acid was modified in
  `PeptideSequence`.

- localization_scores:

  Boolean indicating if mod id is a localization score. If TRUE,
  `mod_id` will be ignored and localization cutoff will be used to
  determine sites. Default is FALSE.

- localization_cutoff:

  Default is .75. Localization probabilities below cutoffs will be
  removed. `localization_scores` must be TRUE.

- remove_unlocalized_peptides:

  Default is TRUE. If `localization_scores` is TRUE and probabilities
  are below `localization_cutoff`, the modification site will not be
  able to be determined. These unlocalized peptides can be kept or
  removed. If FALSE the unlocalized peptides will still be used in
  modeling the sites that could be localized.

- terminus_included:

  Boolean indicating if the `PeptideSequence` includes the terminus
  amino acid.

- terminus_id:

  String that indicates what the terminus amino acid is. Default is '.'.

- mod_id_is_numeric:

  Boolean indicating if modification identifier is a number instead of a
  character (i.e. +80 vs \*).

- remove_underscores:

  Boolean indicating if underscores around peptide exist. These should
  be removed to properly count where in sequence the modification
  occurred.

- remove_other_mods:

  keeping mods that are not of interest can mess up the amino acid
  count. Remove them if they are causing issues.

- bracket:

  bracket type that encompasses PTM (usually `[` or `(`). Always pass
  opening bracket (there is a function to grab the close bracket).
  Default is FALSE (i.e. no bracket).

- replace_text:

  If PTM is noted by text (i.e. `Phospho`) and needs to be replaced by
  an indicator (`*`)

## Value

`data.table` with site location added into `Protein` column.

## Examples

``` r
##TODO
```
