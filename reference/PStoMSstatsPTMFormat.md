# Convert Peaks Studio output into MSstatsPTM format

Currently only supports label-free quantification.

## Usage

``` r
PStoMSstatsPTMFormat(
  input,
  annotation,
  input_protein = NULL,
  annotation_protein = NULL,
  use_unmod_peptides = FALSE,
  target_modification = NULL,
  remove_oxidation_peptides = FALSE,
  remove_multi_mod_types = FALSE,
  summaryforMultipleRows = max,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = NULL
)
```

## Arguments

- input:

  name of Peaks Studio PTM output

- annotation:

  name of annotation file which includes Raw.file, Condition,
  BioReplicate, Run. For example annotation see example below.

- input_protein:

  name of Peaks Studio unmodified protein output (optional)

- annotation_protein:

  name of annotation file which includes Raw.file, Condition,
  BioReplicate, Run for unmodified protein output.

- use_unmod_peptides:

  Boolean if the unmodified peptides in the input file should be used to
  construct the unmodified protein output. Only used if `input_protein`
  is not provided. Default is `FALSE`

- target_modification:

  Character name of modification of interest. To use all mod types,
  leave as `NULL`. Default is `NULL`. Note that if the name includes
  special characters, you must include "\\ before the characters. Ex.
  "Phosphorylation \\STY\\"

- remove_oxidation_peptides:

  Boolean if Oxidation (M) modifications should be removed. Default is
  `FALSE`

- remove_multi_mod_types:

  Used if `target_modification` is not `NULL`. `TRUE` will remove
  peptides with multiple types of modifications (ie acetylation and
  phosphorylation). `FALSE` will keep these peptides and summarize them
  seperately.

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

`list` of `data.table`

## Examples

``` r
# The output should be in the following format.
head(raw.input$PTM)
#> # A tibble: 6 × 10
#>   ProteinName PeptideSequence Condition BioReplicate Run        Intensity
#>   <chr>       <chr>           <chr>     <chr>        <chr>          <dbl>
#> 1 Q9UHD8_K262 DAGLK*QAPASR    CCCP      BCH1         CCCP-B1T1   1423906.
#> 2 Q9UHD8_K262 DAGLK*QAPASR    CCCP      BCH1         CCCP-B1T2    877045.
#> 3 Q9UHD8_K262 DAGLK*QAPASR    CCCP      BCH2         CCCP-B2T1    384418.
#> 4 Q9UHD8_K262 DAGLK*QAPASR    CCCP      BCH2         CCCP-B2T2    454858.
#> 5 Q9UHD8_K262 DAGLK*QAPASR    Combo     BCH1         Combo-B1T1  1603377.
#> 6 Q9UHD8_K262 DAGLK*QAPASR    Combo     BCH1         Combo-B1T2   676555.
#> # ℹ 4 more variables: PrecursorCharge <chr>, FragmentIon <lgl>,
#> #   ProductCharge <lgl>, IsotopeLabelType <chr>
head(raw.input$PROTEIN)
#> # A tibble: 6 × 10
#>   ProteinName PeptideSequence Condition BioReplicate Run           Intensity
#>   <chr>       <chr>           <chr>     <chr>        <chr>             <dbl>
#> 1 Q9UHD8      STLINTLFK       CCCP      BCH2         CCCP-B2T1       367944.
#> 2 Q9UHD8      STLINTLFK       CCCP      BCH2         CCCP-B2T2       341207.
#> 3 Q9UHD8      STLINTLFK       Combo     BCH2         Combo-B2T1      185843.
#> 4 Q9UHD8      STLINTLFK       Ctrl      BCH2         Ctrl-B2T1       529224.
#> 5 Q9UHD8      STLINTLFK       Ctrl      BCH2         Ctrl-B2T2       483355.
#> 6 Q9UHD8      STLINTLFK       USP30_OE  BCH2         USP30_OE-B2T1   447795.
#> # ℹ 4 more variables: PrecursorCharge <chr>, FragmentIon <lgl>,
#> #   ProductCharge <lgl>, IsotopeLabelType <chr>
```
