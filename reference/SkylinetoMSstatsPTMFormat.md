# Convert Skyline output into MSstatsPTM format

Currently only supports label-free quantification.

## Usage

``` r
SkylinetoMSstatsPTMFormat(
  input,
  fasta_path,
  fasta_protein_name = "uniprot_iso",
  annotation = NULL,
  input_protein = NULL,
  annotation_protein = NULL,
  use_unmod_peptides = FALSE,
  removeiRT = TRUE,
  filter_with_Qvalue = TRUE,
  qvalue_cutoff = 0.01,
  use_unique_peptide = TRUE,
  remove_few_measurements = FALSE,
  remove_oxidation_peptides = FALSE,
  removeProtein_with1Feature = FALSE,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = NULL
)
```

## Arguments

- input:

  name of Skyline PTM output

- fasta_path:

  A string of path to a FASTA file, used to match PTM peptides.

- fasta_protein_name:

  Name of fasta column that matches with protein name in evidence file.
  Default is `uniprot_iso`.

- annotation:

  name of 'annotation.txt' data which includes Condition, BioReplicate,
  Run. If annotation is already complete in Skyline, use annotation=NULL
  (default). It will use the annotation information from input.

- input_protein:

  name of Skyline unmodified protein output (optional)

- annotation_protein:

  name of 'annotation.txt' data which includes Condition, BioReplicate,
  Run for unmodified protein output. This can be the same as
  `annotation`.

- use_unmod_peptides:

  Boolean if the unmodified peptides in the input file should be used to
  construct the unmodified protein output. Only used if `input_protein`
  is not provided. Default is `FALSE`.

- removeiRT:

  TRUE (default) will remove the proteins or peptides which are labeld
  'iRT' in 'StandardType' column. FALSE will keep them.

- filter_with_Qvalue:

  TRUE(default) will filter out the intensities that have greater than
  qvalue_cutoff in DetectionQValue column. Those intensities will be
  replaced with zero and will be considered as censored missing values
  for imputation purpose.

- qvalue_cutoff:

  Cutoff for DetectionQValue. default is 0.01.

- use_unique_peptide:

  TRUE (default) removes peptides that are assigned for more than one
  proteins. We assume to use unique peptide for each protein.

- remove_few_measurements:

  TRUE will remove the features that have 1 or 2 measurements across
  runs. FALSE is default.

- remove_oxidation_peptides:

  TRUE will remove the peptides including 'oxidation (M)' in
  modification. FALSE is default.

- removeProtein_with1Feature:

  TRUE will remove the proteins which have only 1 feature, which is the
  combination of peptide, precursor charge, fragment and charge. FALSE
  is default.

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
