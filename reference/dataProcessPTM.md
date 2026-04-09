# Data processing and summarization of peptide-level quantification to PTM and protein level quantification

Function to perform data processing and summarization on an experiment
targeting post-translational modifications. Performs normalization,
missing value imputation, feature selection, and summarization. Can
optionally take an additional global protein quantification experiment
for protein-level correction of PTM changes. Can take either label free
or tandem mass tag (TMT) labeled data.

## Usage

``` r
dataProcessPTM(
  data,
  ptm_label_type,
  protein_label_type,
  MBimpute_ptm = FALSE,
  MBimpute_protein = TRUE,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = NULL,
  ...
)
```

## Arguments

- data:

  Name of the output of MSstatsPTM converter function or peptide-level
  quantified data from other tools. It should be a list containing one
  or two data tables, named PTM and PROTEIN for modified and unmodified
  datasets. The list must at least contain the PTM dataset, however the
  PROTEIN dataset is optional.

- ptm_label_type:

  Indicator of labeling type for PTM dataset. Must be one of `LF` or
  `TMT`

- protein_label_type:

  Indicator of labeling type for PROTEIN dataset. Must be one of `LF` or
  `TMT`

- MBimpute_ptm:

  Missing value imputation on the feature-level for the PTM dataset.
  TRUE (default) imputes missing values by Accelated failure model.
  FALSE uses minimum value to impute the missing value for each peptide
  precursor ion.

- MBimpute_protein:

  Missing value imputation on the feature-level for the PROTEIN dataset.
  TRUE (default) imputes missing values by Accelated failure model.
  FALSE uses minimum value to impute the missing value for each peptide
  precursor ion.

- ...:

  Additional parameters passed to either the
  [dataProcess](https://rdrr.io/pkg/MSstats/man/dataProcess.html)
  function from `MSstats` or the
  [proteinSummarization](https://rdrr.io/pkg/MSstatsTMT/man/proteinSummarization.html)
  function from `MSstatsTMT`.
