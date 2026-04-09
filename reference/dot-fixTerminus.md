# Fix terminus location adjustments

Fix terminus location adjustments

## Usage

``` r
.fixTerminus(data, terminus_id, unmod_pep_col)
```

## Arguments

- data:

  data.table containing peptide data

- terminus_id:

  character string identifying the terminus (e.g. N-terminus)

- unmod_pep_col:

  character string specifying the column name containing unmodified
  peptide sequences

## Value

data.table with corrected Start positions
