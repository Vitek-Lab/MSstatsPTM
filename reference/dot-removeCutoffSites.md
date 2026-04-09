# Remove sites below cutoff probability

Remove sites below cutoff probability

## Usage

``` r
.removeCutoffSites(data, mod_pep_col, cutoff, remove_unlocalized_peptides)
```

## Arguments

- data:

  data.table

- mod_pep_col:

  column in data with modified sites

- cutoff:

  numeric cutoff. Default is .75.

- remove_unlocalized_peptides:

  Boolean if to remove peptides that could not be fully localized.

## Value

data.table with modifications below cutoff removed
