# Locate modified sites with a peptide

`locateMod` locates modified sites with a peptide.

## Usage

``` r
locateMod(peptide, aaStart, residueSymbol)
```

## Arguments

- peptide:

  A string. Peptide sequence.

- aaStart:

  An integer. Starting index of the peptide.

- residueSymbol:

  A string. Modification residue and denoted symbol.

## Value

A string.

## Examples

``` r
locateMod("P*EP*TIDE", 3, "\\*")
#> [1] 4 6
```
