# Annotate modification site

`annotSite` annotates modified sites as their residues and locations.

## Usage

``` r
annotSite(aaIndex, residue, lenIndex = NULL)
```

## Arguments

- aaIndex:

  An integer vector. Location of the sites.

- residue:

  A string vector. Amino acid residue.

- lenIndex:

  An integer. Default is `NULL`

## Value

A string.

## Examples

``` r
annotSite(10, "K")
#> [1] "K10"
annotSite(10, "K", 3L)
#> [1] "K010"
```
