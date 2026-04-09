# Annotate modified sites with associated peptides

`PTMlocate` annotates modified sites with associated peptides.

## Usage

``` r
locatePTM(peptide, uniprot, fasta, modResidue, modSymbol, rmConfound = FALSE)
```

## Arguments

- peptide:

  A string vector of peptide sequences. The peptide sequence does not
  include its preceding and following AAs.

- uniprot:

  A string vector of Uniprot identifiers of the peptides' originating
  proteins. UniProtKB entry isoform sequence is used.

- fasta:

  A data.table with FASTA information. Output of `tidyFasta`.

- modResidue:

  A string. Modifiable amino acid residues.

- modSymbol:

  A string. Symbol of a modified site.

- rmConfound:

  A logical. `TRUE` removes confounded unmodified sites, `FALSE`
  otherwise. Default is `FALSE`.

## Value

A data frame with three columns: `uniprot_iso`, `peptide`, `site`.

## Examples

``` r
fasta = tidyFasta(system.file("extdata", "O13297.fasta", package="MSstatsPTM"))
locatePTM("DRVSYIHNDSC*TR", "O13297", fasta, "C", "\\*")
#> Key: <uniprot_iso>
#>    uniprot_iso        peptide   site
#>         <char>         <char> <char>
#> 1:      O13297 DRVSYIHNDSC*TR   C467
```
