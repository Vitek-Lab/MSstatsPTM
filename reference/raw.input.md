# Example of input PTM dataset for LabelFree/DDA/DIA experiments.

It can be the output of MSstatsPTM converter
ProgenesistoMSstatsPTMFormat or other MSstats converter functions
(Please see MSstatsPTM_LabelFree_Workflow vignette). The dataset is
formatted as a list with two data.tables named PTM and PROTEIN. In each
data.table the variables are as follows:

## Usage

``` r
raw.input
```

## Format

A list of two data.tables named PTM and PROTEIN with 1745 and 478 rows
respectively.

## Details

- ProteinName : Name of protein with modification site mapped in with an
  underscore. ie "Protein_4_Y474"

- PeptideSequence

- Condition : Condition (ex. Healthy, Cancer, Time0)

- BioReplicate : Unique ID for biological subject.

- Run : MS run ID.

- Intensity

- PrecursorCharge

- FragmentIon

- ProductCharge

- IsotopeLabelType

## Examples

``` r
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
