# Read and tidy a FASTA file

`tidyFasta` reads and tidys FASTA file. Use this function as the first
step in identifying modification sites.

## Usage

``` r
tidyFasta(path)
```

## Arguments

- path:

  A string of path to a FASTA file.

## Value

A data.table with columns named `header`, `sequence`, `uniprot_ac`,
`uniprot_iso`, `entry_name`.

## Examples

``` r
tidyFasta(system.file("extdata", "O13297.fasta", package="MSstatsPTM"))
#>                                                                                                                                        header
#>                                                                                                                                        <char>
#> 1: sp|O13297|CET1_YEAST mRNA-capping enzyme subunit beta OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=CET1 PE=1 SV=2
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 sequence
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   <char>
#> 1: MSYTDNPPQTKRALSLDDLVNHDENEKVKLQKLSEAANGSRPFAENLESDINQTETGQAAPIDNYKESTGHGSHSQKPKSRKSSNDDEETDTDDEMGASGEINFDSEMDFDYDKQHRNLLSNGSPPMNDGSDANAKLEKPSDDSIHQNSKSDEEQRIPKQGNEGNIASNYITQVPLQKQKQTEKKIAGNAVGSVVKKEEEANAAVDNIFEEKATLQSKKNNIKRDLEVLNEISASSKPSKYRNVPIWAQKWKPTIKALQSINVKDLKIDPSFLNIIPDDDLTKSVQDWVYATIYSIAPELRSFIELEMKFGVIIDAKGPDRVNPPVSSQCVFTELDAHLTPNIDASLFKELSKYIRGISEVTENTGKFSIIESQTRDSVYRVGLSTQRPRFLRMSTDIKTGRVGQFIEKRHVAQLLLYSPKDSYDVKISLNLELPVPDNDPPEKYKSQSPISERTKDRVSYIHNDSCTRIDITKVENHNQNSKSRQSETTHEVELEINTPALLNAFDNITNDSKEYASLIRTFLNNGTIIRRKLSSLSYEIFEGSKKVM
#>    uniprot_ac uniprot_iso entry_name
#>        <char>      <char>     <char>
#> 1:     O13297      O13297 CET1_YEAST
```
