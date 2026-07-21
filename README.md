# MSstatsPTM

<!-- badges: start -->
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/MSstatsPTM.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MSstatsPTM/)
[![Codecov test coverage](https://codecov.io/gh/Vitek-Lab/MSstatsPTM/branch/devel/graph/badge.svg)](https://codecov.io/gh/Vitek-Lab/MSstatsPTM/branch/devel)
[![Bioconductor Downloads Rank](https://bioconductor.org/shields/downloads/release/MSstatsPTM.svg)](https://bioconductor.org/packages/stats/bioc/MSstatsPTM/)
[![Years in Bioconductor](https://bioconductor.org/shields/years-in-bioc/MSstatsPTM.svg)](https://bioconductor.org/packages/release/bioc/html/MSstatsPTM.html#since)
[![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
<!-- badges: end -->

MSstatsPTM is an R/Bioconductor package for statistical relative quantification
of post-translational modifications (PTMs) in mass spectrometry-based proteomics
experiments. It supports label-free DDA, DIA, and SRM workflows as well as
tandem mass tag (TMT) labeling. The analysis jointly models the abundance of
modified peptides (PTM sites) and their corresponding unmodified proteins, then
adjusts the PTM-level changes for changes in overall protein abundance so that
detected differences reflect genuine changes in modification stoichiometry
rather than changes in protein expression. MSstatsPTM provides functions for
summarization, estimation of PTM-site abundance, and detection of changes in
PTMs across experimental conditions.

MSstatsPTM is part of the [MSstats](https://github.com/Vitek-Lab/MSstats)
family of packages, developed and maintained by the
[Vitek Lab](https://olga-vitek-lab.khoury.northeastern.edu/) at Northeastern
University. The package and its documentation are also available at
[msstats.org](http://msstats.org).

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsPTM")
```

The development version can be installed directly from this repository:

```r
remotes::install_github("Vitek-Lab/MSstatsPTM")
```

## Quick Start

```r
library(MSstatsPTM)

# Example label-free PTM + unmodified-protein data
data(raw.input)

# Summarize PTM and protein data to run level
summary <- dataSummarizationPTM(raw.input, use_log_file = FALSE)

# Define the comparison of interest
comparison <- matrix(c(-1, 0, 1, 0), nrow = 1)
row.names(comparison) <- "CCCP-Ctrl"
colnames(comparison) <- c("CCCP", "Combo", "Ctrl", "USP30_OE")

# Model PTM and protein levels, adjusting PTMs for protein abundance
model <- groupComparisonPTM(summary,
                            ptm_label_type = "LF",
                            protein_label_type = "LF",
                            contrast.matrix = comparison,
                            use_log_file = FALSE)

head(model$ADJUSTED.Model)
```

Real experiments typically start by converting a search tool's output into
MSstatsPTM format with one of the `*toMSstatsPTMFormat()` converters below, then
proceed with `dataSummarizationPTM()` and `groupComparisonPTM()` as above. See
the vignettes for complete, tool-specific examples.

## Supported Converters

MSstatsPTM does not read raw search-tool output directly. Instead, a converter
translates each tool's PTM and protein reports into MSstatsPTM format before
`dataSummarizationPTM()` is called:

| Search tool / format | Converter function |
| --- | --- |
| Skyline | `SkylinetoMSstatsPTMFormat()` |
| MaxQuant | `MaxQtoMSstatsPTMFormat()` |
| Progenesis | `ProgenesistoMSstatsPTMFormat()` |
| Spectronaut | `SpectronauttoMSstatsPTMFormat()` |
| Proteome Discoverer | `PDtoMSstatsPTMFormat()` |
| DIA-NN | `DIANNtoMSstatsPTMFormat()` |
| FragPipe | `FragPipetoMSstatsPTMFormat()` |
| Metamorpheus | `MetamorpheusToMSstatsPTMFormat()` |
| Protein Prospector | `ProteinProspectortoMSstatsPTMFormat()` |
| Peak Studio | `PStoMSstatsPTMFormat()` |

See the [label-free](vignettes/MSstatsPTM_LabelFree_Workflow.Rmd) and
[TMT](vignettes/MSstatsPTM_TMT_Workflow.Rmd) workflow vignettes for the required
input files and options for each converter.

## Documentation

- [MSstatsPTM label-free workflow](vignettes/MSstatsPTM_LabelFree_Workflow.Rmd) — worked example for label-free DDA/DIA/SRM
- [MSstatsPTM TMT workflow](vignettes/MSstatsPTM_TMT_Workflow.Rmd) — worked example for TMT labeling
- [Official website: msstats.org](http://msstats.org)
- [Bioconductor package page and reference manual](https://bioconductor.org/packages/MSstatsPTM)

## Getting Help / Reporting Bugs

- **Questions about usage, statistical methods, or troubleshooting:** please
  post to the [MSstats Google Group](https://groups.google.com/forum/#!forum/msstats).
  This is monitored by the development team and searchable, so it's the fastest
  way to get help and to see if your question has already been answered.
- **Bug reports and feature requests for this repository:** please open a
  [GitHub issue](https://github.com/Vitek-Lab/MSstatsPTM/issues).

## References

If you use MSstatsPTM, please cite:

1. Kohler D, Tsai TH, Verschueren E, Huang T, Hinkle T, Phu L, Choi M, Vitek O.
   **MSstatsPTM: Statistical Relative Quantification of Posttranslational
   Modifications in Bottom-Up Mass Spectrometry-Based Proteomics.**
   *Mol Cell Proteomics*. 2022;22(1):100477.
   [DOI: 10.1016/j.mcpro.2022.100477](https://doi.org/10.1016/j.mcpro.2022.100477)

## Funding

MSstats development has been supported by the Chan Zuckerberg Initiative's
[Essential Open Source Software for Science](https://chanzuckerberg.com/eoss/proposals/).

## License

MSstatsPTM is released under the [Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
