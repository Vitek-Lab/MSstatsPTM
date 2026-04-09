# MSstatsPTM

MSstatsPTM provides general statistical methods for quantitative
characterization of post-translational modifications (PTMs). Supports
DDA, DIA, and tandem mass tag (TMT) labeling. Typically, the analysis
involves the quantification of PTM sites (i.e., modified residues) and
their corresponding proteins, as well as the integration of the
quantification results. MSstatsPTM provides functions for summarization,
estimation of PTM site abundance, and detection of changes in PTMs
across experimental conditions.

## Installation

The most recent version of the package is available on Bioconductor:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("MSstatsPTM")

Optionally, the package can be installed directly from Github:

    devtools::install_github("Vitek-Lab/MSstatsPTM", build_vignettes = TRUE)

## Contributing

We welcome contributions from the community. For details on how to
contribute to the development of MSstatsPTM, please refer to the
[CONTRIBUTING](https://github.com/Vitek-Lab/MSstatsPTM/blob/master/.github/CONTRIBUTING.md)
file.

## License

[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0)
