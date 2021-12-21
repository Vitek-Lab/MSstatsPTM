# MSstatsPTM

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/Vitek-Lab/MSstatsPTM.svg?branch=master)](https://travis-ci.org/Vitek-Lab/MSstatsPTM)
[![Codecov test coverage](https://codecov.io/gh/Vitek-Lab/MSstatsPTM/branch/master/graph/badge.svg)](https://codecov.io/gh/Vitek-Lab/MSstatsPTM?branch=master)
<!-- badges: end -->


The goal of MSstatsPTM is to provide the implementation of general statistical methods for characterization of quantitative changes in global post-translational modification profiling experiments.

## Installation 

The most recent version of the package is available on Bioconductor:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsPTM")
```

Optionally, the package can be installed directly from Github:

```
devtools::install_github("Vitek-Lab/MSstatsPTM", build_vignettes = TRUE)
```

## Contributing

We welcome contributions from the community. For details on how to contribute to the
development of MSstatsPTM, please refer to the [CONTRIBUTING](https://github.com/Vitek-Lab/MSstatsPTM/blob/master/.github/CONTRIBUTING.md) file.

## License

[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0)
