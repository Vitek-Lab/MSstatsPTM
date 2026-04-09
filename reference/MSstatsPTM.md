# MSstatsPTM: A package for detecting differentially abundant post-translational modifications (PTM) in mass spectrometry-based proteomic experiments.

A set of tools for detecting differentially abundant PTMs and proteins
in shotgun mass spectrometry-based proteomic experiments. The package
can handle a variety of acquisition types, including label free and TMT
experiments, acquired with DDA, DIA, SRM or PRM acquisition methods. The
package includes tools to convert raw data from different spectral
processing tools, summarize feature intensities, and fit a linear mixed
effects model. A major advantage of the package is to leverage a
separate global profiling run and adjust the PTM fold change for changes
in the unmodified protein, showing the unconvoluted PTM fold change.
Finally, the package includes functionality to plot a variety of data
visualizations.

## functions

- [`FragPipetoMSstatsPTMFormat`](https://vitek-lab.github.io/MSstatsPTM/reference/FragPipetoMSstatsPTMFormat.md)
  : Generates MSstatsPTM required input format for TMT FragePipe
  outputs.

- [`MaxQtoMSstatsPTMFormat`](https://vitek-lab.github.io/MSstatsPTM/reference/MaxQtoMSstatsPTMFormat.md)
  : Generates MSstatsPTM required input format for label-free and TMT
  MaxQuant outputs.

- [`ProgenesistoMSstatsPTMFormat`](https://vitek-lab.github.io/MSstatsPTM/reference/ProgenesistoMSstatsPTMFormat.md)
  : Generates MSstatsPTM required input format for label-free Progenesis
  outputs.

- [`SpectronauttoMSstatsPTMFormat`](https://vitek-lab.github.io/MSstatsPTM/reference/SpectronauttoMSstatsPTMFormat.md)
  : Generates MSstatsPTM required input format for label-free
  Spectronaut outputs.

- [`SkylinetoMSstatsPTMFormat`](https://vitek-lab.github.io/MSstatsPTM/reference/SkylinetoMSstatsPTMFormat.md)
  : Generates MSstatsPTM required input format for Skyline outputs.

- [`PStoMSstatsPTMFormat`](https://vitek-lab.github.io/MSstatsPTM/reference/PStoMSstatsPTMFormat.md)
  : Generates MSstatsPTM required input format for PEAKS outputs.

- [`PDtoMSstatsPTMFormat`](https://vitek-lab.github.io/MSstatsPTM/reference/PDtoMSstatsPTMFormat.md)
  : Generates MSstatsPTM required input format for Proteome Discoverer
  outputs.

- [`dataSummarizationPTM`](https://vitek-lab.github.io/MSstatsPTM/reference/dataSummarizationPTM.md)
  : Summarizes PSM level quantification to peptide (modification) and
  protein level quantification. For use in label-free analysis

- [`dataSummarizationPTM_TMT`](https://vitek-lab.github.io/MSstatsPTM/reference/dataSummarizationPTM_TMT.md)
  : Summarizes PSM level quantification to peptide (modification) and
  protein level quantification. For use in TMT analysis.

- [`dataProcessPlotsPTM`](https://vitek-lab.github.io/MSstatsPTM/reference/dataProcessPlotsPTM.md)
  : Visualization for explanatory data analysis. Specifically gives
  ability to plot Profile and Quality Control plots.

- [`groupComparisonPTM`](https://vitek-lab.github.io/MSstatsPTM/reference/groupComparisonPTM.md)
  : Tests for significant changes in PTM and protein abundance across
  conditions. Adjusts PTM fold change for changes in protein abundance.

- [`groupComparisonPlotsPTM`](https://vitek-lab.github.io/MSstatsPTM/reference/groupComparisonPlotsPTM.md)
  : Visualization for model-based analysis and summarization

## See also

Useful links:

- Report bugs at <https://github.com/Vitek-Lab/MSstatsPTM/issues>

## Author

**Maintainer**: Anthony Wu <wu.anthon@northeastern.edu>

Authors:

- Devon Kohler <kohler.d@northeastern.edu>

- Tsung-Heng Tsai <tsai.tsungheng@gmail.com>

- Deril Raju <raju.d@northeastern.edu>

- Ting Huang <thuang0703@gmail.com>

- Mateusz Staniak <mtst@mstaniak.pl>

- Meena Choi <mnchoi67@gmail.com>

- Olga Vitek <o.vitek@northeastern.edu>
