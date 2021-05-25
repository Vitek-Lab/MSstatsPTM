#' MSstatsPTM: A package for detecting differencially abundant post
#' translational modifications (PTM) in shotgun mass spectrometry-bsed proteomic
#'  experiments.
#'
#' A set of tools for detecting differentially abundant PTMs and proteins in
#' shotgun mass spectrometry-based proteomic experiments. The package can handle
#' a variety of acquisition types, including label free, DDA, DIA, and TMT. The 
#' package includes tools to convert raw data from different spectral processing
#' tools, summarize feature intensities, and fit a linear mixed effects model.
#' Additionally the package includes functionality to plot a variety of data
#' visualizations.
#'
#' @section functions :
#' \itemize{
#'   \item \code{\link{MaxQtoMSstatsPTMFormat}} : Generates MSstatsPTM required
#'   input format for TMT MaxQuant outputs.
#'   \item \code{\link{ProgenesistoMSstatsPTMFormat}} : Generates MSstatsPTM 
#'   required input format for non-TMT Proteoviz outputs.
#'   \item \code{\link{SpectronauttoMSstatsPTMFormat}} : Generates MSstatsPTM 
#'   required input format for non-TMT Spectronaut outputs.
#'   \item \code{\link{dataSummarizationPTM}} : Summarizes PSM level 
#'   quantification to peptide (modification) and protein level quantification.
#'   For use in non-TMT analysis
#'   \item \code{\link{dataSummarizationPTM_TMT}} : Summarizes PSM level 
#'   quantification to peptide (modification) and protein level quantification.
#'   For use in TMT analysis.
#'   \item \code{\link{dataProcessPlotsPTM}} : Visualization for explanatory
#'   data analysis. Specifically gives ability to plot Profile and Quality 
#'   Control plots.
#'   \item \code{\link{groupComparisonPTM}} : Tests for significant changes in 
#'   PTM and protein abundance across conditions. Adjusts PTM fold change for 
#'   changes in protein abundance.
#'   \item \code{\link{groupComparisonPlotsPTM}} : Visualization for model-based
#'    analysis and summarization
#' }
#'
#' @docType package
#' @name MSstatsPTM
NULL
