#' MSstatsPTM: A package for statistical characterization of PTMs.
#'
#' The MSstatsPTM package provides three categories of important functions for
#' quantitative analysis of PTMs: summarization, estimation and comparison
#'
#' @section Summarization functions:
#' The \code{PTMsummarize} and related functions summarize log2-intensities of
#' spectral features for each PTM site into one value per run.
#'
#' @section Estimation functions:
#' The \code{PTMestimate} and related functions take as input the summarized
#' log2-intensities for each PTM site, performs statistical modeling for the
#' abundance of the site, and returns the estimates of model parameters for all
#' sites in all experimental conditions.
#'
#' @section Comparison functions:
#' The \code{PTMcompareMeans} and related functions perform significance
#' analysis for detecting changes in PTM mean abundances between conditions.
#'
#' @name MSstatsPTM
#' @docType package
#' @importFrom magrittr %>%
#' @importFrom tidyr nest unnest
#' @importFrom tidyselect one_of
#' @importFrom dplyr left_join inner_join bind_rows group_by ungroup summarise mutate
#' @importFrom tibble tibble as_tibble
NULL
#> NULL
