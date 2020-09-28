#' MSstatsPTM: A package for statistical characterization of PTMs
#'
#' The MSstatsPTM package provides four main functions for quantitative analysis
#' of PTMs
#'
#' Quantitative analyses of PTMs are supported by four main functions of
#' _MSstatsPTM_:
#'
#' @section Normalization:
#' \code{PTMnormalize()} normalizes the quantified peak intensities to correct
#' systematic variation across MS runs.
#'
#' @section Summarization:
#' \code{PTMsummarize()} summarizes log2-intensities of spectral features
#' (i.e., precursor ions in DDA, fragments in DIA, or transitions in SRM) into
#' one value per PTM site per run or one value per protein per run.
#'
#' @section Estimation:
#' \code{PTMestimate()} takes as input the summarized log2-intensities for
#' each PTM site, performs statistical modeling for the log2-abundance of the
#' site, and returns the estimates of model parameters for all PTM sites in all
#' experimental conditions.
#'
#' @section Comparison:
#' \code{PTMcompareMeans()} performs statistical testing for detecting changes
#' in PTM mean abundances between conditions.
#'
#' @name MSstatsPTM
#' @docType package
#' @importFrom stats df.residual lm median medpolish p.adjust pt
#' @importFrom rlang .data
#' @importFrom tidyr nest unnest
#' @importFrom tidyselect one_of
#' @importFrom dplyr left_join inner_join semi_join bind_rows group_by ungroup
#' summarise mutate count
#' @importFrom tibble tibble as_tibble
NULL
#> NULL
