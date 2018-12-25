#' Summarize feature log-intensities for each PTM site per run.
#'
#' \code{summarize_ptm} summarizes feature log-intensities for each PTM site.
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr nest unnest
#' @importFrom dplyr filter select mutate group_by left_join distinct
#' @importFrom purrr map
#' @param df_site A data frame that includes columns of \code{protein},
#'   \code{site}, \code{group}, \code{run}, \code{feature}, and \code{log2inty}.
#' @param method A string defining the summarization method. Default is
#'   \code{"tmp"}, which applies Tukey's median polish.
#' @param impute_missing A logical. \code{TRUE} performs missing value
#'   imputation confounded with accelerated failure time model, \code{FALSE}
#'   does not handle missing values. Default is \code{FALSE}.
#' @return A data frame of five columns: \code{protein}, \code{site},
#'   \code{group}, \code{run}, and \code{log2inty}.
#' @export
#'
#' @examples
#' summarize_ptm(df_site)
summarize_ptm <- function(df_site, method = "tmp", impute_missing = FALSE) {

    if (!all(c("protein", "site", "group", "run", "feature", "log2inty") %in% names(df_site))) {
        stop("Missing columns! Please make sure the input data frame contains
             the columns of 'protein', 'site', 'group', 'run', 'feature', 'log2inty'")
    }

    df_site <- df_site %>%
        filter(!is.na(log2inty))

    # Experimental design
    design <- df_site %>%
        distinct(run, group)

    # Nested data frame
    if ("batch" %in% names(df_site)) {
        nested_site <- df_site %>%
            filter(!is.na(log2inty)) %>%
            select(protein, site, batch, run, feature, log2inty) %>%
            group_by(protein, site, batch) %>%
            nest()
    } else {
        nested_site <- df_site %>%
            filter(!is.na(log2inty)) %>%
            select(protein, site, run, feature, log2inty) %>%
            group_by(protein, site) %>%
            nest()
    }

    if (!impute_missing) {
        # No missing value imputation
        nested_site <- nested_site %>%
            mutate(sumdata = map(data, summarize_feature, method))
    } else {
        nested_site <- nested_site %>%
            mutate(data = map(data, ~ complete(., feature, run)))
        # Missing value imputation with AFT model
        nested_site <- nested_site %>%
            mutate(aftdata = map(data, fill_censored_aft))
        nested_site <- nested_site %>%
            mutate(sumdata = map(aftdata, summarize_feature, method))
    }
    # Run-level summaries with grouping information
    df_sum <- nested_site %>%
        unnest(sumdata) %>%
        left_join(design)

    return(df_sum)
}
