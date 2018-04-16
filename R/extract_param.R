
#' Extract estimated parameters.
#'
#' \code{extract_param} extracts the estimated model parameters.
#'
#' @param nested_data An instance of nested data frame.
#' @return A data frame restoring the estimated model parameters.
#'
#' @examples
#' extract_param(nested_data)
#'
#' @export
#'
extract_param <- function(nested_data) {
    if ("batch" %in% names(nested_data)) {
        # per-batch model
        nested_param <- nested_data %>%
            dplyr::select(uniprot_iso, site_str, batch, param, df_res)
    } else {
        # all-batch model
        nested_param <- nested_data %>%
            dplyr::select(uniprot_iso, site_str, param, df_res)
    }
    param_mod <- nested_param %>%
        dplyr::filter(site_str != "None") %>%
        tidyr::unnest(param)
    param_unmod <- nested_param %>%
        dplyr::filter(site_str == "None") %>%
        dplyr::select(-site_str) %>%
        tidyr::unnest(param) %>%
        dplyr::rename(df_unmod = df_res, est_unmod = estimate,
                      se_unmod = std.error)

    return(dplyr::left_join(param_mod, param_unmod) %>%
               tidyr::unite(protsite, uniprot_iso, site_str, sep = "--"))
}
