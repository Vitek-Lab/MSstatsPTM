#' Extract estimated parameters.
#'
#' \code{extract_param} extracts the estimated model parameters from nested data
#'   frame.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select rename left_join
#' @importFrom tidyr unnest
#' @param nested An instance of nested data frame.
#' @return A data frame restoring the estimated model parameters.
#' @export
#'
#' @examples
#' extract_param(nested)
extract_param <- function(nested) {
    if ("batch" %in% names(nested)) {
        # per-batch model
        nested_param <- nested %>%
            select(protein, site, batch, param, df_res)
    } else {
        # all-batch model
        nested_param <- nested %>%
            select(protein, site, param, df_res)
    }
    param_mod <- nested_param %>%
        filter(site != "None") %>%
        unnest(param)
    param_unmod <- nested_param %>%
        filter(site == "None") %>%
        select(-site) %>%
        unnest(param) %>%
        rename(df_unmod = df_res, est_unmod = estimate,
                      se_unmod = std.error)

    return(left_join(param_mod, param_unmod))
}
