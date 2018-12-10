
#' Extract estimated parameters.
#'
#' \code{extract_param2} extracts the estimated model parameters from nested data
#'   frame.
#'
#' @param nested An instance of nested data frame.
#' @return A data frame restoring the estimated model parameters.
#'
#' @examples
#' extract_param(nested)
#'
#' @export
#'
extract_param2 <- function(nested) {
    if ("batch" %in% names(nested)) {
        # per-batch model
        nested_param <- nested %>%
            dplyr::select(protein, site, batch, param, df_res)
    } else {
        # all-batch model
        nested_param <- nested %>%
            dplyr::select(protein, site, param, df_res)
    }
    param_mod <- nested_param %>%
        tidyr::unnest(param)

    return(param_mod)

}
