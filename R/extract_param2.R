#' Extract estimated parameters.
#'
#' \code{extract_param2} extracts the estimated model parameters from nested data
#'   frame.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom tidyr unnest
#' @param nested An instance of nested data frame.
#' @return A data frame restoring the estimated model parameters.
#' @export
#'
#' @examples
#' extract_param(nested)
extract_param2 <- function(nested) {
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
        unnest(param)

    return(param_mod)
}
