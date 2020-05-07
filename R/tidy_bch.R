#' Extract estimate of group effect.
#'
#' \code{tidy_bch} extracts the estimate of group effect.
#'
#' @importFrom magrittr %>%
#' @param fit An object of class \code{lm}.
#' @param df_fit A number. Degrees of freedom.
#' @return A data frame restoring the estimated model parameters.
#' @export
#'
#' @examples
#' tidy_bch(fit, df_fit)
tidy_bch <- function(fit, df_fit) {
    params <- broom::tidy(fit) %>%
        dplyr::filter(!stringr::str_detect(term, "batch"))
    if (dplyr::n_distinct(df_fit$group) == 1) {
        params <- params %>%
            dplyr::mutate(group = unique(df_fit$group)) %>%
            dplyr::select(-term, -statistic, -p.value)
    } else {
        params <- params %>%
            dplyr::mutate(group = stringr::str_replace(term, "group", "")) %>%
            dplyr::select(-term, -statistic, -p.value)
    }

    return(params)
}
