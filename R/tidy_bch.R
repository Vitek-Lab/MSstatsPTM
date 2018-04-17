
#' Extract estimate of group effect.
#'
#' \code{tidy_bch} extracts the estimate of group effect.
#'
#' @param fit An object of class \code{lm}.
#' @param df_fit A number. Degrees of freedom.
#' @return A data frame restoring the estimated model parameters.
#'
#' @examples
#' tidy_bch(fit, df_fit)
#'
#' @export
#'
tidy_bch <- function(fit, df_fit) {
    params <- broom::tidy(fit) %>%
        filter(!stringr::str_detect(term, "batch"))
    if (n_distinct(df_fit$group) == 1) {
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
