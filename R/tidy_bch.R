
#' Extract estimate of group effect.
#'
#' \code{tidy_bch} extracts the estimate of group effect.
#'
#' @param bch_fit An object of class \code{lm}.
#' @param df_bch A number. Degrees of freedom.
#' @return A data frame restoring the estimated model parameters.
#'
#' @examples
#' tidy_bch(bch_fit, df_bch)
#'
#' @export
#'
tidy_bch <- function(bch_fit, df_bch) {
    params <- broom::tidy(bch_fit) %>%
        filter(!stringr::str_detect(term, "batch"))
    if (n_distinct(df_bch$group) == 1) {
        params <- params %>%
            dplyr::mutate(group = unique(df_bch$group)) %>%
            dplyr::select(-term, -statistic, -p.value)
    } else {
        params <- params %>%
            dplyr::mutate(group = stringr::str_replace(term, "group", "")) %>%
            dplyr::select(-term, -statistic, -p.value)
    }

    return(params)
}
