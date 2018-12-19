#' Extract estimate of group effect.
#'
#' \code{tidy_bch} extracts the estimate of group effect.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate n_distinct
#' @importFrom stringr str_detect str_replace
#' @param fit An object of class \code{lm}.
#' @param df_fit A number. Degrees of freedom.
#' @return A data frame restoring the estimated model parameters.
#' @export
#'
#' @examples
#' tidy_bch(fit, df_fit)
tidy_bch <- function(fit, df_fit) {
    params <- broom::tidy(fit) %>%
        filter(!str_detect(term, "batch"))
    if (n_distinct(df_fit$group) == 1) {
        params <- params %>%
            mutate(group = unique(df_fit$group)) %>%
            select(-term, -statistic, -p.value)
    } else {
        params <- params %>%
            mutate(group = str_replace(term, "group", "")) %>%
            select(-term, -statistic, -p.value)
    }

    return(params)
}
