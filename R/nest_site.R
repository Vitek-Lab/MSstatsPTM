
#' Whole-plot modeling.
#'
#' \code{nest_site} fits and returns a whole-plot model for each site in
#'   consideration of the experimental design, organized in a nested data frame.
#'
#' @param df_sum A data frame.
#' @param w_batch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return An instance of nested data frame.
#'
#' @examples
#' nest_site(df_sum, w_batch = FALSE)
#'
#' @export
#'
nest_site <- function(df_sum, w_batch = FALSE) {
    if (w_batch) {
        if (!("batch" %in% names(df_sum)))
            stop("There is no information about batch!")

        # One model for all batches (need data with >1 batches)
        nested_sum <- df_sum %>%
            dplyr::group_by(protein, site) %>%
            tidyr::nest() %>%
            dplyr::mutate(nb_bch = purrr::map_int(data, ~ n_distinct(.$batch))) %>%
            dplyr::filter(nb_bch > 1) %>%
            dplyr::select(-nb_bch)
        # Linear model and associated parameter estimates
        nested_sum <- nested_sum %>%
            dplyr::mutate(lm_fit = purrr::map(data, lm_group, w_batch)) %>%
            dplyr::mutate(
                param = purrr::map2(lm_fit, data, tidy_bch),
                df_res = purrr::map_dbl(lm_fit, df.residual)
            )
        # Remove cases not eligible for hypothesis testing (SE is NA)
        # [TODO]: consider another option to report FC alone
        nested_sum <- nested_sum %>%
            dplyr::filter(df_res > 0, !purrr::map_lgl(param, ~any(is.na(.$std.error))))
    } else {
        # One model per site (and potentially batch)
        if ("batch" %in% names(df_sum)) {
            nested_sum <- df_sum %>%
                dplyr::group_by(protein, site, batch) %>%
                tidyr::nest()
        } else {
            nested_sum <- df_sum %>%
                dplyr::group_by(protein, site) %>%
                tidyr::nest()
        }
        # Linear model and associated parameter estimates
        nested_sum <- nested_sum %>%
            dplyr::mutate(lm_fit = purrr::map(data, lm_group, w_batch)) %>%
            dplyr::mutate(
                param = purrr::map2(lm_fit, data, tidy_bch),
                df_res = purrr::map_dbl(lm_fit, df.residual)
            )
        # Remove cases not eligible for hypothesis testing (SE is NA)
        # [TODO]: consider another option to report FC alone
        nested_sum <- nested_sum %>%
            dplyr::filter(df_res > 0, !purrr::map_lgl(param, ~any(is.na(.$std.error))))
    }

    return(nested_sum)
}
