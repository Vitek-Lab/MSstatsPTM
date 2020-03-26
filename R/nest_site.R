#' Whole-plot modeling.
#'
#' \code{nest_site} fits and returns a whole-plot model for each site in
#'   consideration of the experimental design, organized in a nested data frame.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate group_by n_distinct
#' @importFrom tidyr nest
#' @importFrom purrr map map2 map_dbl map_lgl map_int possibly
#' @param df_sum A data frame.
#' @param w_batch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return An instance of nested data frame.
#' @export
#'
#' @examples
#' nest_site(df_sum, w_batch = FALSE)
nest_site <- function(df_sum, w_batch = FALSE) {
    if (w_batch) {
        if (!("batch" %in% names(df_sum)))
            stop("There is no information about batch!")

        # One model for all batches (need data with >1 batches)
        nested_sum <- df_sum %>%
            nest(data = -one_of("protein", "site")) %>%
            mutate(nb_bch = map_int(data, ~ n_distinct(.$batch))) %>%
            filter(nb_bch > 1) %>%
            select(-nb_bch)
        # nested_sum <- df_sum %>%
        #     group_by(protein, site) %>%
        #     nest_legacy() %>%
        #     mutate(nb_bch = map_int(data, ~ n_distinct(.$batch))) %>%
        #     filter(nb_bch > 1) %>%
        #     select(-nb_bch)

        # Linear model and associated parameter estimates
        nested_sum <- nested_sum %>%
            mutate(lm_fit = map(data, possibly(lm_group, otherwise = NA_real_), w_batch)) %>%
            filter(!is.na(lm_fit)) %>%
            mutate(
                param = map2(lm_fit, data, tidy_bch),
                df_res = map_dbl(lm_fit, df.residual)
            )
        # nested_sum <- nested_sum %>%
        #     mutate(lm_fit = map(data, lm_group, w_batch)) %>%
        #     mutate(
        #         param = map2(lm_fit, data, tidy_bch),
        #         df_res = map_dbl(lm_fit, df.residual)
        #     )
        # Remove cases not eligible for hypothesis testing (SE is NA)
        # [TODO]: consider another option to report FC alone
        nested_sum <- nested_sum %>%
            filter(df_res > 0, !map_lgl(param, ~any(is.na(.$std.error))))
    } else {
        # One model per site (and potentially batch)
        if ("batch" %in% names(df_sum)) {
            nested_sum <- df_sum %>%
                nest(data = -one_of("protein", "site", "batch"))
            # nested_sum <- df_sum %>%
            #     group_by(protein, site, batch) %>%
            #     nest_legacy()
        } else {
            nested_sum <- df_sum %>%
                nest(data = -one_of("protein", "site"))
            # nested_sum <- df_sum %>%
            #     group_by(protein, site) %>%
            #     nest_legacy()
        }
        # Linear model and associated parameter estimates
        nested_sum <- nested_sum %>%
            mutate(lm_fit = map(data, possibly(lm_group, otherwise = NA_real_), w_batch)) %>%
            filter(!is.na(lm_fit)) %>%
            mutate(
                param = map2(lm_fit, data, tidy_bch),
                df_res = map_dbl(lm_fit, df.residual)
            )
        # nested_sum <- nested_sum %>%
        #     mutate(lm_fit = map(data, lm_group, w_batch)) %>%
        #     mutate(
        #         param = map2(lm_fit, data, tidy_bch),
        #         df_res = map_dbl(lm_fit, df.residual)
        #     )
        # Remove cases not eligible for hypothesis testing (SE is NA)
        # [TODO]: consider another option to report FC alone
        nested_sum <- nested_sum %>%
            filter(df_res > 0, !map_lgl(param, ~any(is.na(.$std.error))))
    }

    return(nested_sum)
}
