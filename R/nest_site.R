
#' Whole-plot modeling.
#'
#' \code{nest_site} fits and returns a whole-plot model considering experimental
#'   design.
#'
#' @param df_sum A string. Peptide sequence.
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
            group_by(uniprot_iso, site_str) %>%
            nest() %>%
            mutate(nb_bch = map_int(data, ~ n_distinct(.$batch))) %>%
            filter(nb_bch > 1) %>%
            select(-nb_bch)
        # Linear model and associated parameter estimates
        nested_sum <- nested_sum %>%
            mutate(lm_fit = map(data, lm_group, w_batch)) %>%
            mutate(
                param = map2(lm_fit, data, tidy_bch),
                df_res = map_dbl(lm_fit, df.residual)
            )
        # Remove cases not eligible for hypothesis testing (SE is NA)
        # [TODO]: consider another option to report FC alone
        nested_sum <- nested_sum %>%
            filter(df_res > 0, !map_lgl(param, ~any(is.na(.$std.error))))
    } else {
        # One model per site (and potentially batch)
        if ("batch" %in% names(df_sum)) {
            nested_sum <- df_sum %>%
                group_by(uniprot_iso, site_str, batch) %>%
                nest()
        } else {
            nested_sum <- df_sum %>%
                group_by(uniprot_iso, site_str) %>%
                nest()
        }
        # Linear model and associated parameter estimates
        nested_sum <- nested_sum %>%
            mutate(lm_fit = map(data, lm_group, w_batch)) %>%
            mutate(
                param = map2(lm_fit, data, tidy_bch),
                df_res = map_dbl(lm_fit, df.residual)
            )
        # Remove cases not eligible for hypothesis testing (SE is NA)
        # [TODO]: consider another option to report FC alone
        nested_sum <- nested_sum %>%
            filter(df_res > 0, !map_lgl(param, ~any(is.na(.$std.error))))
    }

    return(nested_sum)
}
