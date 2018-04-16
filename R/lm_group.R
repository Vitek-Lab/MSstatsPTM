
#' Fit a linear model with group effect (and batch effect if present).
#'
#' \code{lm_group} fits and returns a linear model with \code{group} effect.
#'
#' @param df_onesite A data frame containing columns \code{log2inty},
#'   \code{group}.
#' @param w_batch A logical scalar. \code{TRUE} includes batch effect,
#'   \code{FALSE} otherwise. Default is \code{TRUE}.
#' @return An object of class \code{lm}.
#'
#' @examples
#' lm_group(df_onesite)
#'
#' @export
#'
lm_group <- function(df_onesite, w_batch = FALSE) {
    if (w_batch) {
        # Model with batch effect
        if (n_distinct(df_onesite$batch) == 1)
            stop("Cannot estimate batch effect with a single batch!")

        if (n_distinct(df_onesite$group) == 1) {
            fit <- lm(log2inty ~ batch, data = df_onesite)
        } else {
            fit <- lm(log2inty ~ 0 + group + batch, data = df_onesite)
        }
    } else {
        if (n_distinct(df_onesite$group) == 1) {
            fit <- lm(log2inty ~ 1, data = df_onesite)
        } else {
            fit <- lm(log2inty ~ 0 + group, data = df_onesite)
        }
    }

    return(fit)
}
