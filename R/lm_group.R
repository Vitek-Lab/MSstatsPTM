
#' Fit a linear model with group effect (and batch effect if present).
#'
#' \code{lm_group} fits and returns a linear model with \code{group} effect.
#'
#' @param onesite A data frame containing columns \code{log2inty},
#'   \code{group} for one modification site.
#' @param w_batch A logical scalar. \code{TRUE} includes batch effect,
#'   \code{FALSE} otherwise. Default is \code{FALSE}.
#' @return An object of class \code{lm}.
#'
#' @examples
#' x <- data.frame(group = rep(c("1", "2"), 3), log2inty = rep(c(10, 12), 3) + rnorm(6))
#' lm_group(x)
#'
#' @export
#'
lm_group <- function(onesite, w_batch = FALSE) {
    if (w_batch) {
        # Model with batch effect
        if (n_distinct(onesite$batch) == 1)
            stop("Cannot estimate batch effect with a single batch!")

        if (n_distinct(onesite$group) == 1) {
            fit <- lm(log2inty ~ batch, data = onesite)
        } else {
            fit <- lm(log2inty ~ 0 + group + batch, data = onesite)
        }
    } else {
        if (n_distinct(onesite$group) == 1) {
            fit <- lm(log2inty ~ 1, data = onesite)
        } else {
            fit <- lm(log2inty ~ 0 + group, data = onesite)
        }
    }

    return(fit)
}
