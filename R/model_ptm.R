#' Whole-plot modeling for all sites.
#'
#' \code{model_ptm} fits a whole-plot model for all sites in consideration of
#'   the experimental design and the underlying protein abundance, and returns
#'   the estimates of model parameters.
#'
#' @param df_sum A data frame.
#' @param w_batch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return A data frame.
#' @export
#'
#' @examples
#' model_ptm(df_sum, w_batch = FALSE)
model_ptm <- function(df_sum, w_batch = FALSE) {
    nested <- nest_site(df_sum, w_batch)
    params <- extract_param(nested)

    return(params)
}
