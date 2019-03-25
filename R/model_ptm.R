#' Whole-plot modeling for all sites.
#'
#' \code{model_ptm} fits a whole-plot model for all sites in consideration of
#'   the experimental design and the underlying protein abundance, and returns
#'   the estimates of model parameters.
#'
#' @param df_sum A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{log2inty}, and possibly, \code{batch}.
#' @param w_batch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return A data frame.
#' @export
#'
#' @examples
#' model_ptm(df_sum, w_batch = FALSE)
model_ptm <- function(df_sum, w_batch = FALSE) {

    if (missing(df_sum))
        stop("Input data frame is missing!")
    if (!is.data.frame(df_sum))
        stop("Please provide protein summaries as a data frame in 'df_sum'!")
    cols <- names(df_sum)
    if (!all(c("protein", "site", "group", "run", "log2inty") %in% cols)) {
        stop("Missing columns! Please make sure the input data frame contains
             columns of 'protein', 'site', 'group', 'run', 'log2inty'")
    }
    if (w_batch && !("batch" %in% cols))
        stop("To account for batch effect, add a batch column in df_sum")

    nested <- nest_site(df_sum, w_batch)
    params <- extract_param(nested)

    return(params)
}
