#' Normalization of log2-intensities across MS runs
#'
#' \code{PTMnormalize} normalizes log2-intensities of spectral features across
#' MS runs.
#'
#' @param df A data frame contains all the following columns: \code{run},
#'   \code{feature}, \code{is_mod}, and \code{log2inty}.
#' @param method A string defining the normalization method. Default is
#'   \code{"median"}, which equalizes the medians of log2-intensities across MS
#'   runs. Other methods include to equalize log2 of intensity summation
#'   (\code{"logsum"}), to equalize the means of log2-intensities
#'   (\code{"mean"}), and to adjust the log2-intensities based on a reference
#'   (\code{"ref"}) given by (\code{reference}).
#' @param calc A string defining parts of log2-intensities used for calculating
#'   the adjustment value for each MS run. Default is \code{"all"}, which uses
#'   all measurements. Other options include to use only log2-intensities from
#'   peptides with PTMs of interest (\code{"mod"}) and log2-intensities from
#'   unmodified peptides (\code{"unmod"}).
#' @param reference A data frame defining the adjustment of log2-intensities for
#' each MS runs, with columns of \code{run} and \code{adjLog2inty}.
#'
#' @return A data frame with same columns as in \code{df}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' PTMnormalize(df)
#' }
PTMnormalize <- function(df, method = "median", calc = "all", reference) {
    cols_req <- c("run", "feature", "is_mod", "log2inty")
    norm_opt <- c("median", "mean", "logsum", "ref")
    calc_opt <- c("all", "mod", "unmod")

    if (missing(df))
        stop("The input ", sQuote("df"), " is missing!")
    if (!is.data.frame(df))
        stop("Provide the peak list as a data frame in ", sQuote("df"))
    if (!all(cols_req %in% names(df))) {
        stop("Please include in the data frame all the following columns: ",
             paste0(sQuote(cols_req), collapse = ", "))
    }

    if (!is.character(method) || length(method) != 1) {
        stop("Define normalization method as a string in ", sQuote("method"))
    }
    if (!(method %in% norm_opt)) {
        stop("Define the normalization method as one of the following: ",
             paste0(sQuote(norm_opt), collapse = ", "))
    }

    if (method != "ref") {
        if (missing(calc) || !is.character(calc) || length(calc) != 1) {
            stop("Define calculation option as a string in ", sQuote("calc"))
        }
        if (!(calc %in% calc_opt)) {
            stop("Define the calculation option as one of the following: ",
                 paste0(sQuote(calc_opt), collapse = ", "))
        }
    }

    if (method == "ref") {
        if (missing(reference) || !is.data.frame(reference)) {
            stop("Define the adjustment as a data frame in ", sQuote("reference"))
        }
        if (!all(c("run", "adjLog2inty") %in% names(reference))) {
            stop("Please include in ", sQuote("reference"),
                 " the following columns: ", sQuote("run"), sQuote("adjLog2inty"))
        }
        if (!all(unique(df$run) %in% reference$run)) {
            stop("Adjustment is not fully defined for all MS runs!")
        }
    }

    if (method != "ref") {
        if (calc == "mod") {
            s <- df[df$is_mod, ]
        } else if (calc == "unmod") {
            s <- df[!df$is_mod, ]
        } else {
            s <- df
        }

        g <- group_by(df, .data$run)
        if (method == "median") {
            gs <- summarise(g, log2inty = stats::median(.data$log2inty, na.rm = TRUE))
        } else if (method == "mean") {
            gs <- summarise(g, log2inty = mean(.data$log2inty, na.rm = TRUE))
        } else {
            gs <- summarise(g, log2inty = log2(sum(2 ^ .data$log2inty, na.rm = TRUE)))
        }
        ref <- stats::median(gs$log2inty)
        reference <- tibble(run = gs$run, adjLog2inty = ref - gs$log2inty)
    }

    df_aug <- left_join(df, reference)
    df_aug$log2inty <- df_aug$log2inty + df_aug$adjLog2inty

    df_aug[, names(df_aug) != "adjLog2inty"]
}

