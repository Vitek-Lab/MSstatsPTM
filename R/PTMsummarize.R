#' Summarize feature log-intensities for each PTM site per run.
#'
#' \code{PTMsummarize} summarizes feature log-intensities for each PTM site.
#'
#' @importFrom tidyselect one_of
#' @param df A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{feature}, \code{log2inty}, and possibly,
#'   \code{batch}.
#' @param method A string defining the summarization method. Default is
#'   \code{"tmp"}, which applies Tukey's median polish.
#' @return A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{log2inty}, and possibly, \code{batch}.
#' @export
#' @rdname summarization
#'
#' @examples
#' PTMsummarize(df)
PTMsummarize <- function(df, method = "tmp") {
    if (missing(df))
        stop(paste0("The input ", sQuote("df"), " is missing!"))
    if (!is.data.frame(df))
        stop(paste0("Provide the peak list as a data frame in ", sQuote("df")))
    cols_peak <- c("protein", "site", "group", "run", "feature", "log2inty")
    cols <- names(df)
    if (!all(cols_peak %in% cols)) {
        stop("Please include in the data frame all the following columns: ",
             paste0(sQuote(cols_peak), collapse = ", "))
    }
    df <- df[!is.na(df$log2inty), ]

    # Experimental design
    design <- unique(df[c("run", "group")])

    # Nested data frame
    cols_nested <- c("protein", "site", "run", "feature", "log2inty")
    if ("batch" %in% cols) {
        nested <-  tidyr::nest(
            df[, c(cols_nested, "batch")],
            data = one_of("run", "feature", "log2inty")
        )
    } else {
        nested <-  tidyr::nest(
            df[, cols_nested], data = one_of("run", "feature", "log2inty")
        )
    }

    # Run-level summaries with grouping information
    nested$res <- lapply(nested$data, summarizeFeatures, method)
    nested <- nested[, names(nested) != "data"]
    dplyr::left_join(tidyr::unnest(nested, one_of("res")), design)
}

#' Summarization of feature intensities.
#'
#' \code{summarizeFeatures} summarizes feature intensities and returns one
#'   summarized value per run. Tukey's median polish is used by default.
#'
#' @param df A data frame.
#' @param method A string. Default is \code{"tmp"}.
#' @return A tibble restoring the run-level summarization.
#' @export
#' @rdname summarization
#'
#' @examples
#' df <- data.frame(
#'   run = c("a", "a", "a", "b", "b"),
#'   feature = c("F1", "F2", "F3", "F1", "F3"),
#'   log2inty = rnorm(5))
#' summarizeFeatures(df, method = "tmp")
summarizeFeatures <- function(df, method = "tmp", ...) {
    if (missing(df))
        stop(paste0("The input ", sQuote("df"), " is missing!"))
    if (!is.data.frame(df))
        stop(paste0("Provide the peak list as a data frame in ", sQuote("df")))
    if (!all(c("run", "feature", "log2inty") %in% names(df))) {
        stop("Missing columns! Please make sure the input data frame contains
             columns of 'run', 'feature', 'log2inty'")
    }
    if (missing(method))
        stop("Please specify a summarization method.")
    method <- match.arg(method,
                        choices = summarizeMethods(),
                        several.ok = FALSE)

    if (method == "tmp") {
        res <- summarize_tmp(df, ...)
    } else if (method == "logsum") {
        res <- summarize_logsum(df, ...)
    } else if (method == "mean") {
        res <- summarize_mean(df, ...)
    } else if (method == "median") {
        res <- summarize_med(df, ...)
    } else if (method == "max") {
        res <- summarize_max(df, ...)
    }
    res
}

#' @rdname summarization
summarizeMethods <- function() {
    c("tmp", "logsum", "mean", "median", "max")
}

#' @rdname summarization
summarize_tmp <- function(df, ...) {
    wd <- tidyr::pivot_wider(df[, c("feature", "run", "log2inty")],
                             names_from = feature, values_from = log2inty)
    m <- data.matrix(wd[, -1])
    res <- stats::medpolish(m, na.rm = TRUE, trace.iter = FALSE)
    dplyr::tibble(run = wd$run, log2inty = res$overall + res$row)
}

#' @rdname summarization
summarize_logsum <- function(df, ...) {
    by_run <- dplyr::group_by(df, run)
    dplyr::summarise(by_run, log2inty = log2(sum(2 ^ log2inty, na.rm = TRUE)))
}

#' @rdname summarization
summarize_mean <- function(df, ...) {
    by_run <- dplyr::group_by(df, run)
    dplyr::summarise(by_run, log2inty = mean(log2inty, na.rm = TRUE))
}

#' @rdname summarization
summarize_med <- function(df, ...) {
    by_run <- dplyr::group_by(df, run)
    dplyr::summarise(by_run, log2inty = median(log2inty, na.rm = TRUE))
}

#' @rdname summarization
summarize_max <- function(df, ...) {
    by_run <- dplyr::group_by(df, run)
    dplyr::summarise(by_run, log2inty = max(log2inty, na.rm = TRUE))
}
