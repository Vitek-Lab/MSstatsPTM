#' Site-level summarization
#'
#' \code{PTMsummarize} summarizes log2-intensities of spectral features for each
#' PTM site into one value per run.
#'
#' @param df A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{feature}, \code{log2inty}, and possibly,
#'   \code{batch}.
#' @param method A string defining the summarization method. Default is
#'   \code{"tmp"}, which applies Tukey's median polish. Other methods include
#'   log2 of intensity summation (\code{"logsum"}), and mean (\code{"mean"}),
#'   median (\code{"median"}) and max (\code{"max"}) of the log-intensities.
#'
#' @return A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{log2inty}, and possibly, \code{batch}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' PTMsummarize(df)
#' }
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

    # Nested data frame with site as the analysis unit
    cols_nested <- c("protein", "site", "run", "feature", "log2inty")
    if ("batch" %in% cols) {
        nested <- nest(df[, c(cols_nested, "batch")],
                       data = one_of("run", "feature", "log2inty"))
    } else {
        nested <- nest(df[, cols_nested],
                       data = one_of("run", "feature", "log2inty"))
    }

    # Summarize features per site per MS run
    nested$res <- lapply(nested$data, summarizeFeatures, method)
    nested <- nested[, names(nested) != "data"]
    design <- unique(df[c("run", "group")])  # Experimental design
    left_join(unnest(nested, one_of("res")), design)
}

#' Summarization for one site
#'
#' \code{summarizeFeatures} summarizes feature log2-intensities for a PTM site
#' and returns one summarized value per run. Tukey's median polish is used by
#' default.
#'
#' @param df A data frame with columns of \code{run}, \code{feature}, and
#'   \code{log2inty}.
#' @param method A string defining the summarization method. Default is
#'   \code{"tmp"}, which applies Tukey's median polish. Other methods include
#'   log2 of the sum of intensity (\code{"logsum"}), and mean (\code{"mean"}),
#'   median (\code{"median"}) and max (\code{"max"}) of the log2-intensities.
#'
#' @return A tibble restoring one summarized value per MS run.
#'
#' @export
#'
#' @examples
#' df <- data.frame(
#'   run = c("a", "a", "a", "b", "b"),
#'   feature = c("F1", "F2", "F3", "F1", "F3"),
#'   log2inty = rnorm(5))
#' summarizeFeatures(df, method = "tmp")
summarizeFeatures <- function(df, method = "tmp") {
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
    method <- match.arg(
        method, choices = summarizeMethods(), several.ok = FALSE
    )

    if (method == "tmp") {
        res <- .summarize_tmp(df)
    } else if (method == "logsum") {
        res <- .summarize_logsum(df)
    } else if (method == "mean") {
        res <- .summarize_mean(df)
    } else if (method == "median") {
        res <- .summarize_med(df)
    } else if (method == "max") {
        res <- .summarize_max(df)
    }
    res
}

#' @keywords internal
summarizeMethods <- function() {
    c("tmp", "logsum", "mean", "median", "max")
}

#' @keywords internal
.summarize_tmp <- function(df, ...) {
    wd <- tidyr::pivot_wider(df[, c("feature", "run", "log2inty")],
                             names_from = .data$feature, values_from = .data$log2inty)
    m <- data.matrix(wd[, -1])
    res <- stats::medpolish(m, na.rm = TRUE, trace.iter = FALSE)
    tibble(run = wd$run, log2inty = res$overall + res$row)
}

#' @keywords internal
.summarize_logsum <- function(df) {
    by_run <- group_by(df, .data$run)
    summarise(by_run, log2inty = log2(sum(2 ^ .data$log2inty, na.rm = TRUE)))
}

#' @keywords internal
.summarize_mean <- function(df, ...) {
    by_run <- group_by(df, .data$run)
    summarise(by_run, log2inty = mean(.data$log2inty, na.rm = TRUE))
}

#' @keywords internal
.summarize_med <- function(df, ...) {
    by_run <- group_by(df, .data$run)
    summarise(by_run, log2inty = stats::median(.data$log2inty, na.rm = TRUE))
}

#' @keywords internal
.summarize_max <- function(df, ...) {
    by_run <- group_by(df, .data$run)
    summarise(by_run, log2inty = max(.data$log2inty, na.rm = TRUE))
}
