#' Normalization of log2-intensities across MS runs.
#'
#' \code{PTMnormalize} normalizes log2-intensities of spectral features across
#' MS runs.
#'
#' @param df A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{feature}, \code{log2inty}, and possibly,
#'   \code{batch}.
#' @param method A string defining the normalization method. Default is
#'   \code{"median"}, which equalizes the medians of log2-intensities across MS
#'   runs. Other methods include to equalize log2 of intensity summation
#'   (\code{"logsum"}), to equalize the means of log2-intensities
#'   (\code{"mean"}), and to adjust the log2-intensities based on a reference
#'   (\code{"ref"}) given by (\code{reference}).
#' @param reference A data frame defining the adjustment of log2-intensities for
#' each MS runs, with columns of \code{run} and \code{adjLog2inty}
#'
#' @return A data frame with same columns as in \code{df}.
#'
#' @export
#'
#' @examples
#' PTMsummarize(df)
PTMsummarize <- function(df, method = "median", reference) {
    if (missing(df))
        stop("The input ", sQuote("df"), " is missing!")
    if (!is.data.frame(df))
        stop("Provide the peak list as a data frame in ", sQuote("df"))
    cols_peak <- c("protein", "site", "group", "run", "feature", "log2inty")
    cols <- names(df)
    if (!all(cols_peak %in% cols)) {
        stop("Please include in the data frame all the following columns: ",
             paste0(sQuote(cols_peak), collapse = ", "))
    }
    norm_meth <- c("median", "mean", "logsum", "ref")
    if (!is.character(method) || length(method) != 1) {
        stop("Define normalization method as a string in ", sQuote("method"))
    } else if (!(method %in% norm_meth)) {
        stop("Define the normalization method as one of the following: ",
             paste0(sQuote(norm_meth), collapse = ", "))
    }
    if (method == "ref") {
        if (missing(reference) || !is.data.frame(reference)) {
            stop("Define the adjustment as a data frame in ", sQuote("reference"))
        } else if (!all(c("run", "adjLog2inty") %in% names(reference))) {
            stop("Please include in ", sQuote("reference"),
                 " the following columns: ", sQuote("run"), sQuote("adjLog2inty"))
        } else if (!all(unique(df$run) %in% reference$run)) {
            stop("Adjustment is not fully defined for all MS runs!")
        }
    }

    if (method != "ref") {
        r <- group_by(df, run)
        if (method == "median") {
            s <- summarise(r, log2inty = stats::median(log2inty, na.rm = TRUE))
        } else if (method == "mean") {
            s <- summarise(r, log2inty = mean(log2inty, na.rm = TRUE))
        } else {
            s <- summarise(r, log2inty = log2(sum(2 ^ log2inty, na.rm = TRUE)))
        }
        ref <- stats::median(s$log2inty)
        reference <- tibble(run = s$run, adjLog2inty = ref - s$log2inty)
    }

    df_aug <- left_join(df, reference)
    df_aug$log2inty <- df_aug$log2inty + df_aug$adjLog2inty

    df_aug[, names(df_aug) != "adjLog2inty"]
}


#' Normalize feature intensities across runs.
#'
#' \code{normalize_ptm} normalizes feature intensities across runs in a PTM
#'   dataset, using unpaired unmodified peptides.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate group_by summarise ungroup left_join
#' @param data A data frame.
#' @return A data frame.
#' @export
#'
#' @examples
#' normalize_ptm(data)
normalize_ptm <- function(data) {

    if (!is.data.frame(data))
        stop("Input data should be a data frame!")
    cols <- names(data)
    if (!all(c("is_mod", "run", "log2inty") %in% cols))
        stop("Input data frame should contain columns of is_mod, run, log2inty!")

    # Based on unmodified peptides
    if ("batch" %in% cols) {
        medians <- data %>%
            filter(!is_mod) %>%
            group_by(batch, run) %>%
            summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>%
            mutate(log2inty_bch = median(log2inty_med)) %>%
            ungroup()
    } else {
        medians <- data %>%
            filter(!is_mod) %>%
            group_by(run) %>%
            summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>%
            mutate(log2inty_bch = median(log2inty_med)) %>%
            ungroup()
    }
    normdata <- data %>%
        left_join(medians) %>%
        mutate(log2inty = ifelse(is.na(log2inty_med), log2inty,
                                 log2inty - log2inty_med + log2inty_bch)) %>%
        select(-log2inty_med, -log2inty_bch)

    return(normdata)
}
