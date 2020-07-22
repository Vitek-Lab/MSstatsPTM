#' Site-level summarization
#'
#' \code{PTMsummarize} summarizes the peak log2-intensities for each PTM site
#' into one value per run. If protein peak-intensities are availble, the same
#' summarization procedure is applied to each protein as well.
#'
#' @param data A list of two data frames named \code{PTM} and \code{Protein}.
#'   The \code{PTM} data frame includes columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{feature}, \code{log2inty}, and possibly,
#'   \code{batch}. The \code{Protein} data frame includes all columns as in
#'   \code{PTM} except \code{site}.
#' @param method A string defining the summarization method. Default is
#'   \code{"tmp"}, which applies Tukey's median polish. Other methods include
#'   log2 of the summation of peak intensities (\code{"logsum"}), and mean
#'   (\code{"mean"}), median (\code{"median"}) and max (\code{"max"}) of the
#'   log2-intensities.
#'
#' @return A list of two data frames named \code{PTM} and \code{Protein}. The
#'   \code{PTM} data frame has columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{log2inty}, and possibly, \code{batch}. The
#'   \code{Protein} data frame includes all as in \code{PTM}, except \code{site}.
#'
#' @examples
#' sim <- PTMsimulateExperiment(nGroup = 2, nRep = 2, nProtein = 1, nSite = 1, nFeature = 5,
#' list(PTM = 25, Protein = 25), list(PTM = c(0, 1), Protein = c(0, 1)),
#' list(PTM = 0.2, Protein = 0.2), list(PTM = 0.05, Protein = 0.05))
#' PTMsummarize(sim)
#'
#' @export
PTMsummarize <- function(data, method = "tmp") {
    # Check the PTM data
    if (is.null(data[["PTM"]]))
        stop("PTM peak list is missing!")
    if (!is.data.frame(data[["PTM"]]))
        stop(paste0("Provide a data frame of peak log2-intensity for",
                    " each site in each run in ", sQuote("data$PTM")))
    cols_peak <- c("protein", "site", "group", "run", "feature", "log2inty")
    if (!all(cols_peak %in% names(data[["PTM"]]))) {
        stop("Please include in the PTM data frame all the following columns: ",
             paste0(sQuote(cols_peak), collapse = ", "))
    }

    # Check the Protein data
    if (is.null(data[["Protein"]])) {
        wo_prot <- TRUE
    } else {
        wo_prot <- FALSE
        if (!is.data.frame(data[["Protein"]]))
            stop(paste0("Provide a data frame of peak log2-intensity for",
                        " each protein in each run in ", sQuote("data$Protein")))
        cols_prot <- setdiff(cols_peak, "site")
        if (!all(cols_prot %in% names(data[["Protein"]]))) {
            stop("Please include in the protein data frame all the following columns: ",
                 paste0(sQuote(cols_prot), collapse = ", "))
        }
    }

    # Summarize for the PTM data
    df <- data[["PTM"]]
    df <- df[!is.na(df$log2inty), ]

    # Nested data frame with site as the analysis unit
    cols_nested <- setdiff(cols_peak, "group")
    if ("batch" %in% names(df)) {
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
    summ <- left_join(unnest(nested, one_of("res")), design)

    if (wo_prot) {
        res <- list(PTM = summ)
    } else {
        # Summarize for the Protein data
        df_prot <- data[["Protein"]]
        df_prot <- df_prot[!is.na(df_prot$log2inty), ]

        # Nested data frame with protein as the analysis unit
        cols_nested <- setdiff(cols_prot, "group")
        if ("batch" %in% names(df_prot)) {
            nested_prot <- nest(df_prot[, c(cols_nested, "batch")],
                                data = one_of("run", "feature", "log2inty"))
        } else {
            nested_prot <- nest(df_prot[, cols_nested],
                                data = one_of("run", "feature", "log2inty"))
        }

        # Summarize features per protein per MS run
        nested_prot$res <- lapply(nested_prot$data, summarizeFeatures, method)
        nested_prot <- nested_prot[, names(nested_prot) != "data"]
        design_prot <- unique(df_prot[c("run", "group")])  # Experimental design
        summ_prot <- left_join(unnest(nested_prot, one_of("res")), design_prot)

        res <- list(PTM = summ, Protein = summ_prot)
    }
    res
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
#' @examples
#' df <- data.frame(
#'   run = c("a", "a", "a", "b", "b"),
#'   feature = c("F1", "F2", "F3", "F1", "F3"),
#'   log2inty = rnorm(5))
#' summarizeFeatures(df, method = "tmp")
#'
#' @export
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
        method, choices = .summarizeMethods(), several.ok = FALSE
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
.summarizeMethods <- function() {
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
