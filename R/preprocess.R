#' Normalization of log2-intensities across MS runs
#'
#' \code{PTMnormalize} normalizes log2-intensities of spectral features across
#' MS runs using a reference, or by equalizing a chosen summary (the log2
#' intensity summation, median, or mean of log2-intensities) from all features,
#' features of modified peptides or features of unmodified peptides.
#'
#' @param data A list of two data frames named \code{PTM} and \code{PROTEIN}.
#'   Both the \code{PTM} data frame and the \code{PROTEIN} data frame include
#'   columns of \code{run}, \code{feature}, and \code{log2inty}.
#' @param method A string defining the normalization method. Default is
#'   \code{"median"}, which equalizes the medians of log2-intensities across MS
#'   runs. Other methods include to equalize log2 of intensity summation
#'   (\code{"logsum"}), to equalize the means of log2-intensities
#'   (\code{"mean"}), and to adjust the log2-intensities based on a reference
#'   (\code{"ref"}) given by (\code{refs}).
#' @param refs A list of two data frames named \code{PTM} and \code{PROTEIN}.
#'   Each defines the adjustment of log2-intensities for the MS runs in its
#'   corresponding data.
#'
#' @return Normalized data stored as in \code{data}.
#'
#' @examples
#' sim <- PTMsimulateExperiment(
#'     nGroup=2, nRep=2, nProtein=1, nSite=1, nFeature=5,
#'     list(PTM=25, PROTEIN=25), list(PTM=c(0, 1), PROTEIN=c(0, 1)),
#'     list(PTM=0.2, PROTEIN=0.2), list(PTM=0.05, PROTEIN=0.05)
#' )
#' PTMnormalize(sim)
#'
#' @export
PTMnormalize <- function(data, method="median", refs) {
    cols_req <- c("run", "feature", "log2inty")
    norm_opt <- c("median", "mean", "logsum", "ref")

    # Check the PTM data
    if (is.null(data[["PTM"]]))
        stop("PTM peak list is missing!")
    if (!is.data.frame(data[["PTM"]]))
        stop(paste0("Provide a data frame of peak log2-intensity for",
                    " the PTM data in ", sQuote("data$PTM")))
    if (!all(cols_req %in% names(data[["PTM"]]))) {
        stop(
            "Include in the PTM data frame the following columns: ",
            paste0(sQuote(cols_req), collapse = ", ")
        )
    }

    # Check the PROTEIN data
    if (is.null(data[["PROTEIN"]])) {
        wo_prot <- TRUE
    } else {
        wo_prot <- FALSE
        if (!is.data.frame(data[["PROTEIN"]]))
            stop(paste0(
                "Provide a data frame of peak log2-intensity for",
                " the PROTEIN data in ", sQuote("data$PROTEIN")
            ))
        if (!all(cols_req %in% names(data[["PROTEIN"]]))) {
            stop(
                "Include in the PROTEIN data frame the following columns: ",
                paste0(sQuote(cols_req), collapse = ", ")
            )
        }
    }

    # Check input method
    if (!is.character(method) || length(method) != 1) {
        stop("Define normalization method as a string in ", sQuote("method"))
    }
    if (!(method %in% norm_opt)) {
        stop(
            "Define the normalization method as one of the following: ",
            paste0(sQuote(norm_opt), collapse = ", ")
        )
    }

    if (method == "ref") {
        # Check and obtain reference
        if (missing(refs))
            stop("Define the adjustment in ", sQuote("refs"))

        # Check the reference for the PTM data
        if (!is.data.frame(refs[["PTM"]])) {
            stop(
                "Define the adjustment for PTM data as a data frame in ",
                sQuote("refs$PTM")
            )
        }
        if (!all(c("run", "adjLog2inty") %in% names(refs[["PTM"]]))) {
            stop(
                "Please include in ", sQuote("refs$PTM"),
                " the following columns: ", sQuote("run"),
                sQuote("adjLog2inty")
            )
        }
        if (!all(unique(data[["PTM"]]$run) %in% refs[["PTM"]]$run)) {
            stop("Adjustment is not fully defined for all MS runs!")
        }
        ref_PTM <- refs[["PTM"]]

        if (!wo_prot) {
            # Check the reference for the PROTEIN data
            if (!is.data.frame(refs[["PROTEIN"]])) {
                stop(
                    "Define adjustment for PROTEIN data as a data frame in ",
                    sQuote("refs$PROTEIN")
                )
            }
            if (!all(c("run", "adjLog2inty") %in% names(refs[["PROTEIN"]]))) {
                stop(
                    "Please include in ", sQuote("refs$PROTEIN"),
                    " the following columns: ", sQuote("run"),
                    sQuote("adjLog2inty")
                )
            }
            wa <- all(unique(data[["PROTEIN"]]$run) %in% refs[["PROTEIN"]]$run)
            if (!wa) {
                stop("Adjustment is not fully defined for all MS runs!")
            }
            ref_prot <- refs[["PROTEIN"]]
        }
    } else {
        # Compute reference
        ref_PTM <- .getReference(data[["PTM"]], method)
        if (!wo_prot) {
            ref_prot <- .getReference(data[["PROTEIN"]], method)
        }
    }

    # Adjust based on the reference
    if (wo_prot) {
        res <- list(PTM = .byReference(data[["PTM"]], ref_PTM))
    } else {
        res <- list(
            PTM = .byReference(data[["PTM"]], ref_PTM),
            PROTEIN = .byReference(data[["PROTEIN"]], ref_prot)
        )
    }
    res
}


#' @keywords internal
.getReference <- function(df, method="median") {
    g <- group_by(df, .data$run)
    if (method == "median") {
        gs <- summarise(g, log2inty = stats::median(.data$log2inty, na.rm=TRUE))
    } else if (method == "mean") {
        gs <- summarise(g, log2inty = mean(.data$log2inty, na.rm=TRUE))
    } else {
        gs <- summarise(g, log2inty = log2(sum(2 ^ .data$log2inty, na.rm=TRUE)))
    }
    gbl <- stats::median(gs$log2inty)
    tibble(run = gs$run, adjLog2inty = gbl - gs$log2inty)
}

#' @keywords internal
.byReference <- function(df, ref) {
    df_aug <- left_join(df, ref)
    df_aug$log2inty <- df_aug$log2inty + df_aug$adjLog2inty
    df_aug[, names(df_aug) != "adjLog2inty"]
}
