#' Estimate log2-abundances of PTM sites and proteins
#'
#' \code{PTMestimate} takes as input the summarized log2-intensities for each
#' PTM site, performs statistical modeling for the abundance of the site, and
#' returns the estimates of model parameters for all sites in all experimental
#' conditions. If protein log2-intensities are availble, the same estimation
#' procedure is applied to each protein as well.
#'
#' @param data A list of two data frames named \code{PTM} and \code{PROTEIN}.
#'   The \code{PTM} data frame includes columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{log2inty}, and possibly, \code{batch}. The
#'   \code{PROTEIN} data frame includes all columns as in \code{PTM} except
#'   \code{site}.
#' @param fctBatch A logical defining the handling of batch effect for all data
#'   or two logicals for the \code{PTM} and \code{PROTEIN} (if provided) data
#'   separately. \code{TRUE} considers a fixed batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return A list of two lists named \code{PTM} and \code{PROTEIN}. The
#'   \code{PTM} list has four elements: \code{protein} (a string vector of
#'   protein names), \code{site} (a string vector of PTM sites), \code{param}
#'   (a list of model parameter estimates for each site), and \code{df} (a
#'   numeric vector of degrees of freedom for each model). The \code{PROTEIN}
#'   list includes all as in \code{PTM}, except \code{site}.
#'
#' @examples
#' sim <- PTMsimulateExperiment(
#'     nGroup=2, nRep=2, nProtein=1, nSite=1, nFeature=5,
#'     logAbundance=list(
#'         PTM=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05),
#'         PROTEIN=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05)
#'     )
#' )
#' s <- PTMsummarize(sim)
#' PTMestimate(s)
#'
#' @export
PTMestimate <- function(data, fctBatch=FALSE) {
    # Check PTM data
    if (is.null(data[["PTM"]]))
        stop("PTM information is missing!")
    if (!is.data.frame(data[["PTM"]]))
        stop("Provide a data frame of summarized values in 'data$PTM'")
    cols_site <- c("protein", "site", "group", "run", "log2inty")
    if (!all(cols_site %in% names(data[["PTM"]]))) {
        stop(
            "Please include in the PTM data frame all the following columns: ",
            paste0(sQuote(cols_site), collapse = ", ")
        )
    }
    # Check PROTEIN data
    if (is.null(data[["PROTEIN"]])) {
        wo_prot <- TRUE
    } else {
        wo_prot <- FALSE
        if (!is.data.frame(data[["PROTEIN"]]))
            stop("Provide a data frame of summarized values in 'data$PROTEIN'")
        cols_prot <- setdiff(cols_site, "site")
        if (!all(cols_prot %in% names(data[["PROTEIN"]]))) {
            stop(
                "Include in the PROTEIN data frame the following columns: ",
                paste0(sQuote(cols_prot), collapse = ", ")
            )
        }
    }
    # Define batch handling
    if (!is.logical(fctBatch))
        stop("'fctBatch' should be logical value(s)")
    if (!(length(fctBatch) %in% c(1L, 2L)))
        stop(paste0(
            "Provide one logical value for all data or two logical",
            " values for PTM and PROTEIN separately in 'fctBatch'"
        ))
    if (length(fctBatch) == 2L && wo_prot)
        stop("'fctBatch' is defined for PROTEIN, but data is missing!")
    if (length(fctBatch) == 1L && !wo_prot) fctBatch <- rep(fctBatch, 2)
    est_ptm <- estimateAbundance(data[["PTM"]], fctBatch[1], perProtein = FALSE)
    if (wo_prot) {
        res <- list(PTM = est_ptm)
    } else {
        est_prot <- estimateAbundance(
            data[["PROTEIN"]], fctBatch[2], perProtein = TRUE
        )
        res <- list(PTM = est_ptm, PROTEIN = est_prot)
    }
    res
}


#' Estimate log2-abundances of PTM sites or proteins
#'
#' \code{estimateAbundance} takes as input the summarized log2-intensities for
#' each PTM site, performs statistical modeling for the abundance of the site,
#' and returns the estimates of model parameters for all sites in all
#' experimental conditions.
#'
#' @param df A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{log2inty}, and possibly, \code{batch}.
#' @param fctBatch A logical. \code{TRUE} considers a fixed batch effect,
#'   \code{FALSE} otherwise. Default is \code{FALSE}.
#' @param perProtein A logical. \code{TRUE} ignores the site-level information
#'   for PTM and considers protein as a whole, \code{FALSE} otherwise. Default
#'   is \code{FALSE}.
#'
#' @return A list of two elements named \code{PTM} and \code{PROTEIN}. The
#'   \code{PTM} list has four elements: \code{protein} (a string vector of
#'   protein names), \code{site} (a string vector of PTM sites), \code{param}
#'   (a list of model parameter estimates for each site), and \code{df} (a
#'   numeric vector of degrees of freedom for each model). The \code{PROTEIN}
#'   list includes all as in \code{PTM}, except \code{site}.
#'
#' @examples
#' sim <- PTMsimulateExperiment(
#'     nGroup=2, nRep=2, nProtein=1, nSite=1, nFeature=5,
#'     logAbundance=list(
#'         PTM=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05),
#'         PROTEIN=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05)
#'     )
#' )
#' s <- PTMsummarize(sim)
#' estimateAbundance(s[["PTM"]])
#' estimateAbundance(s[["PROTEIN"]], perProtein=TRUE)
#'
#' @export
estimateAbundance <- function(df, fctBatch=FALSE, perProtein=FALSE) {
    cols_summarized <- c("protein", "group", "run", "log2inty")
    if (!perProtein) cols_summarized <- c(cols_summarized, "site")
    cols <- names(df)
    if (!all(cols_summarized %in% cols)) {
        stop(
            "Please include in the data frame all the following columns: ",
            paste0(sQuote(cols_summarized), collapse = ", ")
        )
    }
    if ("batch" %in% cols) {
        df <- df[, c(cols_summarized, "batch")]
        if (fctBatch && length(unique(df$batch)) == 1)
            stop("Cannot estimate batch effect with only one unique batch ID!")
    } else {
        df <- df[, cols_summarized]
        if (fctBatch)
            stop("To account for batch effect, add a 'batch' column in df")
    }

    cols_data <- c("group", "run", "log2inty")
    if (fctBatch) {
        nested <- nest(df, data = one_of(c(cols_data, "batch")))
        singles <- vapply(
            nested$data, function(x) length(unique(x$batch)) == 1, logical(1)
        )
        nested <- nested[!singles, ]
    } else {
        nested <- nest(df, data = one_of(cols_data))
    }
    # Fit linear models
    nested$lm_fit <- lapply(nested$data, function(x) {
        try(fitLinearModel(x, fctBatch), silent = TRUE)
    })
    errors <- vapply(
        nested$lm_fit, function(res) inherits(res, "try-error"), logical(1)
    )
    if (any(errors))
        warning("There were 1 or more errors while fitting models")
    nested <- nested[!errors, ]

    # Remove cases not eligible for hypothesis testing (SE is NA)
    nested$param <- Map(tidyEstimates, nested$lm_fit, nested$data)
    nested$df <- vapply(nested$lm_fit, df.residual, FUN.VALUE = double(1))
    nas <- vapply(
        nested$param, function(x) any(is.na(x$std.error)), logical(1)
    )
    nested <- nested[!nas, ]
    as.list(nested[, !(names(nested) %in% c("data", "lm_fit"))])
}


#' Fit linear model
#'
#' \code{fitLinearModel} fits and returns a linear model with \code{log2inty}
#' as response, and \code{group} and possibly \code{batch} as fixed effects.
#'
#' @param df A data frame with columns \code{log2inty}, \code{group}, and
#'   \code{batch} for one PTM site.
#' @param fctBatch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#'
#' @return An \code{lm} model object.
#'
#' @examples
#' x1 <- data.frame(
#'     batch=rep(c("1", "2"), each=4),
#'     group=rep(c("1", "2"), 4),
#'     log2inty=rep(c(10, 12), 4) + rnorm(8)
#' )
#' fitLinearModel(x1, fctBatch=TRUE)
#'
#' x2 <- data.frame(
#'     group=rep(c("1", "2"), 3),
#'     log2inty=rep(c(10, 12), 3) + rnorm(6)
#' )
#' fitLinearModel(x2)
#'
#' @export
fitLinearModel <- function(df, fctBatch=FALSE) {
    if (missing(df))
        stop(paste0("The input ", sQuote("df"), " is missing!"))
    if (!is.data.frame(df))
        stop(
            "Provide summarized log2-intensities as a data frame in ",
            sQuote("df")
        )
    if (fctBatch) {
        fit <- fixedGroupBatch(df)
    } else {
        fit <- fixedGroup(df)
    }
    fit
}


#' Linear model with group and batch effects
#'
#' \code{fixedGroupBatch} fits and returns a linear model with \code{log2inty}
#' as response, and \code{group} and \code{batch} as fixed effects.
#'
#' @param df A data frame with columns \code{log2inty}, \code{group}, and
#'   \code{batch} for one PTM site.
#'
#' @return An \code{lm} model object.
#'
#' @examples
#' x <- data.frame(
#'     batch=rep(c("1", "2"), each=4),
#'     group=rep(c("1", "2"), 4),
#'     log2inty=rep(c(10, 12), 4) + rnorm(8)
#' )
#' fixedGroupBatch(x)
#'
#' @export
fixedGroupBatch <- function(df) {
    if (length(unique(df$batch)) == 1)
        stop("Cannot estimate batch effect with a single batch!")
    if (length(unique(df$group)) == 1) {
        fit <- lm(log2inty ~ batch, data = df)
    } else {
        fit <- lm(log2inty ~ 0 + group + batch, data = df)
    }
    fit
}


#' Linear model with group effect
#'
#' \code{fixedGroup} fits and returns a linear model with \code{group} as a
#' fixed effect.
#'
#' @param df A data frame with columns \code{log2inty} and \code{group} for one
#'   PTM site.
#'
#' @return An \code{lm} model object.
#'
#' @examples
#' x <- data.frame(
#'     group=rep(c("1", "2"), 3),
#'     log2inty=rep(c(10, 12), 3) + rnorm(6)
#' )
#' fixedGroup(x)
#'
#' @export
fixedGroup <- function(df) {
    if (length(unique(df$group)) == 1) {
        fit <- lm(log2inty ~ 1, data = df)
    } else {
        fit <- lm(log2inty ~ 0 + group, data = df)
    }
    fit
}


#' Extract estimate of group effect
#'
#' \code{tidyEstimates} extracts the estimate of group effect from a fitted
#' linear model.
#'
#' @param fit An \code{lm} model object.
#' @param data A data frame used to derive the model object \code{fit}.
#'
#' @return A data frame restoring the estimated model parameters.
#'
#' @examples
#' x <- data.frame(
#'     group=rep(c("1", "2"), 3),
#'     log2inty=rep(c(10, 12), 3) + rnorm(6)
#' )
#' fit <- fitLinearModel(x)
#' tidyEstimates(fit, x)
#'
#' @export
tidyEstimates <- function(fit, data) {
    param <- broom::tidy(fit)
    batches <- grepl("batch", param$term)
    param <- param[!batches, ]
    if (length(unique(data$group)) == 1) {
        param$group <- data$group[1]
    } else {
        param$group <- gsub("group", "", param$term)
    }
    param[, !(names(param) %in% c("term", "statistic", "p.value"))]
}
