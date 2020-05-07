#' Estimate log2-abundances of PTM sites.
#'
#' \code{PTMestimateAbundance} takes as input the summarized log2-intensities
#' for each PTM site, performs statistical modeling for the abundance of the
#' site, and returns the estimates of model parameters for all sites in all
#' experimental conditions.
#'
#' @param df A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{log2inty}, and possibly, \code{batch}.
#' @param fac_batch A logical. \code{TRUE} considers a fixed batch effect,
#'   \code{FALSE} otherwise. Default is \code{FALSE}.
#' @return A list.
#'
#' @export
#' @examples
#' PTMestimateAbundance(df, fac_batch = FALSE)
PTMestimateAbundance <- function(df, fac_batch = FALSE) {
    if (missing(df))
        stop("Input data frame is missing!")
    if (!is.data.frame(df))
        stop(paste0("Provide protein summaries as a data frame in ", sQuote("df")))
    cols_summarized <- c("protein", "site", "group", "run", "log2inty")
    cols <- names(df)
    if (!all(cols_summarized %in% cols)) {
        stop("Please include in the data frame all the following columns: ",
             paste0(sQuote(cols_summarized), collapse = ", "))
    }
    if (fac_batch && !("batch" %in% cols))
        stop("To account for batch effect, add a 'batch' column in df")

    if (fac_batch) {
        if (length(unique(df$batch)) == 1)
            stop("There is only one unique batch identifier!")

        # One model for all batches (need data with >1 batches)
        nested <- nest(df, data = -one_of("protein", "site"))
        singles <- sapply(nested$data, function(dat) length(unique(dat$batch)) == 1)
        nested <- nested[!singles, ]
    } else {
        # One model per site (and potentially batch)
        if ("batch" %in% names(df)) {
            nested <- nest(df, data = -one_of("protein", "site", "batch"))
        } else {
            nested <- nest(df, data = -one_of("protein", "site"))
        }
    }
    # Fit linear models
    nested$lm_fit <- lapply(nested$data, function(x) {
        try(fitLinearModel(x, fac_batch), silent = TRUE)
    })
    errors <- sapply(nested$lm_fit, function(res) inherits(res, "try-error"))
    if (any(errors)) {
        warning("There were 1 or more errors while fitting models")
    }
    nested <- nested[!errors, ]

    # Remove cases not eligible for hypothesis testing (SE is NA)
    nested$param <- purrr::map2(nested$lm_fit, nested$data, tidyEstimates)
    nas <- sapply(nested$param, function(res) any(is.na(res$std.error)))
    nested <- nested[!nas, ]
    as.list(nested[, !(names(nested) %in% c("data", "lm_fit"))])
}


#' Fit linear model.
#'
#' \code{fitLinearModel} fits and returns a linear model with \code{log2inty}
#' as response, and \code{group} and possibly \code{batch} as fixed effects.
#'
#' @param df A data frame with columns \code{log2inty}, \code{group}, and
#'   \code{batch} for one PTM site.
#' @param fac_batch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return An \code{lm} model object.
#'
#' @export
#' @examples
#' x1 <- data.frame(
#'   batch = rep(c("1", "2"), each = 4),
#'   group = rep(c("1", "2"), 4),
#'   log2inty = rep(c(10, 12), 4) + rnorm(8)
#' )
#' fitLinearModel(x1, fac_batch = TRUE)
#'
#' x2 <- data.frame(
#'   group = rep(c("1", "2"), 3),
#'   log2inty = rep(c(10, 12), 3) + rnorm(6)
#' )
#' fitLinearModel(x2)
fitLinearModel <- function(df, fac_batch = FALSE) {
    if (missing(df))
        stop(paste0("The input ", sQuote("df"), " is missing!"))
    if (!is.data.frame(df))
        stop("Provide summarized log2-intensities as a data frame in ", sQuote("df"))

    if (fac_batch) {
        fit <- fixedGroupBatch(df)
    } else {
        fit <- fixedGroup(df)
    }
    fit
}


#' Linear model with group and batch effects.
#'
#' \code{fixedGroupBatch} fits and returns a linear model with \code{log2inty}
#' as response, and \code{group} and \code{batch} as fixed effects.
#'
#' @param df A data frame with columns \code{log2inty}, \code{group}, and
#'   \code{batch} for one PTM site.
#' @return An \code{lm} model object.
#'
#' @export
#' @examples
#' x <- data.frame(
#'   batch = rep(c("1", "2"), each = 4),
#'   group = rep(c("1", "2"), 4),
#'   log2inty = rep(c(10, 12), 4) + rnorm(8)
#' )
#' fixedGroupBatch(x)
fixedGroupBatch <- function(df) {
    if (length(unique(df$batch)) == 1)
        stop("Cannot estimate batch effect with a single batch!")
    if (length(unique(df$group)) == 1) {
        fit <- stats::lm(log2inty ~ batch, data = df)
    } else {
        fit <- stats::lm(log2inty ~ 0 + group + batch, data = df)
    }
    fit
}


#' Linear model with group effect.
#'
#' \code{fixedGroup} fits and returns a linear model with \code{group} as a
#' fixed effect.
#'
#' @param df A data frame with columns \code{log2inty} and \code{group} for one
#'   PTM site.
#' @return An \code{lm} model object.
#'
#' @export
#' @examples
#' x <- data.frame(
#'   group = rep(c("1", "2"), 3),
#'   log2inty = rep(c(10, 12), 3) + rnorm(6)
#' )
#' fixedGroup(x)
fixedGroup <- function(df) {
    if (length(unique(df$group)) == 1) {
        fit <- stats::lm(log2inty ~ 1, data = df)
    } else {
        fit <- stats::lm(log2inty ~ 0 + group, data = df)
    }
    fit
}


#' Extract estimate of group effect.
#'
#' \code{tidyEstimates} extracts the estimate of group effect from a fitted
#' linear model.
#'
#' @param fit An \code{lm} model object.
#' @param data A data frame used to derive the model object \code{fit}.
#' @return A data frame restoring the estimated model parameters.
#'
#' @export
#' @examples
#' x <- data.frame(
#'   group = rep(c("1", "2"), 3),
#'   log2inty = rep(c(10, 12), 3) + rnorm(6)
#' )
#' fit <- fitLinearModel(x)
#' tidyEstimates(fit, x)
tidyEstimates <- function(fit, data) {
    param <- broom::tidy(fit)
    batches <- grepl("batch", param$term)
    param <- param[!batches, ]
    # param <- broom::tidy(fit) %>%
    #     dplyr::filter(!stringr::str_detect(term, "batch"))
    if (length(unique(data$group)) == 1) {
        param$group <- data$group[1]
        # param <- param %>%
        #     dplyr::mutate(group = unique(data$group)) %>%
        #     dplyr::select(-term, -statistic, -p.value)
    } else {
        param$group <- gsub("group", "", param$term)
        # param <- param %>%
        #     dplyr::mutate(group = stringr::str_replace(term, "group", "")) %>%
        #     dplyr::select(-term, -statistic, -p.value)
    }
    param$df <- (sigma(fit) ^ 2 / param$std.error ^ 2) - 1
    param[, !(names(param) %in% c("term", "statistic", "p.value"))]
    # param %>%
    #     dplyr::mutate(df = (sigma(fit) ^ 2 / std.error ^ 2) - 1)
}
