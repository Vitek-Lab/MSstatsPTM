#' Estimate log-abundances of PTM sites.
#'
#' \code{PTMestimateAbundance} fits a linear model with summarized
#'   log-intensities for each PTM site, in consideration of the experimental
#'   design and the underlying protein abundance, and returns the estimates of
#'   model parameters.
#'
#' @param df A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{log2inty}, and possibly, \code{batch}.
#' @param fac_batch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return A data frame.
#' @export
#' @rdname estimation
#'
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

    res <- modelAbundance(df, fac_batch)
    if (("batch" %in% cols) && !fac_batch) {
        res <- aggregateBatchEstimates(res)
    }
    res
}

#' Estimate log-abundance of a PTM site per group.
#'
#' \code{modelAbundance} fits and returns a whole-plot model for each site in
#'   consideration of the experimental design.
#'
#' @importFrom tidyr nest
#' @importFrom tidyselect one_of
#' @param df A data frame.
#' @param fac_batch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return A list.
#' @export
#' @rdname estimation
#'
#' @examples
#' modelAbundance(df, fac_batch = FALSE)
modelAbundance <- function(df, fac_batch = FALSE) {
    if (fac_batch) {
        if (!("batch" %in% names(df)))
            stop("There is no information about batch!")
        if (length(unique(df$batch)) == 1)
            stop("There is only one unique batch identifier!")

        # One model for all batches (need data with >1 batches)
        nested <- tidyr::nest(df, data = -one_of("protein", "site"))
        singles <- sapply(nested$data, function(dat) length(unique(dat$batch)) == 1)
        nested <- nested[!singles, ]
        # nested <- df %>%
        #     tidyr::nest(data = -one_of("protein", "site")) %>%
        #     dplyr::filter(purrr::map_int(data, ~ dplyr::n_distinct(.$batch)) > 1)
    } else {
        # One model per site (and potentially batch)
        if ("batch" %in% names(df)) {
            nested <- tidyr::nest(df, data = -one_of("protein", "site", "batch"))
        } else {
            nested <- tidyr::nest(df, data = -one_of("protein", "site"))
        }
    }
    # Fit linear models
    nested$lm_fit <- lapply(nested$data, function(x) {
        try(fitLinearModels(x, fac_batch), silent = TRUE)
    })
    errors <- sapply(nested$lm_fit, function(res) inherits(res, "try-error"))
    if (any(errors)) {
        warning("There were 1 or more errors while fitting models")
    }
    nested <- nested[!errors, ]
    # nested <- nested %>%
    #     dplyr::mutate(lm_fit = purrr::map(
    #         data, purrr::possibly(fitLinearModels, otherwise = NA_real_), fac_batch
    #     ))
    # nested <- nested %>%
    #     dplyr::filter(!is.na(lm_fit))

    # Remove cases not eligible for hypothesis testing (SE is NA)
    nested$param <- purrr::map2(nested$lm_fit, nested$data, tidyEstimates)
    # nested$param <- mapply(tidyEstimates, nested$lm_fit, nested$data)
    nas <- sapply(nested$param, function(res) any(is.na(res$std.error)))
    nested <- nested[!nas, ]
    # nested <- nested %>%
    #     dplyr::mutate(param = purrr::map2(lm_fit, data, tidyEstimates)) %>%
    #     dplyr::filter(!purrr::map_lgl(param, ~any(is.na(.$std.error))))

    # nested %>%
    #     select(-one_of("data", "lm_fit")) %>%
    #     as.list()
    as.list(nested[, !(names(nested) %in% c("data", "lm_fit"))])
}

#' Extract estimate of group effect.
#'
#' \code{tidyEstimates} extracts the estimate of group effect.
#'
#' @param fit An \code{lm} model object.
#' @param data A data frame.
#' @return A data frame restoring the estimated model parameters.
#' @export
#'
#' @examples
#' tidyEstimates(fit, data)
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

#' Aggregate parameter estimates across batches.
aggregateBatchEstimates <- function(data) {
    if (!("batch" %in% names(data)))
        stop("There is no information about batch!")

    n_batch <- length(unique(data$batch))
    if (n_batch == 1) {
        res <- list(protein = data[["protein"]],
                    site = data[["site"]],
                    param = data[["param"]])
        return(res)
    }

    # Params from all batches per site in a row
    df <- dplyr::tibble(
        protein = data[["protein"]],
        site = data[["site"]],
        batch = data[["batch"]],
        param = data[["param"]]
    )
    cnt <- dplyr::count(df, protein, site)
    df <- dplyr::semi_join(df, cnt[cnt$n == n_batch, ])
    nested <- tidyr::nest(df, params = one_of("batch", "param"))
    # nested <- dplyr::tibble(
    #     protein = data[["protein"]],
    #     site = data[["site"]],
    #     batch = data[["batch"]],
    #     param = data[["param"]]
    # ) %>%
    #     dplyr::group_by(protein, site) %>%
    #     filter(dplyr::n() == n_batch) %>%
    #     dplyr::ungroup() %>%
    #     tidyr::nest(params = one_of("batch", "param"))

    aggre_params <- function(params, nb) {
        unnested <- tidyr::unnest(params, one_of("param"))
        grps <- unique(unnested$group)
        valid <- rep(TRUE, length(grps))
        hat <- se <- df <- vector("double", length(grps))
        for (i in seq_along(grps)) {
            sub <- unnested[unnested$group == grps[i], ]
            if (nrow(sub) == nb) {
                s2 <- sub$std.error ^ 2
                hat[i] <- mean(sub$estimate)
                se[i] <- sqrt(sum(s2)) / nb
                df[i] <- sum(s2) ^ 2 / sum(s2 ^ 2 / sub$df)
            } else {
                valid[i] <- FALSE
            }
        }
        data.frame(group = grps[valid], estimate = hat[valid],
                   std.error = se[valid], df = df[valid])
        # params %>%
        #     tidyr::unnest(one_of("param")) %>%
        #     dplyr::group_by(group) %>%
        #     dplyr::filter(dplyr::n() == nb) %>%
        #     dplyr::mutate(s2 = std.error ^ 2) %>%
        #     dplyr::summarise(
        #         estimate = mean(estimate),
        #         std.error = sqrt(sum(s2)) / nb,
        #         df = sum(s2) ^ 2 / sum(s2 ^ 2 / df)
        #     )
    }

    nested$param <- lapply(nested$params, aggre_params, n_batch)
    as.list(nested[, c("protein", "site", "param")])
    # nested %>%
    #     select(one_of("protein", "site", "param")) %>%
    #     as.list()
}

#' Protein-level adjustment.
adjustProteinEstimates <- function(data, refdata) {
    nested <- dplyr::tibble(
        protein = data[["protein"]],
        site = data[["site"]],
        param = data[["param"]]
    )
    nested_ref <- dplyr::tibble(
        protein = refdata[["protein"]],
        param_ref = refdata[["param"]]
    )
    joined <- dplyr::inner_join(nested, nested_ref)

    adj_params <- function(param, param_ref) {
        names(param_ref)[names(param_ref) == "estimate"] <- "estimate_ref"
        names(param_ref)[names(param_ref) == "std.error"] <- "std.error_ref"
        names(param_ref)[names(param_ref) == "df"] <- "df_ref"
        # param_ref <- param_ref %>%
        #     dplyr::rename(estimate_ref = estimate,
        #                   std.error_ref = std.error,
        #                   df_ref = df)

        res <- dplyr::inner_join(param, param_ref)
        s2 <- res$std.error ^ 2
        s2_ref <- res$std.error_ref ^ 2
        res$estimate <- res$estimate - res$estimate_ref
        res$std.error <- sqrt(s2 + s2_ref) / 2
        numer <- (s2 + s2_ref) ^ 2
        denom <- (s2 ^ 2 / res$df + s2_ref ^ 2 / res$df_ref)
        res$df <- numer / denom
        res[, c("group", "estimate", "std.error", "df")]
        # dplyr::inner_join(param, param_ref) %>%
        #     dplyr::mutate(s2 = std.error ^ 2, s2_ref = std.error_ref ^ 2) %>%
        #     dplyr::transmute(
        #         group = group,
        #         estimate = estimate - estimate_ref,
        #         std.error = sqrt(s2 + s2_ref) / 2,
        #         df = (s2 + s2_ref) ^ 2 / (s2 ^ 2 / df + s2_ref ^ 2 / df_ref)
        #     )
    }

    joined$param <- mapply(adj_params, joined$param, joined$param_ref)
    as.list(joined[, c("protein", "site", "param")])
    # joined$param <- purrr::map2(joined$param, joined$param_ref, adj_params)
    # joined %>%
    #     select(one_of("protein", "site", "param")) %>%
    #     as.list()
}

#' Fit linear models.
#'
#' \code{fitLinearModels} fits linear model based on summarized data.
#'
#' @param df_prot A data frame.
#' @param fac_batch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return An object of class \code{lm}.
#' @export
#' @rdname estimation
#'
fitLinearModels <- function(df, fac_batch = FALSE) {
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

#' Fit a linear model with group effect and batch effect.
#'
#' \code{fixedGroupBatch} fits and returns a linear model with \code{group} effect.
#'
#' @param df A data frame containing columns \code{log2inty},
#'   \code{group} for one modification site.
#' @return An object of class \code{lm}.
#' @export
#'
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

#' Fit a linear model with group effect.
#'
#' \code{fixedGroup} fits and returns a linear model with \code{group} effect.
#'
#' @param df A data frame containing columns \code{log2inty},
#'   \code{group} for one modification site.
#' @return An object of class \code{lm}.
#' @export
#'
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

#' Combine estimated parameters.
#'
#' \code{combineParameters} combines the estimated model parameters from nested data
#'   frame. This function will soon be retired.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select rename left_join
#' @importFrom tidyr unnest
#' @param data A list.
#' @return A data frame restoring the estimated model parameters.
#' @export
#' @rdname estimation
#'
#' @examples
#' combineParameters(lparams)
combineParameters <- function(data) {
    if ("batch" %in% names(data)) {
        # per-batch model
        nested_param <- dplyr::tibble(
            protein = data[["protein"]],
            site = data[["site"]],
            batch = data[["batch"]],
            param = data[["param"]],
            df_res = data[["df.residual"]]
        )
    } else {
        # all-batch model
        nested_param <- dplyr::tibble(
            protein = data[["protein"]],
            site = data[["site"]],
            param = data[["param"]],
            df_res = data[["df.residual"]]
        )
    }
    param_mod <- nested_param %>%
        dplyr::filter(site != "None") %>%
        tidyr::unnest(one_of("param"))
    param_unmod <- nested_param %>%
        dplyr::filter(site == "None") %>%
        dplyr::select(-one_of("site")) %>%
        tidyr::unnest(one_of("param")) %>%
        dplyr::rename(df_unmod = df_res, est_unmod = estimate, se_unmod = std.error)

    return(dplyr::left_join(param_mod, param_unmod))
}
