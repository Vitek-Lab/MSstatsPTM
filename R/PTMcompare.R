#' Compare mean abundances for all PTM sites across conditions.
#'
#' \code{PTMcompareMeans} performs significance analysis for detection of
#'   changes in PTM mean abundances between conditions.
#'
#' @param data A data frame.
#' @param controls Control group as reference in all the comparisons.
#' @param cases Case groups in all the comparisons.
#' @return A data frame.
#' @export
#'
#' @examples
#' PTMcompareMeans(df_mod, controls, cases)
PTMcompareMeans <- function(data, controls, cases) {
    nested <- dplyr::tibble(
        protein = data[["protein"]],
        site = data[["site"]],
        param = data[["param"]]
    )
    unnested <- tidyr::unnest(nested, one_of("param"))

    groups <- unique(c(controls, cases))
    in_grp <- unnested$group %in% groups
    nna <- !(is.na(unnested$estimate) | is.na(unnested$std.error))
    unnested <- unnested[in_grp & nna, ]
    # unnested <- nested %>%
    #     tidyr::unnest(one_of("param")) %>%
    #     dplyr::filter(group %in% groups) %>%
    #     dplyr::filter(!is.na(estimate), !is.na(std.error))

    cols_res <- c("protein", "site", "Label", "log2FC", "SE", "Tvalue", "DF",
                  "pvalue", "adj.pvalue")
    res <- vector("list", length(cases))
    for (i in seq_along(cases)) {
        case <- cases[i]
        ctrl <- controls[i]
        twogroup <- unnested[unnested$group %in% c(case, ctrl), ]
        twogroup$group <- ifelse(twogroup$group == case, "case", "ctrl")
        # twogroup <- unnested %>%
        #     dplyr::filter(group %in% c(case, ctrl)) %>%
        #     dplyr::mutate(group = ifelse(group == case, "case", "ctrl"))

        # Estimated in both groups
        cnt <- dplyr::count(twogroup, protein, site)
        cnt2 <- cnt[cnt$n == 2, ]
        full <- dplyr::semi_join(twogroup, cnt2)
        single <- dplyr::anti_join(twogroup, cnt2)
        # full <- twogroup %>%
        #     dplyr::group_by(protein, site) %>%
        #     dplyr::filter(dplyr::n() == 2) %>%
        #     dplyr::ungroup()
        # # Missing in one group
        # single <- twogroup %>%
        #     anti_join(full, by = c("protein", "site"))

        # Hypothesis testing
        full$s2 <- full$std.error ^ 2
        full_wd <- tidyr::pivot_wider(
            full[, names(full) != "std.error"],
            names_from = group,
            values_from = c(estimate, s2, df)
        )
        s2_sum <- full_wd$s2_case + full_wd$s2_ctrl
        full_wd$Label <- paste(case, ctrl, sep = " vs ")
        full_wd$log2FC <- full_wd$estimate_case - full_wd$estimate_ctrl
        full_wd$SE <- sqrt(s2_sum)
        full_wd$Tvalue <- full_wd$log2FC / full_wd$SE
        numer <- s2_sum ^ 2
        denom <- full_wd$s2_case ^ 2 / full_wd$df_case +
            full_wd$s2_ctrl ^ 2 / full_wd$df_ctrl
        full_wd$DF <- numer / denom
        full_wd$pvalue <- 2 * pt(abs(full_wd$Tvalue), df = full_wd$DF,
                                 lower.tail = FALSE)
        full_wd$adj.pvalue <- stats::p.adjust(full_wd$pvalue, method = "BH")

        single$Label <- paste(case, ctrl, sep = " vs ")
        single$log2FC <- ifelse(single$group == case, Inf, -Inf)
        single$SE <- NA
        single$Tvalue <- NA
        single$DF <- NA
        single$pvalue <- NA
        single$adj.pvalue <- 0
        # full <- full %>%
        #     dplyr::mutate(s2 = std.error ^ 2) %>%
        #     dplyr::select(-std.error) %>%
        #     tidyr::pivot_wider(names_from = group, values_from = c(estimate, s2, df)) %>%
        #     dplyr::mutate(
        #         s2_sum = s2_case + s2_ctrl,
        #         Label = paste(case, ctrl, sep = " vs "),
        #         log2FC = estimate_case - estimate_ctrl,
        #         SE = sqrt(s2_sum),
        #         Tvalue = log2FC / SE,
        #         DF = s2_sum ^ 2 / (s2_case ^ 2 / df_case + s2_ctrl ^ 2 / df_ctrl),
        #         pvalue = 2 * pt(abs(Tvalue), df = DF, lower.tail = FALSE),
        #         adj.pvalue = stats::p.adjust(pvalue, method = "BH")
        #     )
        # single <- single %>%
        #     dplyr::mutate(
        #         Label = paste(case, ctrl, sep = " vs "),
        #         log2FC = ifelse(group == case, Inf, -Inf),
        #         SE = NA, Tvalue = NA, DF = NA, pvalue = NA, adj.pvalue = 0
        #     )
        res[[i]] <- dplyr::bind_rows(full[, cols_res], single[, cols_res])
    }
    dplyr::bind_rows(res)
}




#' Compare mean abundances for all PTM sites across conditions.
#'
#' \code{PTMcompareMeans} performs significance analysis for detection of
#'   changes in PTM mean abundances between conditions.
#'
#' @param data A list that contain estimates for PTM site abundance
#'   (\code{data$PTM}), and possibly, protein abundance (\code{data$Protein}).
#' @param controls A string vector of control groups as reference.
#' @param cases A string vector of case groups.
#' @param adjProtein A logical. \code{TRUE} performs protein-level adjustment,
#'   \code{FALSE} otherwise. Default is \code{FALSE}.
#' @return A data frame.
#' @export
#'
#' @examples
#' PTMcompareMeans(df_mod, controls, cases)
PTMcompareMeans <- function(data, controls, cases, adjProtein = FALSE) {
    if (is.null(data[["PTM"]]))
        stop("PTM estimates are not provided")
    if (is.null(data[["Protein"]]) && adjProtein)
        stop("To adjust for protein level, protein estimates are required")

    res <- extractMeanDiff(data[["PTM"]], controls, cases)
    if (adjProteinLevel) {
        res_ref <- extractMeanDiff(data[["Protein"]], controls, cases)
        res <- adjustProteinEstimates(res, res_ref)
    }

    res
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
    if (length(unique(data$group)) == 1) {
        param$group <- data$group[1]
    } else {
        param$group <- gsub("group", "", param$term)
    }
    res <- list(param = param[, c("group", "estimate", "std.error")], df = df.residual(fit))
    res
}

extractMeanDiff <- function(data, controls, cases) {
    # Params per site in a row
    params <- dplyr::tibble(
        protein = data[["protein"]],
        site = data[["site"]],
        param = data[["param"]],
        df = data[["df"]]
    )

    res <- vector("list", length(controls))
    for (i in seq_along(controls)) {
        onectrx <- vector("list", nrow(params))
        for (j in seq_along(params$param)) {
            onectrx[[j]] <- onetest(params$param[[j]], params$df[j],
                                    controls[i], cases[i])
        }
        res[[i]] <- dplyr::bind_rows(onectrx)
    }
    dplyr::bind_rows(res)
}

onetest <- function(param, df, ctrl, case) {
    if (!any(c(ctrl, case) %in% param$group)) {
        return(NULL)
    } else if (sum(c(ctrl, case) %in% param$group) == 1) {
        res <- data.frame(Label = paste(case, ctrl, sep = " vs "),
                          log2FC = ifelse(case %in% param$group, Inf, -Inf),
                          SE = NA, Tvalue = NA, DF = NA, pvalue = NA)
    } else {
        param_case <- param[param$group == case, ]
        param_ctrl <- param[param$group == ctrl, ]
        log2fc <- param_case$estimate - param_ctrl$estimate
        stderr <- sqrt(param_case$std.error ^ 2 + param_ctrl$std.error ^ 2)
        tval <- log2fc / stderr
        pval <- 2 * pt(abs(tval), df, lower.tail = FALSE)
        res <- data.frame(Label = paste(case, ctrl, sep = " vs "),
                          log2FC = log2fc, SE = stderr, Tvalue = tval, DF = df,
                          pvalue = pval)
    }
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


    cols_res <- c("protein", "site", "Label", "log2FC", "SE", "Tvalue", "DF",
                  "pvalue", "adj.pvalue")
    res <- vector("list", length(cases))
    for (i in seq_along(cases)) {
        case <- cases[i]
        ctrl <- controls[i]
        twogroup <- unnested[unnested$group %in% c(case, ctrl), ]
        twogroup$group <- ifelse(twogroup$group == case, "case", "ctrl")

        # Estimated in both groups
        cnt <- dplyr::count(twogroup, protein, site)
        cnt2 <- cnt[cnt$n == 2, ]
        full <- dplyr::semi_join(twogroup, cnt2)
        single <- dplyr::anti_join(twogroup, cnt2)

        # Hypothesis testing
        full$s2 <- full$std.error ^ 2
        full_wd <- tidyr::pivot_wider(
            full[, names(full) != "std.error"],
            names_from = group,
            values_from = c(estimate, s2, df)
        )
        s2_sum <- full_wd$s2_case + full_wd$s2_ctrl
        full_wd$Label <- paste(case, ctrl, sep = " vs ")
        full_wd$log2FC <- full_wd$estimate_case - full_wd$estimate_ctrl
        full_wd$SE <- sqrt(s2_sum)
        full_wd$Tvalue <- full_wd$log2FC / full_wd$SE
        numer <- s2_sum ^ 2
        denom <- full_wd$s2_case ^ 2 / full_wd$df_case +
            full_wd$s2_ctrl ^ 2 / full_wd$df_ctrl
        full_wd$DF <- numer / denom
        full_wd$pvalue <- 2 * pt(abs(full_wd$Tvalue), df = full_wd$DF,
                                 lower.tail = FALSE)
        full_wd$adj.pvalue <- stats::p.adjust(full_wd$pvalue, method = "BH")

        single$Label <- paste(case, ctrl, sep = " vs ")
        single$log2FC <- ifelse(single$group == case, Inf, -Inf)
        single$SE <- NA
        single$Tvalue <- NA
        single$DF <- NA
        single$pvalue <- NA
        single$adj.pvalue <- 0
        res[[i]] <- dplyr::bind_rows(full[, cols_res], single[, cols_res])
    }
    dplyr::bind_rows(res)
}
