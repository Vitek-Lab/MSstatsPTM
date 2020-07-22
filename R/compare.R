#' Compare mean abundances for all PTM sites across conditions
#'
#' \code{PTMcompareMeans} performs significance analysis for detection of
#' changes in PTM mean abundances between conditions.
#'
#' @param data A list of two elements named \code{PTM} and \code{Protein}. The
#'   \code{PTM} list has four elements: \code{protein} (a string vector of
#'   protein names), \code{site} (a string vector of PTM sites), \code{param}
#'   (a list of model parameter estimates for each site), and \code{df} (a
#'   numeric vector of degrees of freedom for each model). The \code{Protein}
#'   list includes all as in \code{PTM}, except the element \code{site}.
#' @param controls A string vector of control groups in the comparisons.
#' @param cases A string vector of case groups.
#' @param adjProtein A logical. \code{TRUE} performs protein-level adjustment,
#'   \code{FALSE} otherwise. Default is \code{FALSE}.
#'
#' @return A data frame.
#'
#' @examples
#' \dontrun{
#' PTMcompareMeans(df_mod, controls, cases)
#' }
#'
#' @export
PTMcompareMeans <- function(data, controls, cases, adjProtein = FALSE) {
    # Check PTM estimates
    if (is.null(data[["PTM"]]))
        stop("PTM estimates are missing!")
    elm_site <- c("protein", "site", "param", "df")
    if (!all(elm_site %in% names(data[["PTM"]]))) {
        stop("Please include in the PTM list all the following elements: ",
             paste0(sQuote(elm_site), collapse = ", "))
    }
    if (adjProtein) {
        # Check Protein estimates
        if (is.null(data[["Protein"]]))
            stop("To adjust for protein level, protein estimates are required")
        elm_prot <- c("protein", "param", "df")
        if (!all(elm_prot %in% names(data[["Protein"]]))) {
            stop("Please include in the Protein list all the following elements: ",
                 paste0(sQuote(elm_site), collapse = ", "))
        }
    }
    if (!(is.character(controls) && is.character(cases)))
        stop("Provide the control and case groups as character vectors")
    if (length(controls) != length(cases))
        stop(paste0("The lengths of ", sQuote("controls"), " and ",
                    sQuote("cases"), " should be identical"))

    res <- extractMeanDiff(data[["PTM"]], controls, cases, per_protein = FALSE)
    if (adjProtein) {
        res_prot <- extractMeanDiff(data[["Protein"]], controls, cases,
                                    per_protein = TRUE)
        res <- .adjustProteinLevel(res, res_prot)
    }
    res
}


#' Compare mean abundances for all PTM sites (or proteins) across conditions
#'
#' \code{extractMeanDiff} performs significance analysis for detection of
#' changes in PTM mean abundances between conditions.
#'
#' @param data A list of abundance estimates with the following elements:
#'   \code{protein}, \code{site}, \code{param}, and \code{df}. \code{site} may
#'   be excluded when performing per-protein analysis (\code{per_protein} is
#'   \code{TRUE}).
#' @param controls A string vector of control groups in the comparisons.
#' @param cases A string vector of case groups.
#' @param per_protein A logical. \code{TRUE} ignores the site-level information
#'   for PTM and considers protein as a whole, \code{FALSE} performs site-level
#'   analysis. Default is \code{FALSE}.
#'
#' @return A data frame.
#'
#' @export
extractMeanDiff <- function(data, controls, cases, per_protein = FALSE) {
    if (!per_protein && is.null(data[["site"]]))
        stop("Site-level analysis requires the information of PTM site")

    w_batch <- "batch" %in% names(data)
    if (per_protein) {
        # One row per protein
        params <- tibble(
            protein = data[["protein"]],
            param = data[["param"]],
            df = data[["df"]]
        )
        if (w_batch) {
            nb <- length(unique(data[["batch"]]))
            params$batch <- data[["batch"]]
            cnt <- count(params, .data$protein)
            cnt <- cnt[cnt$n == nb, ]
            params <- semi_join(params, cnt)
        }
    } else {
        # One row per site
        params <- tibble(
            protein = data[["protein"]],
            site = data[["site"]],
            param = data[["param"]],
            df = data[["df"]]
        )
        if (w_batch) {
            nb <- length(unique(data[["batch"]]))
            params$batch <- data[["batch"]]
            cnt <- count(params, .data$protein, .data$site)
            cnt <- cnt[cnt$n == nb, ]
            params <- semi_join(params, cnt)
        }
    }

    res <- vector("list", length(controls))
    for (i in seq_along(controls)) {
        # Test for one contrast
        ctrl <- controls[i]
        case <- cases[i]
        tests <- Map(.onetest, params$param, params$df, ctrl, case)
        nores <- sapply(tests, is.null)
        onectrx <- bind_rows(tests)
        onectrx$Protein <- params$protein[!nores]
        if (!per_protein) {
            onectrx$Site <- params$site[!nores]
        }
        if (w_batch) {
            if (per_protein) {
                cnt <- count(onectrx, .data$Protein)
            } else {
                cnt <- count(onectrx, .data$Protein, .data$Site)
            }
            cnt <- cnt[cnt$n == nb, ]
            onectrx <- semi_join(onectrx, cnt)
            onectrx <- onectrx[, !(names(onectrx) %in% c("Tvalue", "pvalue"))]
            nested <- nest(onectrx, data = one_of(c("log2FC", "SE", "DF")))
            nested$ave <- lapply(nested$data, .aggregateFC)
            nested <- nested[, !(names(nested) %in% c("data"))]
            onectrx <- unnest(nested, one_of("ave"))
        }
        cols <- names(onectrx)
        is_first <- cols %in% c("Protein", "Site")
        onectrx <- onectrx[, c(cols[is_first], cols[!is_first])]
        res[[i]] <- onectrx
    }
    bind_rows(res)
}


#' Test for one contrast in one site (or protein)
#' @keywords internal
.onetest <- function(param, df, ctrl, case) {
    if (!any(c(ctrl, case) %in% param$group)) {
        return(NULL)
    } else if (sum(c(ctrl, case) %in% param$group) == 1) {
        res <- tibble(Label = paste(case, ctrl, sep = " vs "),
                      log2FC = ifelse(case %in% param$group, Inf, -Inf),
                      SE = NA, Tvalue = NA, DF = NA, pvalue = NA)
    } else {
        param_case <- param[param$group == case, ]
        param_ctrl <- param[param$group == ctrl, ]
        log2fc <- param_case$estimate - param_ctrl$estimate
        stderr <- sqrt(param_case$std.error ^ 2 + param_ctrl$std.error ^ 2)
        tval <- log2fc / stderr
        pval <- 2 * stats::pt(abs(tval), df, lower.tail = FALSE)
        res <- tibble(Label = paste(case, ctrl, sep = " vs "), log2FC = log2fc,
                      SE = stderr, Tvalue = tval, DF = df, pvalue = pval)
    }
    res
}


#' Aggregate FC estimates from multiple sources by averaging over them
#' @keywords internal
.aggregateFC <- function(data) {
    log2fc <- mean(data$log2FC)
    s2 <- data$SE ^ 2
    stderr <- sqrt(sum(s2)) / nrow(data)
    numer <- sum(s2) ^ 2
    denom <- sum(s2 ^ 2 / data$DF)
    df <- numer / denom
    tval <- log2fc / stderr
    pval <- 2 * stats::pt(abs(tval), df, lower.tail = FALSE)

    tibble(log2FC = log2fc, SE = stderr, Tvalue = tval, DF = df, pvalue = pval)
}


#' Protein-level adjustment
#' @keywords internal
.adjustProteinLevel <- function(diffSite, diffProtein) {
    diffRef <- diffProtein[, c("Protein", "Label", "log2FC", "SE", "DF")]
    names(diffRef)[names(diffRef) == "log2FC"] <- "log2FC_ref"
    names(diffRef)[names(diffRef) == "SE"] <- "SE_ref"
    names(diffRef)[names(diffRef) == "DF"] <- "DF_ref"
    joined <- inner_join(diffSite, diffRef)

    missing_ctrl <- joined[joined$log2FC == Inf, ]  # PTM missing in control
    missing_case <- joined[joined$log2FC == -Inf, ] # PTM missing in case
    res_mctrl <- tibble(Protein = missing_ctrl$Protein,
                        Site = missing_ctrl$Site, Label = missing_ctrl$Label,
                        log2FC = Inf, SE = NA, Tvalue = NA, DF = NA, pvalue = NA)
    res_mcase <- tibble(Protein = missing_case$Protein,
                        Site = missing_case$Site, Label = missing_case$Label,
                        log2FC = -Inf, SE = NA, Tvalue = NA, DF = NA, pvalue = NA)

    idx_full <- abs(joined$log2FC) != Inf & abs(joined$log2FC_ref) != Inf
    full <- joined[idx_full, ]
    log2fc <- full$log2FC - full$log2FC_ref
    s2 <- full$SE ^ 2
    s2_ref <- full$SE_ref ^ 2
    stderr <- sqrt(s2 + s2_ref)
    numer <- (s2 + s2_ref) ^ 2
    denom <- (s2 ^ 2 / full$DF + s2_ref ^ 2 / full$DF_ref)
    df <- numer / denom
    tval <- log2fc / stderr
    pval <- 2 * stats::pt(abs(tval), df, lower.tail = FALSE)
    res_full <- tibble(Protein = full$Protein,
                       Site = full$Site,
                       Label = full$Label,
                       log2FC = log2fc,
                       SE = stderr,
                       Tvalue = tval,
                       DF = df,
                       pvalue = pval)

    bind_rows(res_full, res_mctrl, res_mcase)
}
