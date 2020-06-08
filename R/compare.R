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
#' @return A data frame.
#'
#' @export
#' @examples
#' \dontrun{
#' PTMcompareMeans(df_mod, controls, cases)
#' }
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
#' \code{PTMcompareMeans} performs significance analysis for detection of
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
#' @return A data frame.
#'
#' @export
extractMeanDiff <- function(data, controls, cases, per_protein = FALSE) {
    if (!per_protein && is.null(data[["site"]]))
        stop("Site-level analysis requires the information of PTM site")

    if (per_protein) {
        # One row per protein
        params <- tibble(
            protein = data[["protein"]],
            param = data[["param"]],
            df = data[["df"]]
        )
    } else {
        # One row per site
        params <- tibble(
            protein = data[["protein"]],
            site = data[["site"]],
            param = data[["param"]],
            df = data[["df"]]
        )
    }

    res <- vector("list", length(controls))
    for (i in seq_along(controls)) {
        # Test for one contrast
        tests <- vector("list", nrow(params))
        ctrl <- controls[i]
        case <- cases[i]
        for (j in seq_along(params$param)) {
            tests[[j]] <- .onetest(params$param[[j]], params$df[j], ctrl, case)
        }
        onectrx <- bind_rows(tests)
        onectrx$Protein <- params$protein
        if (!per_protein) {
            onectrx$Site <- params$site
        }
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


#' Protein-level adjustment
#' @keywords internal
.adjustProteinLevel <- function(diffSite, diffProtein) {
    diffRef <- diffProtein[, c("Protein", "log2FC", "SE", "DF")]
    names(diffRef)[names(diffRef) == "log2FC"] <- "log2FC_ref"
    names(diffRef)[names(diffRef) == "SE"] <- "SE_ref"
    names(diffRef)[names(diffRef) == "DF"] <- "DF_ref"
    joined <- inner_join(diffSite, diffRef)

    log2fc <- joined$log2FC - joined$log2FC_ref
    s2 <- joined$SE ^ 2
    s2_ref <- joined$SE_ref ^ 2
    stderr <- sqrt(s2 + s2_ref)
    numer <- (s2 + s2_ref) ^ 2
    denom <- (s2 ^ 2 / joined$DF + s2_ref ^ 2 / joined$DF_ref)
    df <- numer / denom
    tval <- log2fc / stderr
    pval <- 2 * stats::pt(abs(tval), df, lower.tail = FALSE)

    tibble(Protein = joined$Protein,
           Site = joined$Site,
           Label = joined$Label,
           log2FC = log2fc,
           SE = stderr,
           Tvalue = tval,
           DF = df,
           pvalue = pval)
}
