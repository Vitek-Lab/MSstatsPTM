#' Compare mean abundances for all PTM sites across conditions.
#'
#' \code{PTMcompareMeans} performs significance analysis for detection of
#' changes in PTM mean abundances between conditions.
#'
#' @param data A list that contain estimates for PTM site abundance
#'   (\code{data$PTM}), and possibly, protein abundance (\code{data$Protein}).
#' @param controls A string vector of control groups as reference.
#' @param cases A string vector of case groups.
#' @param adjProtein A logical. \code{TRUE} performs protein-level adjustment,
#'   \code{FALSE} otherwise. Default is \code{FALSE}.
#' @return A data frame.
#'
#' @export
#' @examples
#' PTMcompareMeans(df_mod, controls, cases)
PTMcompareMeans <- function(data, controls, cases, adjProtein = FALSE) {
    if (is.null(data[["PTM"]]))
        stop("PTM estimates are not provided")
    if (is.null(data[["Protein"]]) && adjProtein)
        stop("To adjust for protein level, protein estimates are required")

    res <- extractMeanDiff(data[["PTM"]], controls, cases)
    if (adjProtein) {
        res_ref <- extractMeanDiff(data[["Protein"]], controls, cases)
        res <- adjustProteinLevel(res, res_ref)
    }

    res
}


extractMeanDiff <- function(data, controls, cases) {
    # Params per site in a row
    params <- tibble(
        protein = data[["protein"]],
        site = data[["site"]],
        param = data[["param"]],
        df = data[["df"]]
    )

    res <- vector("list", length(controls))
    for (i in seq_along(controls)) {
        onectrx <- vector("list", nrow(params))
        for (j in seq_along(params$param)) {
            onectrx[[j]] <- tibble(
                Protein = params$protein[j],
                Site = params$site[j],
                onetest(params$param[[j]], params$df[j], controls[i], cases[i])
            )
        }
        res[[i]] <- bind_rows(onectrx)
    }
    bind_rows(res)
}

onetest <- function(param, df, ctrl, case) {
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
        pval <- 2 * pt(abs(tval), df, lower.tail = FALSE)
        res <- tibble(Label = paste(case, ctrl, sep = " vs "), log2FC = log2fc,
                      SE = stderr, Tvalue = tval, DF = df, pvalue = pval)
    }
    res
}


#' Protein-level adjustment.
adjustProteinLevel <- function(data, refdata) {

    refdata <- dplyr::select(refdata, Protein, log2FC_ref = log2FC, SE_ref = SE, DF_ref = DF)
    joined <- dplyr::inner_join(data, refdata)

    log2fc <- joined$log2FC - joined$log2FC_ref
    s2 <- joined$SE ^ 2
    s2_ref <- joined$SE_ref ^ 2
    stderr <- sqrt(s2 + s2_ref)
    numer <- (s2 + s2_ref) ^ 2
    denom <- (s2 ^ 2 / joined$DF + s2_ref ^ 2 / joined$DF_ref)
    df <- numer / denom

    tval <- log2fc / stderr
    pval <- 2 * pt(abs(tval), df, lower.tail = FALSE)
    res <- tibble(Protein = joined$Protein,
                  Site = joined$Site,
                  Label = joined$Label,
                  log2FC = log2fc,
                  SE = stderr,
                  Tvalue = tval,
                  DF = df,
                  pvalue = pval)
    res
}
