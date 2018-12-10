
#' Compare abundance of a modified site across conditions.
#'
#' \code{compare_mod} performs significance analysis for differentially modified
#'   sites across conditions.
#'
#' @param df_mod A data frame.
#' @param grp_ctrl Control group as reference.
#' @param grp_case Case group.
#' @param protadj A logical. \code{TRUE} performs protein-level adjustment,
#'   \code{FALSE} otherwise. Default is \code{TRUE}.
#' @return A data frame.
#'
#' @examples
#' compare_mod(df_mod, grp_ctrl, grp_case, protadj = TRUE)
#'
#' @export
#'
compare_mod2 <- function(df_mod, grp_ctrl, grp_case, protadj = TRUE) {
    df_mod <- df_mod %>%
        tidyr::unite(protsite, protein, site, sep = "--")

    if (protadj) {
        # With protein-level adjustment
        full_mod <- df_mod %>%
            dplyr::filter(group %in% c(grp_ctrl, grp_case)) %>%
            dplyr::filter(!is.na(df_res) | !is.na(df_unmod)) %>%
            dplyr::mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>%
            tidyr::complete(protsite, key_grp)
        # Model-based inference of log-difference
        diff_mod <- full_mod %>%
            dplyr::group_by(protsite) %>%
            dplyr::filter(all(!is.na(df_res)), all(!is.na(df_unmod))) %>%
            dplyr::group_by(protsite, df_res, df_unmod) %>%
            dplyr::summarise(
                log2fc = diff(estimate), se2 = sum(std.error ^ 2),
                log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
            ) %>%
            dplyr::ungroup()
        res <- diff_mod %>%
            dplyr::mutate(
                log2FC = log2fc - log2fc_unmod, std_error = sqrt(se2 + se2_unmod),
                DF = (se2 + se2_unmod) ^ 2 / (se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod),
                statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)
            ) %>%
            dplyr::select(protsite, log2FC, std_error, DF, statistic, p_value) %>%
            dplyr::mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
        # Missing in one group
        part_mod <- full_mod %>%
            dplyr::anti_join(diff_mod %>% dplyr::select(protsite)) %>%
            dplyr::group_by(protsite) %>%
            dplyr::filter(all(is.na(df_res) == is.na(df_unmod))) %>%
            dplyr::ungroup() %>%
            dplyr::filter(!is.na(df_res)) %>%
            dplyr::distinct(protsite, key_grp)
        if (nrow(part_mod) > 0) {
            untest_mod <- part_mod %>%
                dplyr::mutate(
                    log2FC = ifelse(key_grp == "G1", Inf, -Inf),
                    std_error = NA, DF = NA, statistic = NA, p_value = NA,
                    contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                ) %>%
                dplyr::select(-key_grp)
            res <- res %>% dplyr::bind_rows(untest_mod)
        }
    } else {
        # Without protein-level adjustment
        full_mod <- df_mod %>%
            dplyr::filter(group %in% c(grp_ctrl, grp_case)) %>%
            dplyr::filter(!is.na(df_res)) %>%
            dplyr::mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>%
            tidyr::complete(protsite, key_grp)
        # Model-based inference of log-difference
        diff_mod <- full_mod %>%
            dplyr::group_by(protsite) %>%
            dplyr::filter(!any(is.na(df_res))) %>%
            dplyr::group_by(protsite, df_res) %>%
            dplyr::summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>%
            dplyr::ungroup()
        res <- diff_mod %>%
            dplyr::mutate(
                log2FC = log2fc, std_error = sqrt(se2), DF = df_res, statistic = log2FC / std_error,
                p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)
            ) %>%
            dplyr::select(protsite, log2FC, std_error, DF, statistic, p_value) %>%
            dplyr::mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
        # Missing in one group
        part_mod <- full_mod %>%
            dplyr::anti_join(diff_mod %>% dplyr::select(protsite)) %>%
            dplyr::filter(!is.na(df_res)) %>%
            dplyr::distinct(protsite, key_grp)
        if (nrow(part_mod) > 0) {
            untest_mod <- part_mod %>%
                dplyr::mutate(
                    log2FC = ifelse(key_grp == "G1", Inf, -Inf),
                    std_error = NA, DF = NA, statistic = NA, p_value = NA,
                    contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                ) %>%
                dplyr::select(-key_grp)
            res <- res %>% dplyr::bind_rows(untest_mod)
        }
    }

    return(res %>% tidyr::separate(protsite, into = c("protein", "site"), sep = "--"))
}
