
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
compare_mod <- function(df_mod, grp_ctrl, grp_case, protadj = TRUE) {
    if ("batch" %in% names(df_mod)) {
        # Batch-aggregated testing (with multiple difference estimates)
        if (protadj) {
            # With protein-level adjustment
            full_mod <- df_mod %>%
                filter(group %in% c(grp_ctrl, grp_case)) %>%
                filter(!is.na(df_res) | !is.na(df_unmod)) %>%
                mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>%
                complete(protsite, batch, key_grp)
            # Model-based inference of log-difference
            diff_mod <- full_mod %>%
                group_by(protsite) %>%
                filter(all(!is.na(df_res)), all(!is.na(df_unmod))) %>%
                group_by(protsite, batch, df_res, df_unmod) %>%
                summarise(
                    log2fc = diff(estimate), se2 = sum(std.error ^ 2),
                    log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
                ) %>% ungroup()
            res <- diff_mod %>%
                group_by(protsite) %>%
                summarise(
                    log2FC = mean(log2fc) - mean(log2fc_unmod),
                    std_error = sqrt(sum(se2) + sum(se2_unmod)) / n(),
                    DF = (sum(se2) + sum(se2_unmod)) ^ 2 / sum(se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod),
                    statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)
                ) %>%
                mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
            # Missing in one group
            part_mod <- full_mod %>%
                anti_join(diff_mod %>% select(protsite)) %>%
                group_by(protsite) %>%
                filter(all(is.na(df_res) == is.na(df_unmod))) %>%
                filter(
                    n_distinct(is.na(df_res[key_grp == "G0"])) == 1,
                    n_distinct(is.na(df_res[key_grp == "G1"])) == 1
                ) %>%
                ungroup() %>%
                filter(!is.na(df_res)) %>%
                distinct(protsite, key_grp)
            if (nrow(part_mod) > 0) {
                untest_mod <- part_mod %>%
                    mutate(
                        log2FC = ifelse(key_grp == "G1", Inf, -Inf),
                        std_error = NA, DF = NA, statistic = NA, p_value = NA,
                        contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                    ) %>% select(-key_grp)
                res <- res %>% bind_rows(untest_mod)
            }
        } else {
            # Without protein-level adjustment
            full_mod <- df_mod %>%
                filter(group %in% c(grp_ctrl, grp_case)) %>%
                filter(!is.na(df_res)) %>%
                mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>%
                complete(protsite, batch, key_grp)
            # Model-based inference of log-difference
            diff_mod <- full_mod %>%
                group_by(protsite) %>%
                filter(all(!is.na(df_res))) %>%
                group_by(protsite, batch, df_res) %>%
                summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>%
                ungroup()
            res <- diff_mod %>%
                group_by(protsite) %>%
                summarise(log2FC = mean(log2fc), std_error = sqrt(sum(se2)) / n(),
                          DF = sum(se2) ^ 2 / sum(se2 ^ 2 / df_res),
                          statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)) %>%
                mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
            # Missing in one group
            part_mod <- full_mod %>%
                anti_join(diff_mod %>% select(protsite)) %>%
                group_by(protsite) %>%
                filter(
                    n_distinct(is.na(df_res[key_grp == "G0"])) == 1,
                    n_distinct(is.na(df_res[key_grp == "G1"])) == 1
                ) %>%
                ungroup() %>%
                filter(!is.na(df_res)) %>%
                distinct(protsite, key_grp)
            if (nrow(part_mod) > 0) {
                untest_mod <- part_mod %>%
                    mutate(
                        log2FC = ifelse(key_grp == "G1", Inf, -Inf),
                        std_error = NA, DF = NA, statistic = NA, p_value = NA,
                        contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                    ) %>% select(-key_grp)
                res <- res %>% bind_rows(untest_mod)
            }
        }
    } else {
        # Single difference estimate between groups
        if (protadj) {
            # With protein-level adjustment
            full_mod <- df_mod %>%
                filter(group %in% c(grp_ctrl, grp_case)) %>%
                filter(!is.na(df_res) | !is.na(df_unmod)) %>%
                mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>%
                complete(protsite, key_grp)
            # Model-based inference of log-difference
            diff_mod <- full_mod %>%
                group_by(protsite) %>%
                filter(all(!is.na(df_res)), all(!is.na(df_unmod))) %>%
                group_by(protsite, df_res, df_unmod) %>%
                summarise(
                    log2fc = diff(estimate), se2 = sum(std.error ^ 2),
                    log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
                ) %>% ungroup()
            res <- diff_mod %>%
                mutate(log2FC = log2fc - log2fc_unmod, std_error = sqrt(se2 + se2_unmod),
                       DF = (se2 + se2_unmod) ^ 2 / (se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod),
                       statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)) %>%
                select(protsite, log2FC, std_error, DF, statistic, p_value) %>%
                mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
            # Missing in one group
            part_mod <- full_mod %>%
                anti_join(diff_mod %>% select(protsite)) %>%
                group_by(protsite) %>%
                filter(all(is.na(df_res) == is.na(df_unmod))) %>%
                ungroup() %>%
                filter(!is.na(df_res)) %>%
                distinct(protsite, key_grp)
            if (nrow(part_mod) > 0) {
                untest_mod <- part_mod %>%
                    mutate(
                        log2FC = ifelse(key_grp == "G1", Inf, -Inf),
                        std_error = NA, DF = NA, statistic = NA, p_value = NA,
                        contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                    ) %>% select(-key_grp)
                res <- res %>% bind_rows(untest_mod)
            }
        } else {
            # Without protein-level adjustment
            full_mod <- df_mod %>%
                filter(group %in% c(grp_ctrl, grp_case)) %>%
                filter(!is.na(df_res)) %>%
                mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>%
                complete(protsite, key_grp)
            # Model-based inference of log-difference
            diff_mod <- full_mod %>%
                group_by(protsite) %>%
                filter(!any(is.na(df_res))) %>%
                group_by(protsite, df_res) %>%
                summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>%
                ungroup()
            res <- diff_mod %>%
                mutate(log2FC = log2fc, std_error = sqrt(se2), DF = df_res,
                       statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)) %>%
                select(protsite, log2FC, std_error, DF, statistic, p_value) %>%
                mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
            # Missing in one group
            part_mod <- full_mod %>%
                anti_join(diff_mod %>% select(protsite)) %>%
                filter(!is.na(df_res)) %>%
                distinct(protsite, key_grp)
            if (nrow(part_mod) > 0) {
                untest_mod <- part_mod %>%
                    mutate(
                        log2FC = ifelse(key_grp == "G1", Inf, -Inf),
                        std_error = NA, DF = NA, statistic = NA, p_value = NA,
                        contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                    ) %>% select(-key_grp)
                res <- res %>% bind_rows(untest_mod)
            }
        }
    }

    return(res)
}

