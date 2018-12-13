#' Compare abundance of a modified site across conditions.
#'
#' \code{compare_mod} performs significance analysis for differentially modified
#'   sites across conditions.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate group_by summarise ungroup bind_rows
#'   n_distinct distinct anti_join
#' @importFrom tidyr unite separate complete
#' @param df_mod A data frame.
#' @param controls Control group as reference.
#' @param cases Case group.
#' @param protadj A logical. \code{TRUE} performs protein-level adjustment,
#'   \code{FALSE} otherwise. Default is \code{TRUE}.
#' @return A data frame.
#' @export
#'
#' @examples
#' compare_mod(df_mod, controls, cases, protadj = TRUE)
compare_mod <- function(df_mod, controls, cases, protadj = TRUE) {
    df_mod <- df_mod %>%
        unite(protein_site, protein, site, sep = "--")

    results <- vector("list", length = length(cases))
    for (i in seq_along(cases)) {
        grp_ctrl <- controls[i]
        grp_case <- cases[i]
        if ("batch" %in% names(df_mod)) {
            # Batch-aggregated testing (with multiple difference estimates)
            if (protadj) {
                # With protein-level adjustment
                full_mod <- df_mod %>%
                    filter(group %in% c(grp_ctrl, grp_case)) %>%
                    filter(!is.na(df_res) | !is.na(df_unmod)) %>%
                    mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>%
                    complete(protein_site, batch, key_grp)
                # Model-based inference of log-difference
                diff_mod <- full_mod %>%
                    group_by(protein_site) %>%
                    filter(all(!is.na(df_res)), all(!is.na(df_unmod))) %>%
                    group_by(protein_site, batch, df_res, df_unmod) %>%
                    summarise(
                        log2fc = diff(estimate), se2 = sum(std.error ^ 2),
                        log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
                    ) %>%
                    ungroup()
                res <- diff_mod %>%
                    group_by(protein_site) %>%
                    summarise(
                        log2FC = mean(log2fc) - mean(log2fc_unmod),
                        std_error = sqrt(sum(se2) + sum(se2_unmod)) / n(),
                        DF = (sum(se2) + sum(se2_unmod)) ^ 2 / sum(se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod),
                        statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)
                    ) %>%
                    mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
                # Missing in one group
                part_mod <- full_mod %>%
                    anti_join(diff_mod %>% select(protein_site)) %>%
                    group_by(protein_site) %>%
                    filter(all(is.na(df_res) == is.na(df_unmod))) %>%
                    filter(
                        n_distinct(is.na(df_res[key_grp == "G0"])) == 1,
                        n_distinct(is.na(df_res[key_grp == "G1"])) == 1
                    ) %>%
                    ungroup() %>%
                    filter(!is.na(df_res)) %>%
                    distinct(protein_site, key_grp)
                if (nrow(part_mod) > 0) {
                    untest_mod <- part_mod %>%
                        mutate(
                            log2FC = ifelse(key_grp == "G1", Inf, -Inf),
                            std_error = NA, DF = NA, statistic = NA, p_value = NA,
                            contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                        ) %>%
                        select(-key_grp)
                    res <- res %>% bind_rows(untest_mod)
                }
            } else {
                # Without protein-level adjustment
                full_mod <- df_mod %>%
                    filter(group %in% c(grp_ctrl, grp_case)) %>%
                    filter(!is.na(df_res)) %>%
                    mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>%
                    complete(protein_site, batch, key_grp)
                # Model-based inference of log-difference
                diff_mod <- full_mod %>%
                    group_by(protein_site) %>%
                    filter(all(!is.na(df_res))) %>%
                    group_by(protein_site, batch, df_res) %>%
                    summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>%
                    ungroup()
                res <- diff_mod %>%
                    group_by(protein_site) %>%
                    summarise(
                        log2FC = mean(log2fc), std_error = sqrt(sum(se2)) / n(),
                        DF = sum(se2) ^ 2 / sum(se2 ^ 2 / df_res), statistic = log2FC / std_error,
                        p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)
                    ) %>%
                    mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
                # Missing in one group
                part_mod <- full_mod %>%
                    anti_join(diff_mod %>% select(protein_site)) %>%
                    group_by(protein_site) %>%
                    filter(
                        n_distinct(is.na(df_res[key_grp == "G0"])) == 1,
                        n_distinct(is.na(df_res[key_grp == "G1"])) == 1
                    ) %>%
                    ungroup() %>%
                    filter(!is.na(df_res)) %>%
                    distinct(protein_site, key_grp)
                if (nrow(part_mod) > 0) {
                    untest_mod <- part_mod %>%
                        mutate(
                            log2FC = ifelse(key_grp == "G1", Inf, -Inf),
                            std_error = NA, DF = NA, statistic = NA, p_value = NA,
                            contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                        ) %>%
                        select(-key_grp)
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
                    complete(protein_site, key_grp)
                # Model-based inference of log-difference
                diff_mod <- full_mod %>%
                    group_by(protein_site) %>%
                    filter(all(!is.na(df_res)), all(!is.na(df_unmod))) %>%
                    group_by(protein_site, df_res, df_unmod) %>%
                    summarise(
                        log2fc = diff(estimate), se2 = sum(std.error ^ 2),
                        log2fc_unmod = diff(est_unmod), se2_unmod = sum(se_unmod ^ 2)
                    ) %>%
                    ungroup()
                res <- diff_mod %>%
                    mutate(
                        log2FC = log2fc - log2fc_unmod, std_error = sqrt(se2 + se2_unmod),
                        DF = (se2 + se2_unmod) ^ 2 / (se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod),
                        statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)
                    ) %>%
                    select(protein_site, log2FC, std_error, DF, statistic, p_value) %>%
                    mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
                # Missing in one group
                part_mod <- full_mod %>%
                    anti_join(diff_mod %>% select(protein_site)) %>%
                    group_by(protein_site) %>%
                    filter(all(is.na(df_res) == is.na(df_unmod))) %>%
                    ungroup() %>%
                    filter(!is.na(df_res)) %>%
                    distinct(protein_site, key_grp)
                if (nrow(part_mod) > 0) {
                    untest_mod <- part_mod %>%
                        mutate(
                            log2FC = ifelse(key_grp == "G1", Inf, -Inf),
                            std_error = NA, DF = NA, statistic = NA, p_value = NA,
                            contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                        ) %>%
                        select(-key_grp)
                    res <- res %>% bind_rows(untest_mod)
                }
            } else {
                # Without protein-level adjustment
                full_mod <- df_mod %>%
                    filter(group %in% c(grp_ctrl, grp_case)) %>%
                    filter(!is.na(df_res)) %>%
                    mutate(key_grp = ifelse(group == grp_ctrl, "G0", "G1")) %>%
                    complete(protein_site, key_grp)
                # Model-based inference of log-difference
                diff_mod <- full_mod %>%
                    group_by(protein_site) %>%
                    filter(!any(is.na(df_res))) %>%
                    group_by(protein_site, df_res) %>%
                    summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>%
                    ungroup()
                res <- diff_mod %>%
                    mutate(
                        log2FC = log2fc, std_error = sqrt(se2), DF = df_res, statistic = log2FC / std_error,
                        p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)
                    ) %>%
                    select(protein_site, log2FC, std_error, DF, statistic, p_value) %>%
                    mutate(contrast = paste(grp_case, grp_ctrl, sep = " vs "))
                # Missing in one group
                part_mod <- full_mod %>%
                    anti_join(diff_mod %>% select(protein_site)) %>%
                    filter(!is.na(df_res)) %>%
                    distinct(protein_site, key_grp)
                if (nrow(part_mod) > 0) {
                    untest_mod <- part_mod %>%
                        mutate(
                            log2FC = ifelse(key_grp == "G1", Inf, -Inf),
                            std_error = NA, DF = NA, statistic = NA, p_value = NA,
                            contrast = paste(grp_case, grp_ctrl, sep = " vs ")
                        ) %>%
                        select(-key_grp)
                    res <- res %>% bind_rows(untest_mod)
                }
            }
        }
        results[[i]] <- res
    }
    output <- bind_rows(results) %>%
        separate(protein_site, into = c("protein", "site"), sep = "--")

    return(output)
}
