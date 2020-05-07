
# Summarization -----------------------------------------------------------

#' Summarize feature log-intensities for each PTM site per run.
#'
#' \code{summarize_ptm} summarizes feature log-intensities for each PTM site.
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr nest unnest
#' @importFrom dplyr filter select mutate group_by left_join distinct
#' @importFrom purrr map
#' @param df_site A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{feature}, \code{log2inty}, and possibly,
#'   \code{batch}.
#' @param method A string defining the summarization method. Default is
#'   \code{"tmp"}, which applies Tukey's median polish.
#' @param impute_missing A logical. \code{TRUE} performs missing value
#'   imputation with the accelerated failure time model, \code{FALSE} does not
#'   handle missing values. Default is \code{FALSE}.
#' @return A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{log2inty}, and possibly, \code{batch}.
#' @export
#'
#' @examples
#' summarize_ptm(df_site)
summarize_ptm <- function(df_site, method = "tmp", impute_missing = FALSE) {

    if (missing(df_site))
        stop("Input data frame is missing!")
    if (!is.data.frame(df_site))
        stop("Please provide the input peak list as a data frame in 'df_site'!")
    cols <- names(df_site)
    if (!all(c("protein", "site", "group", "run", "feature", "log2inty") %in% cols)) {
        stop("Missing columns! Please make sure the input data frame contains
             columns of 'protein', 'site', 'group', 'run', 'feature', 'log2inty'")
    }

    df_site <- df_site %>%
        filter(!is.na(log2inty))

    # Experimental design
    design <- df_site %>%
        distinct(run, group)

    # Nested data frame
    if ("batch" %in% cols) {
        nested_site <- df_site %>%
            select(protein, site, batch, run, feature, log2inty) %>%
            nest(data = one_of("run", "feature", "log2inty"))
    } else {
        nested_site <- df_site %>%
            select(protein, site, run, feature, log2inty) %>%
            nest(data = one_of("run", "feature", "log2inty"))
    }
    # if ("batch" %in% cols) {
    #     nested_site <- df_site %>%
    #         select(protein, site, batch, run, feature, log2inty) %>%
    #         group_by(protein, site, batch) %>%
    #         nest()
    # } else {
    #     nested_site <- df_site %>%
    #         select(protein, site, run, feature, log2inty) %>%
    #         group_by(protein, site) %>%
    #         nest()
    # }

    if (!impute_missing) {
        # No missing value imputation
        nested_site <- nested_site %>%
            mutate(sumdata = map(data, summarize_feature, method))
    } else {
        nested_site <- nested_site %>%
            mutate(data = map(data, ~ complete(., feature, run)))
        # Missing value imputation with AFT model
        nested_site <- nested_site %>%
            mutate(aftdata = map(data, fill_censored_aft))
        nested_site <- nested_site %>%
            mutate(sumdata = map(aftdata, summarize_feature, method))
    }
    # Run-level summaries with grouping information
    if (!impute_missing) {
        df_sum <- nested_site %>%
            select(-data) %>%
            unnest(one_of("sumdata")) %>%
            left_join(design)
    } else {
        df_sum <- nested_site %>%
            select(-data, aftdata) %>%
            unnest(one_of("sumdata")) %>%
            left_join(design)
    }
    # df_sum <- nested_site %>%
    #     unnest(sumdata) %>%
    #     left_join(design)

    return(df_sum)
}


#' Summarization of feature intensities.
#'
#' \code{summarize_feature} summarizes feature intensities and returns run-level
#'   summarization. Tukey's median polish is used by default, log (base 2) of
#'   the summed feature intensity is used otherwise.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select group_by summarise
#' @importFrom tidyr spread
#' @importFrom tibble data_frame
#' @param df_prot A data frame.
#' @param method A string. Default is \code{"tmp"}.
#' @return A tibble restoring the run-level summarization.
#' @export
#'
#' @examples
#' summarize_feature(df_prot)
#' summarize_feature(df_prot, method = "tmp")
summarize_feature <- function(df_prot, method = "tmp") {
    if (method == "tmp") {
        inty_wide <- df_prot %>%
            select(feature, run, log2inty) %>%
            spread(feature, log2inty)
        inty_mat <- data.matrix(inty_wide[, -1])
        mp_out <- medpolish(inty_mat, na.rm = TRUE, trace.iter = FALSE)
        df_sum <- data_frame(run = inty_wide$run,
                             log2inty = mp_out$overall + mp_out$row)
    } else {
        df_sum <- df_prot %>%
            group_by(run) %>%
            summarise(log2inty = log2(sum(2 ^ log2inty, na.rm = TRUE)))
    }

    return(df_sum)
}


# Modeling ----------------------------------------------------------------

#' Whole-plot modeling for all sites.
#'
#' \code{model_ptm} fits a whole-plot model for all sites in consideration of
#'   the experimental design and the underlying protein abundance, and returns
#'   the estimates of model parameters.
#'
#' @param df_sum A data frame with columns of \code{protein}, \code{site},
#'   \code{group}, \code{run}, \code{log2inty}, and possibly, \code{batch}.
#' @param w_batch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return A data frame.
#' @export
#'
#' @examples
#' model_ptm(df_sum, w_batch = FALSE)
model_ptm <- function(df_sum, w_batch = FALSE) {

    if (missing(df_sum))
        stop("Input data frame is missing!")
    if (!is.data.frame(df_sum))
        stop("Please provide protein summaries as a data frame in 'df_sum'!")
    cols <- names(df_sum)
    if (!all(c("protein", "site", "group", "run", "log2inty") %in% cols)) {
        stop("Missing columns! Please make sure the input data frame contains
             columns of 'protein', 'site', 'group', 'run', 'log2inty'")
    }
    if (w_batch && !("batch" %in% cols))
        stop("To account for batch effect, add a batch column in df_sum")

    nested <- nest_site(df_sum, w_batch)
    params <- extract_param(nested)

    return(params)
}


#' Whole-plot modeling for all sites.
#'
#' \code{model_ptm2} fits a whole-plot model for all sites in consideration of
#'   the experimental design, and returns the estimates of model parameters.
#'
#' @param df_sum A data frame.
#' @param w_batch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return A data frame.
#' @export
#'
#' @examples
#' model_ptm(df_sum, w_batch = FALSE)
model_ptm2 <- function(df_sum, w_batch = FALSE) {
    nested <- nest_site(df_sum, w_batch)
    params <- extract_param2(nested)

    return(params)
}


#' Whole-plot modeling.
#'
#' \code{nest_site} fits and returns a whole-plot model for each site in
#'   consideration of the experimental design, organized in a nested data frame.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate group_by n_distinct
#' @importFrom tidyr nest
#' @importFrom purrr map map2 map_dbl map_lgl map_int possibly
#' @param df_sum A data frame.
#' @param w_batch A logical. \code{TRUE} considers batch effect, \code{FALSE}
#'   otherwise. Default is \code{FALSE}.
#' @return An instance of nested data frame.
#' @export
#'
#' @examples
#' nest_site(df_sum, w_batch = FALSE)
nest_site <- function(df_sum, w_batch = FALSE) {
    if (w_batch) {
        if (!("batch" %in% names(df_sum)))
            stop("There is no information about batch!")

        # One model for all batches (need data with >1 batches)
        nested_sum <- df_sum %>%
            nest(data = -one_of("protein", "site")) %>%
            mutate(nb_bch = map_int(data, ~ n_distinct(.$batch))) %>%
            filter(nb_bch > 1) %>%
            select(-nb_bch)
        # nested_sum <- df_sum %>%
        #     group_by(protein, site) %>%
        #     nest_legacy() %>%
        #     mutate(nb_bch = map_int(data, ~ n_distinct(.$batch))) %>%
        #     filter(nb_bch > 1) %>%
        #     select(-nb_bch)

        # Linear model and associated parameter estimates
        nested_sum <- nested_sum %>%
            mutate(lm_fit = map(data, possibly(lm_group, otherwise = NA_real_), w_batch)) %>%
            filter(!is.na(lm_fit)) %>%
            mutate(
                param = map2(lm_fit, data, tidy_bch),
                df_res = map_dbl(lm_fit, df.residual)
            )
        # nested_sum <- nested_sum %>%
        #     mutate(lm_fit = map(data, lm_group, w_batch)) %>%
        #     mutate(
        #         param = map2(lm_fit, data, tidy_bch),
        #         df_res = map_dbl(lm_fit, df.residual)
        #     )
        # Remove cases not eligible for hypothesis testing (SE is NA)
        # [TODO]: consider another option to report FC alone
        nested_sum <- nested_sum %>%
            filter(df_res > 0, !map_lgl(param, ~any(is.na(.$std.error))))
    } else {
        # One model per site (and potentially batch)
        if ("batch" %in% names(df_sum)) {
            nested_sum <- df_sum %>%
                nest(data = -one_of("protein", "site", "batch"))
            # nested_sum <- df_sum %>%
            #     group_by(protein, site, batch) %>%
            #     nest_legacy()
        } else {
            nested_sum <- df_sum %>%
                nest(data = -one_of("protein", "site"))
            # nested_sum <- df_sum %>%
            #     group_by(protein, site) %>%
            #     nest_legacy()
        }
        # Linear model and associated parameter estimates
        nested_sum <- nested_sum %>%
            mutate(lm_fit = map(data, possibly(lm_group, otherwise = NA_real_), w_batch)) %>%
            filter(!is.na(lm_fit)) %>%
            mutate(
                param = map2(lm_fit, data, tidy_bch),
                df_res = map_dbl(lm_fit, df.residual)
            )
        # nested_sum <- nested_sum %>%
        #     mutate(lm_fit = map(data, lm_group, w_batch)) %>%
        #     mutate(
        #         param = map2(lm_fit, data, tidy_bch),
        #         df_res = map_dbl(lm_fit, df.residual)
        #     )
        # Remove cases not eligible for hypothesis testing (SE is NA)
        # [TODO]: consider another option to report FC alone
        nested_sum <- nested_sum %>%
            filter(df_res > 0, !map_lgl(param, ~any(is.na(.$std.error))))
    }

    return(nested_sum)
}


#' Extract estimated parameters.
#'
#' \code{extract_param} extracts the estimated model parameters from nested data
#'   frame.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select rename left_join
#' @importFrom tidyr unnest
#' @param nested An instance of nested data frame.
#' @return A data frame restoring the estimated model parameters.
#' @export
#'
#' @examples
#' extract_param(nested)
extract_param <- function(nested) {
    if ("batch" %in% names(nested)) {
        # per-batch model
        nested_param <- nested %>%
            select(protein, site, batch, param, df_res)
    } else {
        # all-batch model
        nested_param <- nested %>%
            select(protein, site, param, df_res)
    }
    param_mod <- nested_param %>%
        filter(site != "None") %>%
        unnest(one_of("param"))
    param_unmod <- nested_param %>%
        filter(site == "None") %>%
        select(-site) %>%
        unnest(one_of("param")) %>%
        rename(df_unmod = df_res, est_unmod = estimate, se_unmod = std.error)

    return(left_join(param_mod, param_unmod))
}


#' Extract estimated parameters.
#'
#' \code{extract_param2} extracts the estimated model parameters from nested data
#'   frame.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom tidyr unnest
#' @param nested An instance of nested data frame.
#' @return A data frame restoring the estimated model parameters.
#' @export
#'
#' @examples
#' extract_param(nested)
extract_param2 <- function(nested) {
    if ("batch" %in% names(nested)) {
        # per-batch model
        nested_param <- nested %>%
            select(protein, site, batch, param, df_res)
    } else {
        # all-batch model
        nested_param <- nested %>%
            select(protein, site, param, df_res)
    }
    param_mod <- nested_param %>%
        unnest(param)

    return(param_mod)
}


# Comparison --------------------------------------------------------------

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


#' Compare abundance of a modified site across conditions.
#'
#' \code{compare_mod} performs significance analysis for differentially modified
#'   sites across conditions.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate group_by summarise ungroup bind_rows
#'   distinct anti_join
#' @importFrom tidyr unite separate complete
#' @param df_mod A data frame.
#' @param grp_ctrl Control group as reference.
#' @param grp_case Case group.
#' @param protadj A logical. \code{TRUE} performs protein-level adjustment,
#'   \code{FALSE} otherwise. Default is \code{TRUE}.
#' @return A data frame.
#' @export
#'
#' @examples
#' compare_mod(df_mod, grp_ctrl, grp_case, protadj = TRUE)
compare_mod2 <- function(df_mod, grp_ctrl, grp_case, protadj = TRUE) {
    df_mod <- df_mod %>%
        unite(protsite, protein, site, sep = "--")

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
            ) %>%
            ungroup()
        res <- diff_mod %>%
            mutate(
                log2FC = log2fc - log2fc_unmod, std_error = sqrt(se2 + se2_unmod),
                DF = (se2 + se2_unmod) ^ 2 / (se2 ^ 2 / df_res + se2_unmod ^ 2 / df_unmod),
                statistic = log2FC / std_error, p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)
            ) %>%
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
            complete(protsite, key_grp)
        # Model-based inference of log-difference
        diff_mod <- full_mod %>%
            group_by(protsite) %>%
            filter(!any(is.na(df_res))) %>%
            group_by(protsite, df_res) %>%
            summarise(log2fc = diff(estimate), se2 = sum(std.error ^ 2)) %>%
            ungroup()
        res <- diff_mod %>%
            mutate(
                log2FC = log2fc, std_error = sqrt(se2), DF = df_res, statistic = log2FC / std_error,
                p_value = 2 * pt(abs(statistic), df = DF, lower.tail = FALSE)
            ) %>%
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
                ) %>%
                select(-key_grp)
            res <- res %>% bind_rows(untest_mod)
        }
    }

    return(res %>% separate(protsite, into = c("protein", "site"), sep = "--"))
}
