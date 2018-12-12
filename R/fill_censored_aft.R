#' Annotate censored values and fill with AFT.
#'
#' \code{fill_censored_aft} annotates censored feature intensities and imputes
#'   the censored values using the accelerated failure time model.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate group_by ungroup n_distinct
#' @importFrom survival survreg Surv
#' @param data A data frame.
#' @return A data frame.
#' @export
#'
#' @examples
#' fill_censored_aft(data)
fill_censored_aft <- function(data) {
    if (n_distinct(data$feature) == 1) return(data)  # only one feature
    if (all(!is.na(data$log2inty))) return(data)  # no missing value
    # Annotate observation status and add censored values
    aftdata <- data %>%
        mutate(ind_obs = ifelse(is.na(log2inty), 0L, 1L)) %>%
        group_by(feature) %>%
        mutate(log2inty_aft = ifelse(
            ind_obs == 0, 0.99 * min(log2inty[ind_obs == 1]), log2inty
        )) %>%
        ungroup()
    # AFT model with effects of run and feature
    fit <- tryCatch(
        survreg(Surv(log2inty_aft, ind_obs, type = "left") ~ run + feature,
                data = aftdata, dist = "gaussian"),
        error = function(e) e, warning = function(w) w
    )
    if (is(fit, "warning")) return(data)  # not converged
    aftdata <- aftdata %>%
        mutate(
            log2inty_pred = predict(fit),
            log2inty = ifelse(ind_obs == 0, pmin(log2inty_pred, log2inty_aft), log2inty_aft)
        ) %>%
        select(-log2inty_aft, -log2inty_pred, -ind_obs)

    return(aftdata)
}
