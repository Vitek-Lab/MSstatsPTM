
#' Annotate censored values and fill with AFT.
#'
#' \code{fill_censored_aft} annotates censored feature intensities and imputes
#'   the censored values using the accelerated failure time model.
#'
#' @param data A data frame.
#' @return A data frame.
#'
#' @examples
#' fill_censored_aft(data)
#'
#' @export
#'
fill_censored_aft <- function(data) {
    if (dplyr::n_distinct(data$feature) == 1) return(data)  # only one feature
    if (all(!is.na(data$log2inty))) return(data)  # no missing value
    # Annotate observation status and add censored values
    aftdata <- data %>%
        dplyr::mutate(ind_obs = ifelse(is.na(log2inty), 0L, 1L)) %>%
        dplyr::group_by(feature) %>%
        dplyr::mutate(log2inty_aft = ifelse(
            ind_obs == 0, 0.99 * min(log2inty[ind_obs == 1]), log2inty
        )) %>%
        dplyr::ungroup()
    # AFT model with effects of run and feature
    fit <- tryCatch(
        survival::survreg(survival::Surv(log2inty_aft, ind_obs, type = "left") ~ run + feature,
                data = aftdata, dist = "gaussian"),
        error = function(e) e, warning = function(w) w
    )
    if (is(fit, "warning")) return(data)  # not converged
    aftdata <- aftdata %>%
        dplyr::mutate(
            log2inty_pred = predict(fit),
            log2inty = ifelse(ind_obs == 0, pmin(log2inty_pred, log2inty_aft), log2inty_aft)
        ) %>%
        dplyr::select(-log2inty_aft, -log2inty_pred, -ind_obs)

    return(aftdata)
}
