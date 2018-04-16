
#' Normalize feature intensities across runs.
#'
#' \code{normalize_ptm} normalizes feature intensities across runs in a PTM
#'   dataset, using unpaired unmodified peptides.
#'
#' @param data A data frame.
#' @return A data frame.
#'
#' @examples
#' normalize_ptm(data)
#'
#' @export
#'
normalize_ptm <- function(data) {
    cols <- names(data)
    # Based on unmodified peptides
    if ("batch" %in% cols) {
        medians <- data %>%
            dplyr::filter(!is_mod) %>%
            dplyr::group_by(batch, run) %>%
            dplyr::summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>%
            dplyr::mutate(log2inty_bch = median(log2inty_med)) %>%
            dplyr::ungroup()
    } else {
        medians <- data %>%
            dplyr::filter(!is_mod) %>%
            dplyr::group_by(run) %>%
            dplyr::summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>%
            dplyr::mutate(log2inty_bch = median(log2inty_med)) %>%
            dplyr::ungroup()
    }
    normdata <- data %>%
        dplyr::left_join(medians) %>%
        dplyr::mutate(log2inty = ifelse(is.na(log2inty_med), log2inty,
                                        log2inty - log2inty_med + log2inty_bch)) %>%
        dplyr::select(-log2inty_med, -log2inty_bch)

    return(normdata)
}

