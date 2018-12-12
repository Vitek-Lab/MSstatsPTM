#' Normalize feature intensities across runs.
#'
#' \code{normalize_ptm} normalizes feature intensities across runs in a PTM
#'   dataset, using unpaired unmodified peptides.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select mutate group_by summarise ungroup left_join
#' @param data A data frame.
#' @return A data frame.
#' @export
#'
#' @examples
#' normalize_ptm(data)
normalize_ptm <- function(data) {
    cols <- names(data)
    # Based on unmodified peptides
    if ("batch" %in% cols) {
        medians <- data %>%
            filter(!is_mod) %>%
            group_by(batch, run) %>%
            summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>%
            mutate(log2inty_bch = median(log2inty_med)) %>%
            ungroup()
    } else {
        medians <- data %>%
            filter(!is_mod) %>%
            group_by(run) %>%
            summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>%
            mutate(log2inty_bch = median(log2inty_med)) %>%
            ungroup()
    }
    normdata <- data %>%
        left_join(medians) %>%
        mutate(log2inty = ifelse(is.na(log2inty_med), log2inty,
                                 log2inty - log2inty_med + log2inty_bch)) %>%
        select(-log2inty_med, -log2inty_bch)

    return(normdata)
}

