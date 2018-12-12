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
