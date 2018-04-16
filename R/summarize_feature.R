
#' Summarization of feature intensities.
#'
#' \code{summarize_feature} summarizes feature intensities and returns run-level
#'   summarization. Tukey's median polish is used by default, log (base 2) of
#'   the summed feature intensity is used otherwise.
#'
#' @param df_prot A data frame.
#' @param method A string. Default is \code{tmp}.
#' @return A tibble restoring the run-level summarization.
#'
#' @examples
#' summarize_feature(df_prot)
#' summarize_feature(df_prot, method = "tmp")
#'
#' @export
#'
summarize_feature <- function(df_prot, method = "tmp") {
    if (method == "tmp") {
        inty_wide <- df_prot %>%
            dplyr::select(feature, run, log2inty) %>%
            tidyr::spread(feature, log2inty)
        inty_mat <- data.matrix(inty_wide[, -1])
        mp_out <- medpolish(inty_mat, na.rm = TRUE, trace.iter = FALSE)
        df_sum <- dplyr::data_frame(run = inty_wide$run,
                                    log2inty = mp_out$overall + mp_out$row)
    } else {
        df_sum <- df_prot %>%
            dplyr::group_by(run) %>%
            dplyr::summarise(log2inty = log2(sum(2 ^ log2inty, na.rm = TRUE)))
    }

    return(df_sum)
}
