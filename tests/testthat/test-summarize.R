test_that("summarizeFeatures() uses different methods for summarization", {

    df <- tibble(
        run = c("a", "a", "a", "b", "b"),
        feature = c("F1", "F2", "F3", "F1", "F3"),
        log2inty = c(1, 2, 3, NA, 3)
    )
    res_max    <- tibble(run = c("a", "b"), log2inty = c(3, 3))
    res_mean   <- tibble(run = c("a", "b"), log2inty = c(2, 3))
    res_median <- tibble(run = c("a", "b"), log2inty = c(2, 3))
    res_logsum <- tibble(run = c("a", "b"), log2inty = c(log2(sum(2 ^ c(1, 2, 3))), 3))

    expect_equal(summarizeFeatures(df, method = "max"), res_max)
    expect_equal(summarizeFeatures(df, method = "mean"), res_mean)
    expect_equal(summarizeFeatures(df, method = "median"), res_median)
    expect_equal(summarizeFeatures(df, method = "logsum"), res_logsum)
})
