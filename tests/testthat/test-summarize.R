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

test_that("PTMsummarize() takes and returns the right format", {
    s <- PTMsimulateExperiment(
        nGroup=2, nRep=2, nProtein=1, nSite=2, nFeature=5,
        logAbundance=list(
            PTM=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05),
            PROTEIN=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05)
        )
    )
    ns <- PTMnormalize(s)

    expect_error(PTMnormalize(ns[["PTM"]]), "PTM peak list is missing!")
    expect_output(str(PTMnormalize(ns)), "List of 2")
    expect_identical(names(PTMnormalize(ns)), c("PTM", "PROTEIN"))

    ns[["PROTEIN"]] <- NULL
    expect_output(str(PTMnormalize(ns)), "List of 1")
    expect_identical(names(PTMnormalize(ns)), "PTM")
})

