test_that("simulateSummarization() generates data with correct size", {
    g1r2 <- simulateSummarization(1, 2, 25, 0, 0.2)
    g1r5 <- simulateSummarization(1, 5, 25, 0, 0.2)
    g2r2 <- simulateSummarization(2, 2, 25, c(0, 1), 0.2)
    g2r5 <- simulateSummarization(2, 5, 25, c(0, 1), 0.2)
    g4r2 <- simulateSummarization(4, 2, 25, 1:4, 0.2)
    g4r5 <- simulateSummarization(4, 5, 25, 1:4, 0.2)

    expect_identical(ncol(g1r2), 3L)
    expect_identical(ncol(g1r5), 3L)
    expect_identical(ncol(g2r2), 3L)
    expect_identical(ncol(g2r5), 3L)
    expect_identical(ncol(g4r2), 3L)
    expect_identical(ncol(g4r5), 3L)

    expect_identical(nrow(g1r2), 2L)
    expect_identical(nrow(g1r5), 5L)
    expect_identical(nrow(g2r2), 4L)
    expect_identical(nrow(g2r5), 10L)
    expect_identical(nrow(g4r2), 8L)
    expect_identical(nrow(g4r5), 20L)
})

test_that("PTMsimulateExperiment() returns the right format", {
    s <- PTMsimulateExperiment(
        nGroup=2, nRep=2, nProtein=1, nSite=2, nFeature=5,
        logAbundance=list(
            PTM=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05),
            PROTEIN=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05)
        )
    )
    # s <- PTMsimulateExperiment(
    #     nGroup = 2, nRep = 2, nProtein = 1, nSite = 2, nFeature = 5,
    #     mu = list(PTM = 25, PROTEIN = 25),
    #     delta = list(PTM = c(0, 1), PROTEIN = c(0, 1)),
    #     sRep = list(PTM = 0.2, PROTEIN = 0.2),
    #     sPeak = list(PTM = 0.05, PROTEIN = 0.05)
    # )

    expect_output(str(s), "List of 2")
    expect_identical(names(s), c("PTM", "PROTEIN"))
    expect_is(s[["PTM"]], "data.frame")
    expect_is(s[["PROTEIN"]], "data.frame")
})
