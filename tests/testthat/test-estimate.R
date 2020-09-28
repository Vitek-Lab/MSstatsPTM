test_that("PTMestimate() takes and returns the right format", {
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
    est <- PTMestimate(PTMnormalize(s))

    expect_output(str(est), "List of 2")
    expect_identical(names(est), c("PTM", "PROTEIN"))

    expect_output(str(est[["PTM"]]), "List of 4")
    expect_output(str(est[["PROTEIN"]]), "List of 3")

    expect_identical(names(est[["PTM"]]), c("protein", "site", "param", "df"))
    expect_identical(names(est[["PROTEIN"]]), c("protein", "param", "df"))

    s[["PROTEIN"]] <- NULL
    est <- PTMestimate(PTMnormalize(s))

    expect_output(str(est), "List of 1")
    expect_identical(names(est), "PTM")
})
