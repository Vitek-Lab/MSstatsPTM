test_that("PTMnormalize() takes and returns the right format", {
    s <- PTMsimulateExperiment(
        nGroup=2, nRep=2, nProtein=1, nSite=2, nFeature=5,
        logAbundance=list(
            PTM=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05),
            PROTEIN=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05)
        )
    )

    expect_error(PTMnormalize(s[["PTM"]]), "PTM peak list is missing!")
    expect_output(str(PTMnormalize(s)), "List of 2")
    expect_identical(names(PTMnormalize(s)), c("PTM", "PROTEIN"))

    s[["PROTEIN"]] <- NULL
    expect_output(str(PTMnormalize(s)), "List of 1")
    expect_identical(names(PTMnormalize(s)), "PTM")
})
