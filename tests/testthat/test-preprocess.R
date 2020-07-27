test_that("PTMnormalize() takes and returns the right format", {
    s <- PTMsimulateExperiment(
        nGroup = 2, nRep = 2, nProtein = 1, nSite = 2, nFeature = 5,
        mu = list(PTM = 25, Protein = 25),
        delta = list(PTM = c(0, 1), Protein = c(0, 1)),
        sRep = list(PTM = 0.2, Protein = 0.2),
        sPeak = list(PTM = 0.05, Protein = 0.05)
    )

    expect_error(PTMnormalize(s[["PTM"]]), "PTM peak list is missing!")
    expect_output(str(PTMnormalize(s)), "List of 2")
    expect_identical(names(PTMnormalize(s)), c("PTM", "Protein"))

    s[["Protein"]] <- NULL
    expect_output(str(PTMnormalize(s)), "List of 1")
    expect_identical(names(PTMnormalize(s)), "PTM")
})
