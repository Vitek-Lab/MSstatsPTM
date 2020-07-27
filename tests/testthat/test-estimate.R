test_that("PTMestimate() takes and returns the right format", {
    s <- PTMsimulateExperiment(
        nGroup = 2, nRep = 2, nProtein = 1, nSite = 2, nFeature = 5,
        mu = list(PTM = 25, Protein = 25),
        delta = list(PTM = c(0, 1), Protein = c(0, 1)),
        sRep = list(PTM = 0.2, Protein = 0.2),
        sPeak = list(PTM = 0.05, Protein = 0.05)
    )
    est <- PTMestimate(PTMnormalize(s))

    expect_output(str(est), "List of 2")
    expect_identical(names(est), c("PTM", "Protein"))

    expect_output(str(est[["PTM"]]), "List of 4")
    expect_output(str(est[["Protein"]]), "List of 3")

    expect_identical(names(est[["PTM"]]), c("protein", "site", "param", "df"))
    expect_identical(names(est[["Protein"]]), c("protein", "param", "df"))

    s[["Protein"]] <- NULL
    est <- PTMestimate(PTMnormalize(s))

    expect_output(str(est), "List of 1")
    expect_identical(names(est), "PTM")
})
