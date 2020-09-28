#' Simulate PTM quantification experiments
#'
#' \code{PTMsimulateExperiment} simulates a PTM quantification experiment with
#' a list of log2-intensities of multiple spectral features, PTM sites and the
#' corresponding proteins, in multiple MS runs and conditions.
#'
#' @param nGroup An integer to specify the number of conditions.
#' @param nRep An integer to specify the number of replicates per condition.
#' @param nProtein An integer to specify the number of protein.
#' @param nSite An integer to specify the number of PTM sites per protein.
#' @param nFeature An integer to specify the number of features per site.
#' @param logAbundance A list of two lists named \code{PTM} and \code{PROTEIN}.
#'   Each contains four elements: \code{mu} (a numeric representing the overall
#'   mean log2-abundance), \code{delta} (a numeric vector for the deviation of
#'   the mean log2-abundance for each group from the overall mean), \code{sRep}
#'   (a numeric representing the standard deviation for run-to-run variation),
#'   and \code{sPeak} (a numeric representing the standard deviation in peak
#'   log2-intensities).
#'
#' @return A tibble with columns of \code{protein}, \code{site}, \code{group},
#'   \code{run}, \code{feature}, \code{log2inty}.
#'
#' @examples
#' PTMsimulateExperiment(
#'     nGroup=2, nRep=2, nProtein=1, nSite=1, nFeature=5,
#'     logAbundance=list(
#'         PTM=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05),
#'         PROTEIN=list(mu=25, delta=c(0, 1), sRep=0.2, sPeak=0.05)
#'     )
#' )
#'
#' @export
PTMsimulateExperiment <- function(nGroup, nRep, nProtein, nSite, nFeature,
    logAbundance) {

    # PTM data
    peaks <- vector("list", nProtein)
    for (i in seq_len(nProtein)) {
        peaks[[i]] <- simulateSites(
            nGroup, nRep, nSite, nFeature,
            logAbundance$PTM$mu, logAbundance$PTM$delta,
            logAbundance$PTM$sRep, logAbundance$PTM$sPeak
        )
    }
    sites <- bind_rows(peaks)
    sites$protein <- rep(
        paste0("Protein_", seq_len(nProtein)),
        each = nGroup * nRep * nSite * nFeature
    )

    # PROTEIN data
    peaks <- vector("list", nProtein)
    for (i in seq_len(nProtein)) {
        peaks[[i]] <- simulatePeaks(
            nGroup, nRep, nFeature,
            logAbundance$PROTEIN$mu, logAbundance$PROTEIN$delta,
            logAbundance$PROTEIN$sRep, logAbundance$PROTEIN$sPeak
        )
    }
    prots <- bind_rows(peaks)
    prots$protein <- rep(
        paste0("Protein_", seq_len(nProtein)), each = nGroup * nRep * nFeature
    )

    # Combine PTM and PROTEIN data
    cols <- c("protein", "site", "group", "run", "feature", "log2inty")
    list(PTM = sites[, cols], PROTEIN = prots[, setdiff(cols, "site")])
}

#' Simulate peak log-intensities for PTM sites
#'
#' \code{simulateSites} simulates a list of log2-intensities of multiple
#' spectral features and PTM sites of one protein, in multiple MS runs and
#' conditions.
#'
#' @param nGroup An integer to specify the number of conditions.
#' @param nRep An integer to specify the number of replicates per condition.
#' @param nSite An integer to specify the number of PTM sites per protein.
#' @param nFeature An integer to specify the number of features per site.
#' @param mu A numeric to specify the overall mean log2-intensity.
#' @param delta A numeric to specify the deviation of the mean log2-abundance of
#'   each group from the overall mean.
#' @param sRep A numeric to specify the standard deviation for run-to-run
#'   variation.
#' @param sPeak A numeric to specify the standard deviation in peak
#'   log2-intensities.
#'
#' @return A tibble with columns of \code{site}, \code{group}, \code{run},
#'   \code{feature}, \code{log2inty}.
#'
#' @examples
#' simulateSites(nGroup=2, nRep=2, nSite=2, nFeature=5, 25, c(0, 1), 0.2, 0.05)
#'
#' @export
simulateSites <- function(nGroup, nRep, nSite, nFeature, mu, delta, sRep,
    sPeak) {

    peaks <- vector("list", nSite)
    for (i in seq_len(nSite)) {
        peaks[[i]] <- simulatePeaks(
            nGroup, nRep, nFeature, mu, delta, sRep, sPeak
        )
    }
    sites <- bind_rows(peaks)
    sites$site <- rep(
        paste0("S_", seq_len(nSite)), each = nGroup * nRep * nFeature
    )
    sites
}


#' Simulate peak log2-intensities
#'
#' \code{simulateSites} simulates a list of log2-intensities of multiple
#' spectral features of a PTM site, in multiple MS runs and conditions.
#'
#' @param nGroup An integer to specify the number of conditions.
#' @param nRep An integer to specify the number of replicates per condition.
#' @param nFeature An integer to specify the number of features per site.
#' @param mu A numeric to specify the overall mean log2-intensity.
#' @param delta A numeric to specify the deviation of the mean log2-abundance of
#'   each group from the overall mean.
#' @param sRep A numeric to specify the standard deviation for run-to-run
#'   variation.
#' @param sPeak A numeric to specify the standard deviation in peak
#'   log2-intensities.
#'
#' @return A tibble with columns of \code{group}, \code{run}, \code{feature},
#'   and \code{log2inty}.
#'
#' @examples
#' simulatePeaks(nGroup=2, nRep=3, nFeature=5, 25, c(0, 1), 0.2, 0.05)
#'
#' @export
simulatePeaks <- function(nGroup, nRep, nFeature, mu, delta, sRep, sPeak) {
    summarized <- simulateSummarization(nGroup, nRep, mu, delta, sRep)

    tibble(
        group = rep(summarized$group, each = nFeature),
        run = rep(summarized$run, each = nFeature),
        feature = rep(paste0("F_", seq_len(nFeature)), nGroup * nRep),
        log2inty = unlist(Map(
            stats::rnorm, nFeature, summarized$log2inty, sPeak
        ))
    )
}


#' Simulate site-level summarization for PTM experiment
#'
#' \code{simulateSummarization} simulates the summarized log2-intensity value
#' of a PTM site in each MS run. The value is randomly generated based on a
#' normal distribtuion, where the average log2-intensity is defined for each
#' condition
#'
#' @param nGroup An integer to specify the number of conditions.
#' @param nRep An integer to specify the number of replicates per condition.
#' @param mu A numeric value of the overall mean log2-abundance.
#' @param delta A numeric vector to specify the deviation of the mean
#'   log2-abundance of each group from the overall mean.
#' @param sRep A numeric. Standard deviation of the log2-intensities.
#'
#' @return A tibble with columns of \code{group}, \code{run} and
#'   \code{log2inty}.
#'
#' @examples
#' simulateSummarization(nGroup=2, nRep=3, 25, c(0, 1), 0.2)
#'
#' @export
simulateSummarization <- function(nGroup, nRep, mu, delta, sRep) {
    all_mu <- rep(mu, nGroup) + delta

    tibble(
        group = rep(paste0("G_", seq_len(nGroup)), each = nRep),
        run = paste0("R_", seq_len(nGroup * nRep)),
        log2inty = unlist(Map(stats::rnorm, nRep, all_mu, sRep))
    )
}
