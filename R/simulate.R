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
#' @param mu A list of two elements named \code{PTM} and \code{Protein}. Each
#'   is a numeric representing the overall mean log2-intensity in the PTM or the
#'   Protein data.
#' @param delta A list of two elements named \code{PTM} and \code{Protein}. Each
#'   specifies the deviation of the mean log2-abundance of each group from the
#'   overall mean in the PTM or the Protein data.
#' @param sRep A list of two elements named \code{PTM} and \code{Protein}. Each
#'   is a numeric representing the standard deviation for run-to-run variation.
#' @param sPeak A list of two elements named \code{PTM} and \code{Protein}. Each
#'   is a numeric representing the standard deviation in peak log2-intensities.
#'
#' @return A tibble with columns of \code{protein}, \code{site}, \code{group},
#'   \code{run}, \code{feature}, \code{log2inty}.
#'
#' @export
#' @examples
#' PTMsimulateExperiment(nGroup = 2, nRep = 2, nProtein = 1, nSite = 2, nFeature = 5,
#' list(PTM = 25, Protein = 25), list(PTM = c(0, 1), Protein = c(0, 1)),
#' list(PTM = 0.2, Protein = 0.2), list(PTM = 0.05, Protein = 0.05))
PTMsimulateExperiment <- function(nGroup, nRep, nProtein, nSite, nFeature,
                                  mu, delta, sRep, sPeak) {
    # Site data
    peaks <- vector("list", nProtein)
    for (i in 1:nProtein) {
        peaks[[i]] <- simulateSites(nGroup, nRep, nSite, nFeature,
                                    mu$PTM, delta$PTM, sRep$PTM, sPeak$PTM)
    }
    sites <- bind_rows(peaks)
    sites$protein <- rep(paste0("Protein_", 1:nProtein),
                         each = nGroup * nRep * nSite * nFeature)

    # Protein data
    peaks <- vector("list", nProtein)
    for (i in 1:nProtein) {
        peaks[[i]] <- simulatePeaks(nGroup, nRep, nFeature, mu$Protein,
                                    delta$Protein, sRep$Protein, sPeak$Protein)
    }
    prots <- bind_rows(peaks)
    prots$protein <- rep(paste0("Protein_", 1:nProtein),
                         each = nGroup * nRep * nFeature)
    prots$site <- "None"

    # Combine PTM and Protein data
    cols <- c("protein", "site", "group", "run", "feature", "log2inty")
    bind_rows(sites[, cols], prots[, cols])
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
#' @export
#' @examples
#' simulateSites(nGroup = 2, nRep = 2, nSite = 2, nFeature = 5, 25, c(0, 1), 0.2, 0.05)
simulateSites <- function(nGroup, nRep, nSite, nFeature, mu, delta, sRep, sPeak) {
    peaks <- vector("list", nSite)
    for (i in 1:nSite) {
        peaks[[i]] <- simulatePeaks(nGroup, nRep, nFeature, mu, delta, sRep, sPeak)
    }
    sites <- bind_rows(peaks)
    sites$site <- rep(paste0("S_", 1:nSite), each = nGroup * nRep * nFeature)
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
#' @export
#' @examples
#' simulatePeaks(nGroup = 2, nRep = 3, nFeature = 5, 25, c(0, 1), 0.2, 0.05)
simulatePeaks <- function(nGroup, nRep, nFeature, mu, delta, sRep, sPeak) {
    summarized <- simulateSummarization(nGroup, nRep, mu, delta, sRep)

    tibble(
        group = rep(summarized$group, each = nFeature),
        run = rep(summarized$run, each = nFeature),
        feature = rep(paste0("F_", 1:nFeature), nGroup * nRep),
        log2inty = unlist(Map(stats::rnorm, nFeature, summarized$log2inty, sPeak))
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
#' @return A tibble with columns of \code{group}, \code{run} and \code{log2inty}.
#'
#' @export
#' @examples
#' simulateSummarization(nGroup = 2, nRep = 3, 25, c(0, 1), 0.2)
simulateSummarization <- function(nGroup, nRep, mu, delta, sRep) {
    all_mu <- rep(mu, nGroup) + delta

    tibble(
        group = rep(paste0("G_", 1:nGroup), each = nRep),
        run = paste0("R_", 1:(nGroup * nRep)),
        log2inty = unlist(Map(stats::rnorm, nRep, all_mu, sRep))
    )
}
