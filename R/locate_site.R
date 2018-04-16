
#' Locate potential modification sites.
#'
#' \code{locate_site} locates modifiable sites on a peptide.
#'
#' @param pep_seq A string. Peptide sequence.
#' @param pep_idx A matrix. Starting and ending indices of the peptide.
#' @param mod_res A string. Modification residue.
#' @return A string.
#'
#' @examples
#' locate_site(pep_seq, pep_idx, mod_res)
#'
#' @export
#'
locate_site <- function(pep_seq, pep_idx, mod_res) {
    pep_start <- unname(pep_idx[, "start"])
    site_relidx <- stringr::str_locate_all(pep_seq, mod_res)[[1]]
    site_relstart <- unname(site_relidx[, "start"])
    site_idx <- site_relstart + pep_start - 1

    return(site_idx)
}
