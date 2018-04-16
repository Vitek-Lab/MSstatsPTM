
#' Locate modified sites.
#'
#' \code{locate_mod} locates modified sites on a peptide.
#'
#' @param pep_seq A string. Peptide sequence.
#' @param pep_idx A matrix. Starting and ending indices of the peptide.
#' @param mod_res_symbol A string. Modification residue and denoted symbol.
#' @return A string.
#'
#' @examples
#' locate_mod(pep_seq, pep_idx, mod_res_symbol)
#'
#' @export
#'
locate_mod <- function(pep_seq, pep_idx, mod_res_symbol) {
    pep_start <- unname(pep_idx[, "start"])
    site_relidx <- stringr::str_locate_all(pep_seq, mod_res_symbol)[[1]]
    site_relstart <- unname(site_relidx[, "start"])
    mod_idx <- site_relstart - seq_along(site_relstart) + pep_start

    return(mod_idx)
}
