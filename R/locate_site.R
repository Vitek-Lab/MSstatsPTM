
#' Locate potential modification sites with a peptide.
#'
#' \code{locate_site} locates modifiable sites with a peptide.
#'
#' @param peptide A string. Peptide sequence.
#' @param pep_idx A matrix. Starting and ending indices of the peptide.
#' @param residue A string. Amino acid residue for the modification.
#' @return A string.
#'
#' @examples
#' idx <- matrix(c(2, 5), nrow = 1, dimnames = list(NULL, c("start", "end")))
#' locate_site("ABCA", idx, "A")
#' locate_site("ABCA", idx, "A|C")
#'
#' @export
#'
locate_site <- function(peptide, pep_idx, residue) {
    aa_start <- unname(pep_idx[, "start"])
    aa_rels <- stringr::str_locate_all(peptide, residue)[[1]]
    aa_rel_idx <- unname(aa_rels[, "start"])
    aa_idx <- aa_rel_idx + aa_start - 1

    return(aa_idx)
}
