
#' Locate potential modification sites with a peptide.
#'
#' \code{locate_site} locates modifiable sites with a peptide.
#'
#' @param peptide A string. Peptide sequence.
#' @param aa_start An integer. Starting index of the peptide.
#' @param residue A string. Amino acid residue for the modification.
#' @return A string.
#'
#' @examples
#' locate_site("ABCA", 2, "A")
#' locate_site("ABCA", 2, "A|C")
#'
#' @export
#'
locate_site <- function(peptide, aa_start, residue) {
    aa_rels <- stringr::str_locate_all(peptide, residue)[[1]]
    aa_rel_idx <- unname(aa_rels[, "start"])
    aa_idx <- aa_rel_idx + aa_start - 1

    return(aa_idx)
}
