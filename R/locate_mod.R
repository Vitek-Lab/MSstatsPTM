#' Locate modified sites with a peptide.
#'
#' \code{locate_mod} locates modified sites with a peptide.
#'
#' @importFrom stringr str_locate_all
#' @param peptide A string. Peptide sequence.
#' @param aa_start An integer. Starting index of the peptide.
#' @param residue_symbol A string. Modification residue and denoted symbol.
#' @return A string.
#' @export
#'
#' @examples
#' locate_mod(peptide, aa_start, residue_symbol)
locate_mod <- function(peptide, aa_start, residue_symbol) {
    mod_rels <- str_locate_all(peptide, residue_symbol)[[1]]
    mod_rel_idx <- unname(mod_rels[, "start"])
    mod_idx <- mod_rel_idx - seq_along(mod_rel_idx) + aa_start

    return(mod_idx)
}
