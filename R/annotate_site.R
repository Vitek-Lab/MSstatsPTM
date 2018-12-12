#' Annotate modification site.
#'
#' \code{annotate_site} annotates modified sites as their residues and locations.
#'
#' @param aa_idx An integer vector. Location of the sites.
#' @param aa_idx_all An integer vector. Location of all sites of the protein.
#' @param residue A string vector. Amino acid residue.
#' @return A string.
#' @export
#'
#' @examples
#' annotate_site(10, c(10, 100), "K")
#' annotate_site(c(10, 13), c(10, 13, 100), c("K", "K"))
annotate_site <- function(aa_idx, aa_idx_all, residue) {
    if (purrr::is_empty(aa_idx)) {
        site <- "None"
    } else {
        if (length(residue) != 1 && length(aa_idx) != length(residue))
            stop("Lengths of aa_idx and residue don't match!")

        idx_len <- max(stringr::str_length(aa_idx_all))
        aa_idx_pad <- stringr::str_pad(aa_idx, width = idx_len, pad = "0")
        site <- stringr::str_c(residue, aa_idx_pad, collapse = "-")
    }

    return(site)
}
