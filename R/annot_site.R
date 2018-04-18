
#' Annotate modification site.
#'
#' \code{annot_site} annotates modified sites as their residues and locations.
#'
#' @param aa_idx An integer vector. Location of the sites.
#' @param residue A string vector. Amino acid residue.
#' @param len_idx An integer. Default is \code{NULL}
#' @return A string.
#'
#' @examples
#' annot_site(10, "K")
#' annot_site(10, "K", 3)
#'
#' @export
#'
annot_site <- function(aa_idx, residue, len_idx = NULL) {
    if (purrr::is_empty(aa_idx)) {
        site <- "None"
    } else {
        if (length(residue) != 1 && length(aa_idx) != length(residue))
            stop("Lengths of aa_idx and residue don't match!")

        if (!is.null(len_idx)) {
            if (!is.integer(len_idx))
                stop("len_idx should be an integer!")
            if (max(stringr::str_length(aa_idx)) > len_idx)
                stop("len_idx doesn't reserve enough space for aa_idx!")

            aa_idx_pad <- stringr::str_pad(aa_idx, width = len_idx, pad = "0")
            site <- stringr::str_c(residue, aa_idx_pad, collapse = "-")
        } else {
            site <- stringr::str_c(residue, aa_idx, collapse = "-")
        }
    }

    return(site)
}
