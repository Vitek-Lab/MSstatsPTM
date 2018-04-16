
#' Annotate modification site.
#'
#' \code{annotate_site} annotates modified sites as their residues and locations.
#'
#' @param aa_idx An integer vector. Locations of the sites.
#' @param index_full An integer vector. Locations of all sites of the protein.
#' @param residue A string vector.
#' @return A string.
#'
#' @examples
#' annotate_site(10, c(10, 100), "K")
#' annotate_site(c(10, 13), c(10, 13, 100), c("K", "K"))
#'
#' @export
#'
annotate_site <- function(aa_idx, index_full, residue) {
    if (purrr::is_empty(aa_idx)) {
        site_str <- "None"
    } else {
        ant_len <- max(stringr::str_length(index_full))
        aa_idx_pad <- stringr::str_pad(aa_idx, width = ant_len, pad = "0")
        site_str <- stringr::str_c(residue, aa_idx_pad, collapse = "-")
    }

    return(site_str)
}
