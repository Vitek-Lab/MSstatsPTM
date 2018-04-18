
#' Read and tidy FASTA file.
#'
#' \code{tidy_fasta} reads and tidys FASTA file.
#'
#' @param path A string. Path to the FASTA file.
#' @return A data frame.
#'
#' @examples
#' tidy_fasta(path)
#'
#' @export
#'
tidy_fasta <- function(path) {
    hs_fasta <- Biostrings::readAAStringSet(path)
    hs_fasta <- as.data.frame(hs_fasta) %>%
        dplyr::tbl_df() %>%
        tibble::rownames_to_column(var = "header") %>%
        dplyr::select(header, sequence = x)

    regex_uniprot_ac <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")
    regex_uniprot_iso <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})([-]\\d{1,}){0,1}")

    hs_fasta <- hs_fasta %>%
        dplyr::mutate(
            uniprot_ac = stringr::str_extract(header, pattern = regex_uniprot_ac),
            uniprot_iso = stringr::str_extract(header, pattern = regex_uniprot_iso),
            entry_name = stringr::str_extract(header, pattern = "([^\\s\\|]*)(?=\\s)")
        )

    return(hs_fasta)
}
