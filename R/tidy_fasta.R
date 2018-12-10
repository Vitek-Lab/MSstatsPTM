
#' Read and tidy a FASTA file.
#'
#' \code{tidy_fasta} reads and tidys FASTA file.
#'
#' @param filepath A single string that contains the path to read the FASTA file.
#'
#' @return A tibble with columns named \code{header}, \code{sequence},
#' \code{uniprot_ac}, \code{uniprot_iso}, \code{entry_name}.
#'
#' @examples
#' tidy_fasta(filepath)
#'
#' @export
#'
tidy_fasta <- function(filepath) {
    hs_fasta <- Biostrings::readAAStringSet(filepath)
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

    # hs_fasta <- hs_fasta %>%
    #     dplyr::mutate(
    #         uniprot_ac = stringr::str_extract(header, pattern = regex_uniprot_ac),
    #         uniprot_iso = stringr::str_extract(header, pattern = regex_uniprot_iso)
    #     ) %>%
    #     dplyr::mutate(
    #         entry_name = stringr::str_extract(header, pattern = "([^\\s]*)(?=\\s)") %>%
    #             stringr::str_replace("^[A-Za-z]+\\|", "") %>%
    #             stringr::str_replace(regex_uniprot_iso, "") %>%
    #             stringr::str_replace("\\|", "")
    #     )

    return(hs_fasta)
}
