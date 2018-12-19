#' Read and tidy a FASTA file.
#'
#' \code{tidy_fasta} reads and tidys FASTA file.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate
#' @importFrom stringr str_extract str_remove str_remove_all
#' @importFrom tibble as_data_frame rownames_to_column
#' @param filepath A single string that contains the path to read the FASTA file.
#' @return A tibble with columns named \code{header}, \code{sequence},
#'   \code{uniprot_ac}, \code{uniprot_iso}, \code{entry_name}.
#' @export
#'
#' @examples
#' tidy_fasta(filepath)
tidy_fasta <- function(filepath) {
    hs_fasta <- Biostrings::readAAStringSet(filepath)
    hs_fasta <- as.data.frame(hs_fasta) %>%
        as_data_frame() %>%
        rownames_to_column(var = "header") %>%
        select(header, sequence = x)

    regex_uniprot_ac <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")
    regex_uniprot_iso <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})([-]\\d{1,}){0,1}")

    hs_fasta <- hs_fasta %>%
        mutate(defline = str_extract(header, pattern = "([^\\s]*)(?=\\s)")) %>%
        mutate(
            uniprot_ac = str_extract(defline, pattern = regex_uniprot_ac),
            uniprot_iso = str_extract(defline, pattern = regex_uniprot_iso)
        ) %>%
        mutate(
            entry_name = str_remove(defline, "^(sp|tr)(?=\\|)") %>%
                str_remove(regex_uniprot_iso) %>%
                str_remove_all("\\|")
        ) %>%
        select(-defline)

    return(hs_fasta)
}
