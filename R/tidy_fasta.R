#' Read and tidy a FASTA file.
#'
#' \code{tidy_fasta} reads and tidys FASTA file.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate
#' @importFrom stringr str_extract str_remove str_remove_all
#' @importFrom tibble enframe
#' @param path_fasta A single string that contains the path to read the FASTA file.
#' @return A tibble with columns named \code{header}, \code{sequence},
#'   \code{uniprot_ac}, \code{uniprot_iso}, \code{entry_name}.
#' @export
#'
#' @examples
#' tidy_fasta(path_fasta)
tidy_fasta <- function(path_fasta) {

    if (missing(path_fasta))
        stop("Input path_fasta is missing")
    if (!is.character(path_fasta))
        stop("Provide the path to the FASTA file in string")
    if (length(path_fasta) != 1)
        stop("Provide only one path to the FASTA file at a time")

    aa_set <- Biostrings::readAAStringSet(path_fasta)
    aa_vec <- as.character(aa_set)  # named vector
    tbl_fasta <- enframe(aa_vec, name = "header", value = "sequence")

    uniprot_ac <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")
    uniprot_iso <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})([-]\\d{1,}){0,1}")

    tbl_fasta <- tbl_fasta %>%
        mutate(defline = str_extract(header, pattern = "([^\\s]*)(?=\\s)")) %>%
        mutate(
            uniprot_ac = str_extract(defline, pattern = uniprot_ac),
            uniprot_iso = str_extract(defline, pattern = uniprot_iso)
        ) %>%
        mutate(
            entry_name = str_remove(defline, "^(sp|tr)(?=\\|)") %>%
                str_remove(uniprot_iso) %>%
                str_remove_all("\\|")
        ) %>%
        select(-defline)

    return(tbl_fasta)
}
