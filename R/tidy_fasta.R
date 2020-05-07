#' Read and tidy a FASTA file.
#'
#' \code{tidy_fasta} reads and tidys FASTA file.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate
#' @importFrom stringr str_extract str_remove str_remove_all
#' @importFrom tibble tibble
#' @param aa_file A string of path to a FASTA file.
#' @return A tibble with columns named \code{header}, \code{sequence},
#'   \code{uniprot_ac}, \code{uniprot_iso}, \code{entry_name}.
#' @export
#'
#' @examples
#' tidy_fasta(aa_file)
tidy_fasta <- function(aa_file) {

    ## Check input
    if (missing(aa_file))
        stop(paste0("The input ", sQuote("aa_file"), " is missing"))
    if (!is.character(aa_file))
        stop("Provide the path to the FASTA file as a string")
    if (length(aa_file) != 1)
        stop("Provide only one path to the FASTA file at a time")
    # if (!file.exists(aa_file))
    #     stop(paste0("The file ", sQuote(aa_file), " does not exist"))

    aa <- Biostrings::readAAStringSet(aa_file)
    aa_vec <- as.character(aa)  # named vector
    aa_header <- names(aa_vec)
    aa_seq <- unname(aa_vec)

    pattern_ac <- regex(paste0("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]",
                               "([A-Z][A-Z0-9]{2}[0-9]){1,2})"))
    pattern_iso <- regex(paste0("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]",
                                "([A-Z][A-Z0-9]{2}[0-9]){1,2})([-]\\d{1,}){0,1}"))

    aa_sub <- str_extract(aa_header, pattern = "([^\\s]*)(?=\\s)")
    aa_ac <- str_extract(aa_sub, pattern = pattern_ac)
    aa_iso <- str_extract(aa_sub, pattern = pattern_iso)
    aa_entry <- str_remove(aa_sub, "^(sp|tr)(?=\\|)") %>%
        str_remove(pattern_iso) %>%
        str_remove_all("\\|")

    tibble(header = aa_header, sequence = aa_seq, uniprot_ac = aa_ac,
           uniprot_iso = aa_iso, entry_name = aa_entry)

}


#' Read and tidy a FASTA file.
#'
#' \code{tidy_fasta} reads and tidys FASTA file.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate
#' @importFrom stringr str_extract str_remove str_remove_all
#' @importFrom tibble enframe
#' @param aa_file A string of path to a FASTA file.
#' @return A tibble with columns named \code{header}, \code{sequence},
#'   \code{uniprot_ac}, \code{uniprot_iso}, \code{entry_name}.
#' @export
#'
#' @examples
#' tidy_fasta(aa_file)
tidy_fasta2 <- function(aa_file) {

    if (missing(aa_file))
        stop("Input `aa_file` is missing")
    if (!is.character(aa_file))
        stop("Provide the path to the FASTA file as a string")
    if (length(aa_file) != 1)
        stop("Provide only one path to the FASTA file at a time")

    aa_set <- Biostrings::readAAStringSet(aa_file)
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
