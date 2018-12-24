#' Annotate modified sites with associated peptides.
#'
#' \code{peptide2site} annotate modified sites with associated peptides.
#'
#' @importFrom tibble tibble
#' @importFrom dplyr filter select mutate group_by ungroup left_join bind_rows
#' @importFrom stringr str_detect str_remove_all str_locate_all str_count str_c
#'   str_sub str_length
#' @importFrom purrr map map2 pmap map_int map2_lgl pmap_chr
#' @param peptide A string vector of peptide sequences. The sequences do not
#'   include the preceding and following AAs.
#' @param uniprot A string vector of Uniprot identifiers of the peptides'
#'   originating proteins. UniProtKB entry isoform sequence is used.
#' @param fasta A tibble with FASTA information. Output of \code{tidy_fasta}.
#' @param mod_residue A string. Modifiable amino acid residues.
#' @param mod_symbol A string. Symbol of a modified site.
#' @param remove_confounding A logical. \code{TRUE} removes confounded
#'   unmodified sites, \code{FALSE} otherwise. Default is \code{FALSE}.
#' @return A data frame of three columns: \code{uniprot_iso}, \code{peptide},
#'   \code{site}.
#' @export
#'
peptide2site <- function(peptide, uniprot, fasta, mod_residue, mod_symbol,
                         remove_confounding = FALSE) {

    if (length(peptide) != length(uniprot))
        stop("Length of peptide and uniprot don't match!")
    if (!all(c("uniprot_iso", "sequence") %in% names(fasta)))
        stop("Column of uniprot_iso or sequence is missing from fasta!")

    peptide_seq <- tibble(uniprot_iso = uniprot, peptide = peptide) %>%
        mutate(is_mod = str_detect(peptide, mod_symbol)) %>%
        mutate(peptide_unmod = str_remove_all(peptide, mod_symbol))

    # Locate modifiable sites
    sub_fasta <- fasta %>%
        filter(uniprot_iso %in% peptide_seq$uniprot_iso) %>%
        select(uniprot_iso, sequence) %>%
        mutate(idx_site_full = str_locate_all(sequence, mod_residue) %>%
                   map(~unname(.[, "start"])))

    # Locate peptides
    # [TODO]: use extended AAs for more specific matching
    peptide_fasta <- peptide_seq %>%
        left_join(sub_fasta) %>%
        filter(str_count(sequence, peptide_unmod) == 1) %>%
        mutate(idx_peptide = str_locate_all(sequence, peptide_unmod)) %>%
        mutate(aa_start = map_int(idx_peptide, ~.[, "start"]))

    # Pattern of modified sites
    mod_pattern <- str_c(mod_residue, mod_symbol)

    # Locate modifiable, modified sites (AA residues) associated with peptides
    peptide_fasta <- peptide_fasta %>%
        mutate(
            idx_site = map2(idx_site_full, idx_peptide, ~.x[.x >= .y[1] & .x <= .y[2]]),
            idx_mod = pmap(list(peptide, aa_start, mod_pattern), locate_mod),
            mod_aa = pmap(list(sequence, idx_mod, idx_mod), str_sub)
        )

    # Annotate modified sites
    peptide_fasta <- peptide_fasta %>%
        mutate(
            len_site = map(idx_site_full, ~ str_length(.[length(.)])),
            site = pmap_chr(list(idx_mod, mod_aa, len_site), annot_site)
        ) %>%
        select(uniprot_iso, peptide, peptide_unmod, is_mod, idx_site, idx_mod, site)

    # Handle confounded unmodified sites
    if (remove_confounding) {
        uncfd_unmod <- peptide_fasta %>%
            group_by(uniprot_iso) %>%
            mutate(idx_mod_full = list(unique(unlist(idx_mod)))) %>%
            ungroup() %>%
            filter(!is_mod) %>%
            filter(map2_lgl(idx_site, idx_mod_full, ~ !any(.x %in% .y))) %>%
            select(uniprot_iso, peptide, site)

        res <- peptide_fasta %>%
            filter(is_mod) %>%
            filter(map_int(idx_mod, length) == 1) %>%
            select(uniprot_iso, peptide, site) %>%
            bind_rows(uncfd_unmod)
    } else {
        res <- peptide_fasta %>%
            filter(!is_mod | map_int(idx_mod, length) == 1) %>%
            select(uniprot_iso, peptide, site)
    }

    return(res)
}
