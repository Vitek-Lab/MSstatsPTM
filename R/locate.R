#' Read and tidy a FASTA file
#'
#' \code{tidyFasta} reads and tidys FASTA file.
#'
#' @param path A string of path to a FASTA file.
#'
#' @return A tibble with columns named \code{header}, \code{sequence},
#'   \code{uniprot_ac}, \code{uniprot_iso}, \code{entry_name}.
#'
#' @examples
#' tidyFasta("https://www.uniprot.org/uniprot/O13297.fasta")
#'
#' @export
tidyFasta <- function(path) {

    # Check input
    if (missing(path))
        stop(paste0("The input ", sQuote("path"), " is missing"))
    if (!is.character(path))
        stop("Provide the path to the FASTA file as a string")
    if (length(path) != 1)
        stop("Provide only one path to the FASTA file at a time")
    # if (!file.exists(path))
    #     stop(paste0("The file ", sQuote(path), " does not exist"))

    aa <- Biostrings::readAAStringSet(path)
    aa <- as.character(aa)  # named vector
    header <- names(aa)

    pat_ac <- paste0(
        "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]",
        "([A-Z][A-Z0-9]{2}[0-9]){1,2})"
    )
    pat_iso <- paste0(
        "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]",
        "([A-Z][A-Z0-9]{2}[0-9]){1,2})([-]\\d{1,}){0,1}"
    )

    mch_head <- regexpr(
        pattern = "([^\\s]*)(?=\\s)", text = header, perl = TRUE
    )
    shead <- regmatches(header, m = mch_head)

    mch_uniprot <- regexpr(pattern = pat_ac, text = shead)
    matched <- mch_uniprot != (-1)
    ac <- regmatches(shead, regexpr(pattern = pat_ac, text = shead))
    iso <- regmatches(shead, regexpr(pattern = pat_iso, text = shead))

    trimmed <- sub("^(sp|tr)(?=\\|)", "", shead, perl = TRUE)
    trimmed <- sub(pat_iso, "", trimmed)
    entry <- gsub("\\|", "", trimmed)

    tibble(
        header = header[matched], sequence = as.vector(aa[matched]),
        uniprot_ac = ac, uniprot_iso = iso, entry_name = entry[matched]
    )
}


#' Annotate modified sites with associated peptides
#'
#' \code{PTMlocate} annotates modified sites with associated peptides.
#'
#' @param peptide A string vector of peptide sequences. The peptide sequence
#'   does not include its preceding and following AAs.
#' @param uniprot A string vector of Uniprot identifiers of the peptides'
#'   originating proteins. UniProtKB entry isoform sequence is used.
#' @param fasta A tibble with FASTA information. Output of \code{tidyFasta}.
#' @param modResidue A string. Modifiable amino acid residues.
#' @param modSymbol A string. Symbol of a modified site.
#' @param rmConfound A logical. \code{TRUE} removes confounded
#'   unmodified sites, \code{FALSE} otherwise. Default is \code{FALSE}.
#'
#' @return A data frame with three columns: \code{uniprot_iso}, \code{peptide},
#'   \code{site}.
#'
#' @examples
#' fasta <- tidyFasta("https://www.uniprot.org/uniprot/O13297.fasta")
#' PTMlocate("DRVSYIHNDSC*TR", "O13297", fasta, "C", "\\*")
#'
#' @export
PTMlocate <- function(peptide, uniprot, fasta, modResidue, modSymbol,
    rmConfound=FALSE) {

    good <- .locateCheck(peptide, uniprot, fasta, modResidue, modSymbol)

    peptide_seq <- tibble(uniprot_iso = uniprot, peptide = peptide)
    peptide_seq$is_mod <- grepl(modSymbol, peptide_seq$peptide)
    peptide_seq$peptide_unmod <- gsub(modSymbol, "", peptide_seq$peptide)
    # Locate modifiable sites
    col_seq <- c("uniprot_iso", "sequence")
    sub_fasta <- fasta[fasta$uniprot_iso %in% uniprot, col_seq]
    loc <- gregexpr(modResidue, sub_fasta$sequence, fixed = TRUE)
    sub_fasta$idx_site_full <- lapply(loc, .location_start)
    # Locate peptides (use extended AAs for specific matching when possible)
    peptide_fasta <- left_join(peptide_seq, sub_fasta)
    loc <- Map(function(p, s) gregexpr(p, s, fixed = TRUE)[[1]],
        as.list(peptide_fasta$peptide_unmod), as.list(peptide_fasta$sequence))
    n <- vapply(loc, .num_match, FUN.VALUE = integer(1))
    peptide_fasta <- peptide_fasta[n == 1L, ]
    loc <- loc[n == 1L]
    peptide_fasta$idx_peptide <- lapply(loc, .location)
    peptide_fasta$aa_start <- vapply(loc, .location_start, numeric(1))
    # Pattern of modified sites
    mod_pattern <- paste0(modResidue, modSymbol)
    # Locate modifiable, modified sites (AA residues) associated with peptides
    peptide_fasta$idx_site <- Map(.covered_set,
        peptide_fasta$idx_site_full, peptide_fasta$idx_peptide)
    peptide_fasta$idx_mod <- Map(
        function(p, a) locateMod(p, a, residueSymbol = mod_pattern),
        as.list(peptide_fasta$peptide), as.list(peptide_fasta$aa_start)
    )
    peptide_fasta$mod_aa <- Map(
        function(s, i) if (length(i) == 0) character() else substring(s, i, i),
        as.list(peptide_fasta$sequence), peptide_fasta$idx_mod
    )
    # Annotate modified sites
    peptide_fasta$len_site <- lapply(
        peptide_fasta$idx_site_full, function(x) nchar(x[length(x)])
    )
    peptide_fasta$site <- Map(annotSite,
        peptide_fasta$idx_mod, peptide_fasta$mod_aa, peptide_fasta$len_site)
    peptide_fasta$site <- as.character(peptide_fasta$site)
    col_fasta <- c("uniprot_iso", "peptide", "peptide_unmod", "is_mod",
        "idx_site", "idx_mod", "site")
    peptide_fasta <- peptide_fasta[, col_fasta]

    .rmConfounded(peptide_fasta, rmConfound)
}

#' @keywords internal
.locateCheck <- function(peptide, uniprot, fasta, modResidue, modSymbol) {
    if (missing(peptide))
        stop("Input peptide is missing!")
    if (missing(uniprot))
        stop("Input uniprot is missing!")
    if (missing(fasta))
        stop("Input fasta is missing!")
    if (missing(modResidue))
        stop("Input modResidue is missing!")
    if (missing(modSymbol))
        stop("Input modSymbol is missing!")
    if (!is.character(peptide))
        stop("Please provide peptide sequence as character in peptide!")
    if (!is.character(uniprot))
        stop("Please provide Uniprot protein ID as character in uniprot!")
    if (length(peptide) != length(uniprot))
        stop("peptide and uniprot must be of the same length")
    if (!is.data.frame(fasta))
        stop("Please provide the FASTA information in a data frame!")
    if (!all(c("uniprot_iso", "sequence") %in% names(fasta)))
        stop("Uniprot_iso or sequence is missing from FASTA data frame!")
    TRUE
}

#' @keywords internal
.rmConfounded <- function(peptideFasta, rmConfound) {
    # Handle confounded unmodified sites
    col_res <- c("uniprot_iso", "peptide", "site")
    if (rmConfound) {
        by_prot <- group_by(peptideFasta, .data$uniprot_iso)
        by_prot <- mutate(
            by_prot, idx_mod_full = list(unique(unlist(.data$idx_mod)))
        )
        by_prot <- ungroup(by_prot)
        by_prot <- by_prot[!by_prot$is_mod, ]
        cnfnd <- vector("logical", nrow(by_prot))
        for (i in seq_along(cnfnd)) {
            cnfnds <- by_prot$idx_site[[i]] %in% by_prot$idx_mod_full[[i]]
            cnfnd[i] <- any(cnfnds)
        }
        s_unmod <- by_prot[!cnfnd, col_res]

        n1 <- vapply(peptideFasta$idx_mod, .length_one, FUN.VALUE = logical(1))
        s_mod <- peptideFasta[peptideFasta$is_mod & n1, col_res]
        res <- bind_rows(s_mod, s_unmod)
    } else {
        n1 <- vapply(peptideFasta$idx_mod, .length_one, FUN.VALUE = logical(1))
        res <- peptideFasta[!peptideFasta$is_mod | n1, col_res]
    }
    res
}

#' @keywords internal
.num_match <- function(x) {
    ifelse(
        is.na(x) || identical(as.vector(x), -1L),
        0L, length(attr(x, "match.length"))
    )
}

#' @keywords internal
.location_start <- function(x) {
    start <- as.vector(x)
    if (identical(start, -1L)) {
        return(integer())
    }
    start
}

#' @keywords internal
.location <- function(x) {
    start <- as.vector(x)
    if (identical(start, -1L)) {
        return(cbind(start = integer(), end = integer()))
    }
    end <- as.vector(x) + attr(x, "match.length") - 1

    cbind(start = start, end = end)
}

#' @keywords internal
.covered_set <- function(x, y) {
    x[x >= y[1] & x <= y[2]]
}

#' @keywords internal
.length_one <- function(x) length(x) == 1

#' Locate modified sites with a peptide
#'
#' \code{locateMod} locates modified sites with a peptide.
#'
#' @param peptide A string. Peptide sequence.
#' @param aaStart An integer. Starting index of the peptide.
#' @param residueSymbol A string. Modification residue and denoted symbol.
#'
#' @return A string.
#'
#' @examples
#' locateMod("P*EP*TIDE", 3, "\\*")
#'
#' @export
locateMod <- function(peptide, aaStart, residueSymbol) {
    x <- gregexpr(residueSymbol, peptide)[[1]]
    start <- as.vector(x)
    if (identical(start, -1L)) {
        return(integer())
    }
    start - seq_along(start) + aaStart
}

#' Annotate modification site
#'
#' \code{annotSite} annotates modified sites as their residues and locations.
#'
#' @param aaIndex An integer vector. Location of the sites.
#' @param residue A string vector. Amino acid residue.
#' @param lenIndex An integer. Default is \code{NULL}
#'
#' @return A string.
#'
#' @examples
#' annotSite(10, "K")
#' annotSite(10, "K", 3L)
#'
#' @export
annotSite <- function(aaIndex, residue, lenIndex=NULL) {
    if (length(aaIndex) == 0) {
        return("None")
    } else {
        if (length(residue) != 1 && length(aaIndex) != length(residue))
            stop("Lengths of aaIndex and residue don't match!")
        if (!is.null(lenIndex)) {
            if (!is.integer(lenIndex))
                stop("lenIndex should be an integer!")
            if (max(nchar(aaIndex)) > lenIndex)
                stop("lenIndex doesn't reserve enough space for aaIndex!")

            aa_idx_pad <- formatC(
                aaIndex, width = lenIndex, format = "d", flag = "0"
            )
            site <- paste0(residue, aa_idx_pad, collapse = "-")
        } else {
            site <- paste0(residue, aaIndex, collapse = "-")
        }
    }
    site
}

