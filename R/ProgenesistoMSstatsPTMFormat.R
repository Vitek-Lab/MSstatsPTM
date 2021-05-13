#' Converts non-TMT Progenesis output into the format needed for MSstatsPTM
#'
#' @export
#' @importFrom MSstats ProgenesistoMSstatsFormat
#' @importFrom data.table as.data.table
#'
#' @param ptm_input name of Progenesis output with modified peptides, which is
#' wide-format. 'Accession', Sequence', 'Modification', 'Charge' and one column
#' for each run are required
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which
#' includes Condition, BioReplicate, Run information. It will be matched with
#' the column name of input for MS runs.
#' @param global_protein_input name of Progenesis output with unmodified
#' peptides, which is wide-format. 'Accession', Sequence', 'Modification',
#' 'Charge' and one column for each run are required
#' @param fasta_path string containing path to the corresponding fasta file for 
#' the modified peptide dataset.
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for
#' more than one proteins. We assume to use unique peptide for each protein.
#' @param summaryforMultipleRows max(default) or sum - when there are multiple
#' measurements for certain feature and certain run, use highest or sum of
#' multiple intensities.
#' @param fewMeasurements 'remove'(default) will remove the features that have 1
#'  or 2 measurements across runs.
#' @param removeOxidationMpeptides TRUE will remove the modified peptides
#' including 'Oxidation (M)' sequence. FALSE is default.
#' @param removeProtein_with1Peptide TRUE will remove the proteins which have
#' only 1 peptide and charge. FALSE is default.
#' @param mod.num For modified peptide dataset, must be one of `Single` or 
#' `Total`. The default is `Single`. The number modifications per peptide to be 
#' used. If "Single", only peptides with one modification will be 
#' used. Otherwise "Total" includes peptides with more than one modification.
#' Selecting "Total" may confound the effect of different modifications.
#' @examples
#' 
#' head(raw.input$PTM)
#' head(raw.input$PROTEIN)
#' 
ProgenesistoMSstatsPTMFormat <- function(ptm_input,
                                         annotation,
                                         global_protein_input = FALSE,
                                         fasta_path = FALSE,
                                         useUniquePeptide=TRUE,
                                         summaryforMultipleRows=max, 
                                         fewMeasurements="remove",
                                         removeOxidationMpeptides=FALSE,
                                         removeProtein_with1Peptide=FALSE,
                                         mod.num = 'Single'){

  ## TODO: Add logging
  ## TODO: Add checks

  ptm_input <- as.data.table(ptm_input)
  col_order <- colnames(ptm_input)
  ptm_input$id <- 1:nrow(ptm_input)
  
  X.10 = X.9 = X.8 = NULL
  
  ## Format PTM data
  if (fasta_path == FALSE) {
    ptm_input[,X.10 := paste(X.10, X.8, X.9, sep = "_")]
    ptm_input[2, "X.10"] <- "Accession"
  } else {
    .progensis.add.sites(ptm_input, fasta_path, col_order)
  }

  convert.ptm <- ProgenesistoMSstatsFormat(as.data.frame(ptm_input), annotation,
                                            useUniquePeptide,
                                            summaryforMultipleRows,
                                            fewMeasurements,
                                            removeOxidationMpeptides,
                                            removeProtein_with1Peptide)

  if (global_protein_input != FALSE){

    global_protein_input <- as.data.table(global_protein_input)

    convert.prot <- ProgenesistoMSstatsFormat(global_protein_input, annotation,
                                              useUniquePeptide,
                                              summaryforMultipleRows,
                                              fewMeasurements,
                                              removeOxidationMpeptides,
                                              removeProtein_with1Peptide)
    MSstatsPTM.data <- list("PTM" = convert.ptm,
                            "PROTEIN" = convert.prot)
  } else {
    MSstatsPTM.data <- list("PTM" = convert.ptm,
                            "PROTEIN" = NULL)
  }
  
  return(MSstatsPTM.data)
  
}