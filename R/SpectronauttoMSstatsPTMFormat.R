#' Converts raw PTM MS data from Spectronautt into the format needed for
#' MSstatsPTM.
#'
#' Takes as as input both raw PTM and global protein outputs from Spectronaut.
#'
#' @export
#' @importFrom MSstats SpectronauttoMSstatsFormat
#' @importFrom data.table as.data.table
#' @importFrom stringr str_extract str_count str_locate_all
#' @importFrom checkmate assertChoice assertLogical assertNumeric
#'
#' @param PTM.data name of PTM Spectronaut output, which is long-format.
#' @param fasta A string of path to a FASTA file, used to match PTM peptides.
#' @param Protein.data name of Global Protein Spectronaut output, which is long-format.
#' @param annotation name of 'annotation.txt' data which includes Condition,
#' BioReplicate, Run. If annotation is already complete in Spectronaut, use
#' annotation=NULL (default). It will use the annotation information from input.
#' @param intensity 'PeakArea'(default) uses not normalized peak area.
#' 'NormalizedPeakArea' uses peak area normalized by Spectronaut
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that
#' have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will
#' be replaced with zero and will be considered as censored missing values for
#' imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. default is 0.01.
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for
#' more than one proteins. We assume to use unique peptide for each protein.
#' @param removeFewMeasurements TRUE (default) will remove the features that
#' have 1 or 2 measurements across runs.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have
#' only 1 feature, which is the combination of peptide, precursor charge,
#' fragment and charge. FALSE is default.
#' @param removeNonUniqueProteins TRUE will remove proteins that were not
#' uniquely identified. IE if the protein column contains multiple proteins
#' seperated by ";". TRUE is default
#' @param modificationLabel String of modification name. Default is 'Phospho'. 
#' @param removeiRT TRUE will remove proteins that contain iRT. True is default
#' @param summaryforMultipleRows max(default) or sum - when there are multiple
#' measurements for certain feature and certain run, use highest or sum of
#' multiple intensities.
#' @param which.Conditions list of conditions to format into MSstatsPTM format.
#' If "all" all conditions will be used. Default is "all".
#' @return a list of two data.tables named 'PTM' and 'PROTEIN' in the format 
#' required by MSstatsPTM.
#' @examples
#' 
#' # The output should be in the following format.
#' head(raw.input$PTM)
#' head(raw.input$PROTEIN)
#'
SpectronauttoMSstatsPTMFormat <- function(PTM.data,
                                          fasta,
                                          Protein.data = NULL,
                                          annotation = NULL,
                                          intensity = 'PeakArea',
                                          filter_with_Qvalue = TRUE,
                                          qvalue_cutoff = 0.01,
                                          useUniquePeptide = TRUE,
                                          removeFewMeasurements = TRUE,
                                          removeProtein_with1Feature = FALSE,
                                          removeNonUniqueProteins = TRUE,
                                          modificationLabel = "Phospho",
                                          removeiRT = TRUE,
                                          summaryforMultipleRows=max,
                                          which.Conditions = 'all'){
  
  ## TODO: Add logging
  ## Check variable input
  
  ## Ensure format of input data
  PTM.data <- as.data.table(PTM.data)
  if (!is.null(Protein.data)){
    Protein.data <- as.data.table(Protein.data)
  }
  
  ## Check and filter for available conditions
  if (which.Conditions != 'all') {

    PTM_conditions <- unique(PTM.data$R.Condition)
    if (!is.null(Protein.data)){
      Protein_conditions <- unique(Protein.data$R.Condition)
    } else {
      Protein_conditions <- PTM_conditions
    }
    if((length(setdiff(PTM_conditions, which.Conditions)
    ) == length(PTM_conditions)) |
    (length(setdiff(Protein_conditions, which.Conditions)
    ) == length(Protein_conditions))){
      msg = (paste("None of the conditions specified in which.Conditions are",
                   "available in one or both of LiP/TrP datasets. Please",
                   "ensure the conditions listed in which.Conditions appear",
                   "in the input datasets"))
      stop(msg)
    }

    PTM_conditions <- PTM_conditions[(R.Condition %in% which.Conditions)]
    if (!is.null(Protein_conditions)){
      Protein_conditions <- Protein_conditions[
        (R.Condition %in% which.Conditions)]
    }
  }
  
  ## MSstats process
  df.ptm <- SpectronauttoMSstatsFormat(PTM.data, annotation, intensity,
                                       filter_with_Qvalue, qvalue_cutoff,
                                       useUniquePeptide, removeFewMeasurements,
                                       removeProtein_with1Feature,
                                       summaryforMultipleRows)
  df.ptm <- as.data.table(as.matrix(df.ptm))
  if (!is.null(Protein.data)){
    df.protein <- SpectronauttoMSstatsFormat(as.data.frame(Protein.data), 
                                         annotation, intensity,
                                         filter_with_Qvalue, qvalue_cutoff,
                                         useUniquePeptide, removeFewMeasurements,
                                         removeProtein_with1Feature,
                                         summaryforMultipleRows)
    df.protein <- as.data.table(as.matrix(df.protein))
  }
  
  ## Remove non-unique proteins and modified peptides if requested
  if (removeNonUniqueProteins){
    df.ptm <- df.ptm[!grepl(";", df.ptm$ProteinName),]
  }
  
  if (removeiRT){
    df.ptm <- df.ptm[!grepl("iRT", df.ptm$ProteinName),]
  }
  
  ## Format peptide data for locate peptide function
  df.ptm$PeptideSequence <- gsub("_", "", df.ptm$PeptideSequence)
  df.ptm$PeptideSequence <- gsub(paste0("\\[", modificationLabel, "\\]"),
                                 "*", df.ptm$PeptideSequence)
  
  ## Remove modifications not in modificationLabel
  df.ptm <- df.ptm[!grepl("\\[", df.ptm$PeptideSequence),]
  df.ptm$join_PeptideSequence <- gsub("\\*", "", df.ptm$PeptideSequence)
  
  locate_mod_df <- unique(df.ptm[, c("ProteinName", "PeptideSequence", 
                                     "join_PeptideSequence")])
  
  locate_mod_df <- locate_mod_df[grepl("\\*", locate_mod_df$PeptideSequence),]
  
  ## Load and format FASTA file
  if (identical(typeof(fasta), "character")){
    fasta <- tidyFasta(fasta)
  }
  formated_fasta <- as.data.table(fasta)
  
  min_len_peptide <- 6
  df.fasta.ptm <- merge(locate_mod_df, 
                        formated_fasta[, c("uniprot_iso", "sequence")], 
                        by.x = "ProteinName", by.y = "uniprot_iso")
  
  df.fasta.ptm <- df.fasta.ptm[
    which(nchar(df.fasta.ptm$join_PeptideSequence) > min_len_peptide & str_count(
      df.fasta.ptm$sequence, df.fasta.ptm$join_PeptideSequence) == 1),]
  
  start <- sapply(seq_len(nrow(df.fasta.ptm)),
                  function(i) gregexpr(df.fasta.ptm$PeptideSequence[i], 
                                       df.fasta.ptm$sequence[i])[[1]])
  mod_loc <- sapply(df.fasta.ptm$PeptideSequence, 
                    function(x) {gregexpr("\\*", x)})
  
  mod_index <- sapply(seq_along(mod_loc), function(i){mod_loc[i] <- list(
      as.integer(unlist(mod_loc[i])) + start[i][[1]][1])})
  
  peptide_mod <- mapply(spectro_get_sites, mod_loc, mod_index,
                        df.fasta.ptm$join_PeptideSequence)
  
  df.fasta.ptm$Site <- peptide_mod
  
  df.fasta.join <- unique(df.fasta.ptm[, c("ProteinName", "PeptideSequence", 
                                           "Site")])
  
  if (!removeNonUniqueProteins){
    add_non_unique <- unique(locate_mod_df[grepl(";", locate_mod_df$ProteinName), 
                                    c("ProteinName", "PeptideSequence")])
    add_non_unique$Site <- add_non_unique$PeptideSequence
    df.fasta.join <- rbindlist(list(df.fasta.join, add_non_unique))
  }
  
  #Data formatting for MSstatsLiP analysis
  MSstats_PTM <- merge(df.ptm, df.fasta.join,
                       by = c("ProteinName", "PeptideSequence"))
  MSstats_PTM$ProteinName <- paste(MSstats_PTM$ProteinName,
                                   MSstats_PTM$Site, sep = '_')
  
  MSstats_PTM$Intensity <- ifelse(MSstats_PTM$Intensity <= 1, NA, 
                                  MSstats_PTM$Intensity)
  MSstats_PTM[, join_PeptideSequence := NULL]
  
  if (!is.null(Protein.data)) {
    if (removeNonUniqueProteins){
      df.protein <- df.protein[!grepl(";", df.protein$ProteinName),]
    }
    # df.protein$PeptideSequence <- str_extract(df.protein$PeptideSequence,
    #                                       "([ACDEFGHIKLMNPQRSTVWY]+)")
    df.protein <- df.protein[nchar(PeptideSequence) > min_len_peptide]
    # 
    # df.protein$Intensity <- ifelse(df.protein$Intensity <= 1, NA, 
    #                                df.protein$Intensity)
    
    MSstats_Protein <- df.protein
  } else {
    MSstats_Protein <- NULL
  }
  
  PTMExpt <- list(
    PTM = MSstats_PTM,
    PROTEIN = MSstats_Protein
  )
  
  return(PTMExpt)
  
}