#' Convert output of TMT labeled MaxQuant experiment into MSstatsPTM format
#' 
#' Takes as input TMT experiments from MaxQ and converts the data into the 
#' format needed for MSstatsPTM. Requires only the modified file from MaxQ (for 
#' example Phospho(STY)Sites) and an annotation file for PTM data. To adjust 
#' modified peptides for changes in global protein level, unmodified TMT 
#' experimental data must also be returned.
#'
#' @export
#' @importFrom stringr str_extract regex str_replace fixed str_split
#' @importFrom data.table melt as.data.table `:=` `%like%` setcolorder setDT tstrsplit setcolorder
#' @importFrom MSstatsTMT MaxQtoMSstatsTMTFormat
#' @importFrom checkmate assertChoice assertLogical
#' 
#' @param sites.data modified peptide output from MaxQuant. For example, a
#' phosphorylation experiment would require the Phospho(STY)Sites.txt file
#' @param annotation.ptm data frame annotation file for the ptm level data.
#' Contains column Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition.
#' @param evidence for global protein dataset. name of 'evidence.txt' data, 
#' which includes feature-level data.
#' @param proteinGroups for global protein dataset, name of 'proteinGroups.txt' 
#' data.
#' @param annotation.prot data frame annotation file for the protein level data.
#' Contains column Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition.
#' @param mod.num For modified peptide dataset. The number modifications per 
#' peptide to be used. If "Single", only peptides with one modification will be 
#' used. Otherwise "Total" can be selected which does not cap the number of 
#' modifications per peptide. "Single" is the default. Selecting "Total" may 
#' confound the effect of different modifications.
#' @param TMT.keyword the sub-name of columns in sites.data file. Default is 
#' `TMT`. This corresponds to the columns in the format 
#' `Reporter.intensity.corrected.1.TMT1phos___1`. Specifically, this parameter 
#' indicates the first section of the string `TMT1phos` (Before the mixture 
#' number). If `TMT` is present in the string, set this value to `TMT`. Else if 
#' `TMT` is not there (ie string is in the format `1phos`) leave this parameter 
#' as an empty string ('').
#' @param ptm.keyword the sub-name of columns in the sites.data file. Default is 
#' `phos`. This corresponds to the columns in the format 
#' `Reporter.intensity.corrected.1.TMT1phos___1`. Specifically, this parameter 
#' indicates the second section of the string `TMT1phos` (After the mixture 
#' number). If the string is present, set this parameter. Else if this part of 
#' the string is empty (ie string is in the format `TMT1`) leave this parameter 
#' as an empty string ('').
#' @param which.proteinid.ptm For PTM dataset, which column to use for protein 
#' name. Use 'Proteins'(default) column for protein name. 'Leading.proteins' or 
#' 'Leading.razor.protein' or 'Gene.names' can be used instead to get the 
#' protein ID with single protein. However, those can potentially have the 
#' shared peptides.
#' @param which.proteinid.protein For Protein dataset, which column to use for 
#' protein name. Same options as above.
#' @param removeMpeptides If Oxidation (M) modifications should be removed. 
#' Default is TRUE.
#' @return a list of two data.tables named 'PTM' and 'PROTEIN' in the format 
#' required by MSstatsPTM.
#' @examples
#' 
#' head(raw.input.tmt$PTM)
#' head(raw.input.tmt$PROTEIN)
#' 
MaxQtoMSstatsPTMFormat <- function(sites.data,
                                   annotation.ptm,
                                   evidence = NULL,
                                   proteinGroups = NULL,
                                   annotation.prot = NULL,
                                   mod.num = 'Single',
                                   TMT.keyword = "TMT",
                                   ptm.keyword = "phos",
                                   which.proteinid.ptm = "Protein",
                                   which.proteinid.protein = 
                                     "Leading.razor.protein",
                                   removeMpeptides = FALSE) {
  
  .checkMaxQconverterParams(mod.num,
                            ptm.keyword,
                            which.proteinid.ptm,
                            which.proteinid.protein,
                            removeMpeptides)
  
  pho.data <- as.data.table(sites.data)
  annot.ptm <- as.data.table(annotation.ptm)
  annot <- as.data.table(annotation.prot)
  
  clean.prot <- .check.global.protein(evidence, proteinGroups)

  ## Format annotation for PTM
  setcolorder(annot.ptm, c("Run", "Channel", "Condition", "Mixture", 
                       "TechRepMixture", "Fraction", "BioReplicate"))
  #annot.ptm <- copy(annot)
  annot.ptm <- annot.ptm[, "Fraction" := NULL]
  annot.ptm <- annot.ptm[, "Run" :=NULL]
  annot.ptm <- unique(annot.ptm)
  annot.ptm$Replicate <- paste0("Reporter.intensity.corrected.",
                                gsub('channel.', '', annot.ptm$Channel), 
                                ".", TMT.keyword, 
                                gsub("mixture", "", annot.ptm$Mixture), 
                                ptm.keyword)
                                
    
  # if (keyword == "Plex"){
  #   annot.ptm$Replicate <- paste0("Reporter.intensity.corrected.",
  #                             gsub('channel.', '', annot.ptm$Channel),
  #                             ".",
  #                             paste0('Plex', annot.ptm$Mixture))
  # } else {
  #   annot.ptm$Replicate <- paste0("Reporter.intensity.corrected.",
  #                                 gsub('channel.', '', annot.ptm$Channel),
  #                                 ".",
  #                                 gsub('mixture', 'TMT', annot.ptm$Mixture),
  #                                 keyword)
  # }
  
  

  MSstatsPTMTMT.abun <- .convert.ptm.data(pho.data,
                                          annot.ptm,
                                          mod.num,
                                          ptm.keyword,
                                          which.proteinid.ptm,
                                          removeMpeptides)
  
  ## MaxQ phospho file has duplicate peptides per site 
  ## (if site has equal probability)
  ## Add protein into site to create unique identifier
  MSstatsPTMTMT.abun$PeptideSequence <- paste(
    MSstatsPTMTMT.abun$PeptideSequence, MSstatsPTMTMT.abun$ProteinName, 
    sep = ':')
  
  setDT(MSstatsPTMTMT.abun)[, PeptideSequence := tstrsplit(PeptideSequence, ":",
                                                           keep = 1)]
  
  MSstatsPTMformat <- list('PTM' = MSstatsPTMTMT.abun)

  if (clean.prot == TRUE){

    ## Clean raw data
    #evidence <- as.data.table(evidence)
    evidence <- evidence[!grepl("phos", evidence$Raw.file),]
    #proteinGroups <- as.data.table(proteinGroups)
    
    MSstatsTMT.abun <- MaxQtoMSstatsTMTFormat(evidence = evidence,
                                       proteinGroups = proteinGroups,
                                       annotation = annot,
                                      which.proteinid = which.proteinid.protein)

    MSstatsPTMformat <- list('PTM' = MSstatsPTMTMT.abun, 
                            "PROTEIN" = MSstatsTMT.abun)

  }

  return(MSstatsPTMformat)
}
