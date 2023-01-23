#' Convert Proteome Discoverer output into MSstatsPTM format
#' 
#' Import Proteome Discoverer files, identify modification site location.
#' 
#' @param input PD report corresponding with enriched experimental data.
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which 
#' includes Condition, BioReplicate, Run information. 'Run' will be matched 
#' with 'Spectrum.File'.
#' @param protein_input PD report corresponding with unmodified experimental 
#' data.
#' @param protein_annotation Same format as `annotation` corresponding to 
#' unmodified data.
#' @param mod_id Character that indicates the modification of interest. Default 
#' is `\\(Phospho\\)`. Note `\\` must be included before special characters.
#' @param use_unmod_peptides If `protein_input` is not provided, 
#' unmodified peptides can be extracted from `input` to be used in place of a 
#' global profiling run. Default is `FALSE`.
#' @param useNumProteinsColumn TRUE removes peptides which have more than 1 in 
#' Proteins column of PD output.
#' @param useUniquePeptide TRUE (default) removes peptides that are assigned 
#' for more than one proteins. We assume to use unique peptide for each protein.
#' @param summaryforMultipleRows max(default) or sum - when there are multiple 
#' measurements for certain feature and certain run, use highest or sum of 
#' multiple intensities.
#' @param removeFewMeasurements TRUE (default) will remove the features that 
#' have 1 or 2 measurements across runs.
#' @param removeOxidationMpeptides TRUE will remove the peptides including 
#' 'oxidation (M)' in modification. FALSE is default.
#' @param removeProtein_with1Peptide TRUE will remove the proteins which have 
#' only 1 peptide and charge. FALSE is default.
#' @param which_quantification Use 'Precursor.Area'(default) column for 
#' quantified intensities. 'Intensity' or 'Area' can be used instead.
#' @param which_proteinid Use 'Protein.Accessions'(default) column for protein 
#' name. 'Master.Protein.Accessions' can be used instead.
#' @param which_sequence Use 'Sequence'(default) column for peptide sequence. 
#' 'Annotated.Sequence' can be used instead.
#' @param use_log_file logical. If TRUE, information about data processing will 
#' be saved to a file.
#' @param append logical. If TRUE, information about data processing will be 
#' added to an existing log file.
#' @param verbose logical. If TRUE, information about data processing will be 
#' printed to the console
#' @param log_file_path character. Path to a file to which information about 
#' data processing will be saved. If not provided, such a file will be created 
#' automatically. If 'append = TRUE', has to be a valid path to a file.
#' 
#' @importFrom MSstats PDtoMSstatsFormat 
#' @importFrom stringr str_split str_trim
#' @importFrom data.table setnames
#' @return `list` of `data.table`
#' @export 
#' 
#' @examples
#' # The output should be in the following format.
#' head(raw.input$PTM)
#' head(raw.input$PROTEIN)
PDtoMSstatsPTMFormat = function(input,
                                annotation,
                                protein_input=NULL,
                                protein_annotation=NULL,
                                mod_id="\\(Phospho\\)",
                                use_unmod_peptides=FALSE,
                                useNumProteinsColumn = FALSE,
                                useUniquePeptide = TRUE,
                                summaryforMultipleRows = max,
                                removeFewMeasurements = TRUE,
                                removeOxidationMpeptides = FALSE,
                                removeProtein_with1Peptide = FALSE,
                                which_quantification = "Precursor.Area",
                                which_proteinid = "Protein.Group.Accessions",
                                which_sequence = "Sequence",
                                use_log_file = TRUE,
                                append = FALSE,
                                verbose = TRUE,
                                log_file_path = NULL){
  
  input = as.data.table(input)
  
  ##TODO: add more checks
  if (!is.null(protein_input) & use_unmod_peptides == TRUE){
    stop("Either pass protein_input data or set use_unmod_peptides = TRUE, not both")
  }
  
  ## Pull modifications from input and add into protein name
  input$mods = .extract_pd_mods(input$Modifications, mod_id)
  input[,which_proteinid] = paste(input[, which_proteinid, with=FALSE][[1]], 
                                  input[,"mods"][[1]], 
                                  sep="_")
  
  ptm_input = PDtoMSstatsFormat(input, annotation, useNumProteinsColumn,
                                useUniquePeptide, summaryforMultipleRows,
                                removeFewMeasurements, removeOxidationMpeptides,
                                removeProtein_with1Peptide, 
                                which_quantification, which_proteinid, 
                                which_sequence, use_log_file, append, verbose,
                                log_file_path)
  
  if ("PeptideModifiedSequence" %in% colnames(ptm_input)){
    setnames(ptm_input, c("PeptideModifiedSequence"), c("PeptideSequence"))
  }
  
  msstats_input = list(PTM = ptm_input)
  
  if (!is.null(protein_input)) {
    protein_input = PDtoMSstatsFormat(protein_input, protein_annotation, 
                                      useNumProteinsColumn,
                                      useUniquePeptide, summaryforMultipleRows,
                                      removeFewMeasurements, 
                                      removeOxidationMpeptides,
                                      removeProtein_with1Peptide, 
                                      which_quantification, which_proteinid, 
                                      which_sequence, use_log_file, append, 
                                      verbose, log_file_path)
    
    if ("PeptideModifiedSequence" %in% colnames(protein_input)){
      setnames(protein_input, c("PeptideModifiedSequence"), 
               c("PeptideSequence"))
    }
    
    msstats_input = list(PTM = ptm_input, PROTEIN = protein_input)
  }
  
  if (use_unmod_peptides){
    protein_input = ptm_input[!grepl(mod_id, ptm_input$PeptideSequence),]
    ptm_input = ptm_input[grepl(mod_id, ptm_input$PeptideSequence),]
    
    msstats_input = list(PTM = ptm_input, PROTEIN = protein_input)
  }
  
  return(msstats_input)
  
}

#' Convert Peaks Studio output into MSstatsPTM format
#' 
#' Currently only supports label-free quantification.
#' 
#' @param input name of Peaks Studio PTM output
#' @param annotation name of annotation file which includes Raw.file, Condition,
#' BioReplicate, Run. For example annotation see example below.
#' @param input_protein name of Peaks Studio unmodified protein output 
#' (optional)
#' @param annotation_protein name of annotation file which includes Raw.file, 
#' Condition,
#' BioReplicate, Run for unmodified protein output.
#' @param use_unmod_peptides Boolean if the unmodified peptides in the input 
#' file should be used to construct the unmodified protein output. Only used if
#' `input_protein` is not provided. Default is `FALSE`
#' @param target_modification Character name of modification of interest. To 
#' use all mod types, leave as `NULL`. Default is `NULL`. Note that if the name 
#' includes special characters, you must include "\\" before the characters. Ex.
#' "Phosphorylation \\(STY\\)"
#' @param remove_oxidation_peptides Boolean if Oxidation (M) modifications 
#' should be removed. Default is `FALSE`
#' @param remove_multi_mod_types Used if `target_modification` is not `NULL`. 
#' `TRUE` will remove peptides with multiple types of modifications 
#' (ie acetylation and phosphorylation). `FALSE` will keep these peptides and 
#' summarize them seperately.
#' @param summaryforMultipleRows max(default) or sum - when there are multiple 
#' measurements for certain feature and certain run, use highest or sum of 
#' multiple intensities.
#' @param use_log_file logical. If TRUE, information about data processing will 
#' be saved to a file.
#' @param append logical. If TRUE, information about data processing will be 
#' added to an existing log file.
#' @param verbose logical. If TRUE, information about data processing wil be 
#' printed to the console.
#' @param log_file_path character. Path to a file to which information about 
#' data processing will be saved. If not provided, such a file will be created 
#' automatically. If 'append = TRUE', has to be a valid path to a file.
#' 
#' @importFrom data.table as.data.table melt
#' @importFrom MSstatsConvert MSstatsBalancedDesign
#' 
#' @return `list` of `data.table`
#' @export 
#' 
#' @examples
#' # The output should be in the following format.
#' head(raw.input$PTM)
#' head(raw.input$PROTEIN)
PStoMSstatsPTMFormat = function(
    input, 
    annotation, 
    input_protein = NULL, 
    annotation_protein = NULL,
    use_unmod_peptides = FALSE, 
    target_modification = NULL, 
    remove_oxidation_peptides = FALSE, 
    remove_multi_mod_types = FALSE,
    summaryforMultipleRows = max, 
    use_log_file = TRUE, 
    append = FALSE, 
    verbose = TRUE, 
    log_file_path = NULL
){
  
  input = as.data.table(input)
  if (!is.null(input_protein)){
    input_protein = as.data.table(input_protein)
  }
  
  if (!"Source.File" %in% colnames(input)){
    message("Pivoting input data..")
    input = .pivotPS(input)
  }
  
  message("Merging with annotation..")
  input = merge(input, annotation, all.x = TRUE, by = "Source.File")
  
  ## Filter for required modifications
  message("Filtering modifications based on function arguements..")
  if (remove_oxidation_peptides){
    input = input[!grepl("Oxidation", Mod)]
  }
  
  if (is.null(input_protein) & use_unmod_peptides){
    input_prot = input[is.na(Mod)]
  }
  
  if (!is.null(target_modification)){
    input = input[grepl(target_modification, Mod)]
    
    if (nrow(input) == 0){
      stop("No modifications match target. Please ensure target modification is correctly spelled.")
    }
  } else {
    input = input[!is.na(Mod)]
  }
  
  ## Can be peptides where two mods such as Phosphorylation (STY); Deamidation (NQ)
  ## Setting this will remove these peptides
  if (remove_multi_mod_types){
    input = input[!grepl(";", Mod)]
  }
  
  ## Add mod to protein name
  input$PeptideModifiedSequence = input$PeptideSequence
  input = MSstatsPTMSiteLocator(input, fasta_file=NULL, mod_id="\\*", 
                                mod_id_is_numeric=TRUE, 
                                terminus_included=TRUE, terminus_id="\\.")
  
  ## Finish converter
  input[, c("Raw.File", "Mod", "End", "Start"):=NULL]
  feature_cols = c("ProteinName", "PeptideSequence", "FragmentIon", 
                   "ProductCharge", "PrecursorCharge", "IsotopeLabelType",
                   "BioReplicate", "Condition", "Run" )
  input = MSstatsBalancedDesign(input, c('PeptideSequence', 'ProductCharge'))
  input = as.data.table(input)[, list(
    Intensity = summaryforMultipleRows(Intensity, na.rm = TRUE)), 
    by = feature_cols]

  if (!is.null(input_protein)){
    input_protein[, c("Raw.File", "Mod", "End", "Start"):=NULL]
    input_protein = MSstatsBalancedDesign(input_protein, 
                                       c('PeptideSequence', 'ProductCharge'))
    input_protein = as.data.table(input_protein)[, list(
      Intensity = summaryforMultipleRows(Intensity, na.rm = TRUE)), 
      by = feature_cols]
    
    input_protein = as.data.frame(input_protein)
  }
  
  return(list("PTM" = as.data.frame(input),
              "PROTEIN" = input_protein))
}

#' Convert Skyline output into MSstatsPTM format
#' 
#' Currently only supports label-free quantification.
#' 
#' @param input name of Skyline PTM output
#' @param fasta_path A string of path to a FASTA file, used to match PTM peptides.
#' @param annotation name of 'annotation.txt' data which includes Condition, 
#' BioReplicate, Run. If annotation is already complete in Skyline, use 
#' annotation=NULL (default). It will use the annotation information from input.
#' @param input_protein name of Skyline unmodified protein output (optional)
#' @param annotation_protein name of 'annotation.txt' data which includes Condition, 
#' BioReplicate, Run for unmodified protein output. This can be the same as 
#' `annotation`.
#' @param use_unmod_peptides Boolean if the unmodified peptides in the input 
#' file should be used to construct the unmodified protein output. Only used if
#' `input_protein` is not provided. Default is `FALSE`.
#' @param removeiRT TRUE (default) will remove the proteins or peptides which 
#' are labeld 'iRT' in 'StandardType' column. FALSE will keep them.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that 
#' have greater than qvalue_cutoff in DetectionQValue column. Those intensities 
#' will be replaced with zero and will be considered as censored missing values 
#' for imputation purpose.
#' @param qvalue_cutoff Cutoff for DetectionQValue. default is 0.01.
#' @param use_unique_peptide TRUE (default) removes peptides that are assigned 
#' for more than one proteins. We assume to use unique peptide for each protein.
#' @param remove_few_measurements TRUE will remove the features that 
#' have 1 or 2 measurements across runs. FALSE is default.
#' @param remove_oxidation_peptides TRUE will remove the peptides including 
#' 'oxidation (M)' in modification. FALSE is default.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have 
#' only 1 feature, which is the combination of peptide, precursor charge, 
#' fragment and charge. FALSE is default.
#' @param use_log_file logical. If TRUE, information about data processing will 
#' be saved to a file.
#' @param append logical. If TRUE, information about data processing will be 
#' added to an existing log file.
#' @param verbose logical. If TRUE, information about data processing wil be 
#' printed to the console.
#' @param log_file_path character. Path to a file to which information about 
#' data processing will be saved. If not provided, such a file will be created 
#' automatically. If 'append = TRUE', has to be a valid path to a file.
#' 
#' @importFrom data.table as.data.table setnames
#' @importFrom MSstats SkylinetoMSstatsFormat
#' 
#' @return `list` of `data.table`
#' @export
#' 
#' @examples
#' # The output should be in the following format.
#' head(raw.input$PTM)
#' head(raw.input$PROTEIN)
SkylinetoMSstatsPTMFormat = function(input,
                                     fasta_path,
                                     annotation=NULL,
                                     input_protein=NULL,
                                     annotation_protein=NULL,
                                     use_unmod_peptides = FALSE, 
                                     removeiRT = TRUE,
                                     filter_with_Qvalue = TRUE,
                                     qvalue_cutoff = 0.01,
                                     use_unique_peptide = TRUE,
                                     remove_few_measurements = FALSE,
                                     remove_oxidation_peptides = FALSE,
                                     removeProtein_with1Feature = FALSE, 
                                     use_log_file = TRUE, append = FALSE, 
                                     verbose = TRUE, log_file_path = NULL
){
  
  message("Converting modified data..")
  input = as.data.table(input)
  input_prot = NULL
  if (!is.null(annotation)){
    input = merge(input, annotation, all.x = TRUE, by = "File Name")
  }
  
  ## Filter for required modifications
  message("Filtering modifications based on function arguements..")
  if (is.null(input_protein) & use_unmod_peptides){
    input_prot = input[!grepl("\\[+", input$`Peptide Modified Sequence`), ]
    
    if (nrow(input_prot)){
      stop("No unmodified peptides in PTM dataset. It cannot be used for \
           adjustment.")
    }
  }
  
  ## Ensure only modified peptides are used
  input = input[grepl("\\[+", input$`Peptide Modified Sequence`), ]

  setnames(input, c("Protein Name", "Peptide Sequence", 
                    "Peptide Modified Sequence"), 
           c("ProteinName","PeptideSequence", "PeptideModifiedSequence"))
  
  ## Locate sites
  input = MSstatsPTMSiteLocator(input, fasta_file=fasta_path, mod_id_is_numeric=TRUE, 
                                terminus_included=FALSE)

  
  input = SkylinetoMSstatsFormat(input, annotation = NULL, removeiRT, 
                         filter_with_Qvalue, qvalue_cutoff,
                         use_unique_peptide, remove_few_measurements,
                         remove_oxidation_peptides, removeProtein_with1Feature,
                         use_log_file, append, verbose, log_file_path)
  
  if (!is.null(input_protein)){
    message("Converting unmodified protein data..")
    
    input_prot = SkylinetoMSstatsFormat(input_protein, annotation_protein, 
                                        removeiRT, filter_with_Qvalue, 
                                        qvalue_cutoff, use_unique_peptide, 
                                        remove_few_measurements, 
                                        remove_oxidation_peptides, 
                                        removeProtein_with1Feature, 
                                        use_log_file, append, verbose, 
                                        log_file_path)
  }
  
  return(list("PTM" = as.data.frame(input),
              "PROTEIN" = as.data.frame(input_prot)))
}

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
#' @importFrom MSstats MaxQtoMSstatsFormat
#' @importFrom checkmate assertChoice assertLogical
#' 
#' @param evidence name of 'evidence.txt' data, which includes feature-level 
#' data for enriched (PTM) data.
#' @param annotation_ptm data frame annotation file for the ptm level data.
#' Contains column Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition.
#' @param fasta_path A string of path to a FASTA file, used to match PTM peptides.
#' @param fasta_protein_name Name of fasta column that matches with protein name
#' in evidence file. Default is `uniprot_ac`.
#' @param mod_id Character that indicates the modification of interest. Default 
#' is `\\(Phospho\\)`. Note `\\` must be included before special characters.
#' @param sites_data (Not recommended. Only used if evidence file not provided. 
#' Only works for TMT labeled data) Modified peptide output from MaxQuant. For 
#' example, a phosphorylation experiment would require the Phospho(STY)Sites.txt
#' file
#' @param evidence_prot name of 'evidence.txt' data, which includes 
#' feature-level data for global profiling (unmodified) data.
#' @param proteinGroups name of 'proteinGroups.txt' data. It needs to matching 
#' protein group ID in `evidence_prot`.
#' @param annotation_prot data frame annotation file for the protein level data.
#' Contains column Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition.
#' @param use_unmod_peptides Boolean if the unmodified peptides in the input 
#' file should be used to construct the unmodified protein output. Only used if
#' `input_protein` is not provided. Default is `FALSE`.
#' @param labeling_type Either `TMT` or `LF` (Label-Free) depending on 
#' experimental design. Default is `LF`.
#' @param mod_num (Only if `sites.data` is used) For modified peptide dataset. 
#' The number modifications per peptide to be used. If "Single", only peptides 
#' with one modification will be used. Otherwise "Total" can be selected which 
#' does not cap the number of modifications per peptide. "Single" is the 
#' default. Selecting "Total" may confound the effect of different 
#' modifications.
#' @param TMT_keyword (Only if `sites.data` is used) the sub-name of columns 
#' in sites.data file. Default is `TMT`. This corresponds to the columns in the 
#' format `Reporter.intensity.corrected.1.TMT1phos___1`. Specifically, this 
#' parameter indicates the first section of the string `TMT1phos` (Before the 
#' mixture number). If `TMT` is present in the string, set this value to `TMT`. 
#' Else if `TMT` is not there (ie string is in the format `1phos`) leave this 
#' parameter as an empty string ('').
#' @param ptm_keyword (Only if `sites.data` is used) the sub-name of columns in 
#' the sites.data file. Default is 
#' `phos`. This corresponds to the columns in the format 
#' `Reporter.intensity.corrected.1.TMT1phos___1`. Specifically, this parameter 
#' indicates the second section of the string `TMT1phos` (After the mixture 
#' number). If the string is present, set this parameter. Else if this part of 
#' the string is empty (ie string is in the format `TMT1`) leave this parameter 
#' as an empty string ('').
#' @param which_proteinid_ptm For PTM dataset, which column to use for protein 
#' name. Use 'Proteins'(default) column for protein name. 'Leading.proteins' or 
#' 'Leading.razor.protein' or 'Gene.names' can be used instead to get the 
#' protein ID with single protein. However, those can potentially have the 
#' shared peptides.
#' @param which_proteinid_protein For Protein dataset, which column to use for 
#' protein name. Same options as above.
#' @param removeMpeptides If Oxidation (M) modifications should be removed. 
#' Default is TRUE.
#' @param removeOxidationMpeptides TRUE will remove the peptides including 
#' 'oxidation (M)' in modification. FALSE is default.
#' @param removeProtein_with1Peptide TRUE will remove the proteins which have 
#' only 1 peptide and charge. FALSE is default.
#' @param use_log_file logical. If TRUE, information about data processing will 
#' be saved to a file.
#' @param append logical. If TRUE, information about data processing will be 
#' added to an existing log file.
#' @param verbose logical. If TRUE, information about data processing wil be 
#' printed to the console.
#' @param log_file_path character. Path to a file to which information about 
#' data processing will be saved. If not provided, such a file will be created 
#' automatically. If 'append = TRUE', has to be a valid path to a file.
#' 
#' @return a list of two data.tables named 'PTM' and 'PROTEIN' in the format 
#' required by MSstatsPTM.
#' 
#' @examples
#' head(raw.input.tmt$PTM)
#' head(raw.input.tmt$PROTEIN)
MaxQtoMSstatsPTMFormat = function(evidence=NULL,
                                  annotation_ptm=NULL,
                                  fasta_path=NULL,
                                  fasta_protein_name="uniprot_ac",
                                  mod_id="\\(Phospho \\(STY\\)\\)",
                                  sites_data=NULL,
                                  evidence_prot = NULL,
                                  proteinGroups = NULL,
                                  annotation_prot = NULL,
                                  use_unmod_peptides=FALSE,
                                  labeling_type = "LF",
                                  mod_num = 'Single',
                                  TMT_keyword = "TMT",
                                  ptm_keyword = "phos",
                                  which_proteinid_ptm = "Proteins",
                                  which_proteinid_protein = "Proteins",
                                  removeMpeptides = FALSE,
                                  removeOxidationMpeptides = FALSE,
                                  removeProtein_with1Peptide = FALSE,
                                  use_log_file = TRUE,
                                  append = FALSE,
                                  verbose = TRUE,
                                  log_file_path = NULL) {
  
  ## TODO: add this code to check function
  if (is.null(evidence) & is.null(sites_data)){
    stop("Either evidence and proteinGroups or sites_data files must be included.")
  }
  
  # .checkMaxQconverterParams(mod.num,
  #                           ptm.keyword,
  #                           which.proteinid.ptm,
  #                           which.proteinid.protein,
  #                           removeMpeptides)
  
  if (is.null(sites_data)){
    evidence = as.data.table(evidence)
  } else {
    pho.data = as.data.table(sites_data)
  }
  
  annot.ptm = as.data.table(annotation_ptm)
  
  # clean.prot = .check.global.protein(evidence, proteinGroups)
  
  ## Format annotation for PTM
  if (is.null(sites_data)){
    if(labeling_type == "TMT"){
      
      evidence_sites = MSstatsPTMSiteLocator(evidence, 
                                        protein_name_col = which_proteinid_ptm,
                                        unmod_pep_col = "Sequence",
                                        mod_pep_col = "Modified.sequence",
                                        fasta_file = fasta_path,
                                        fasta_protein_name = "uniprot_ac", 
                                        mod_id=mod_id, 
                                        mod_id_is_numeric = FALSE, 
                                        terminus_included = FALSE)
      
      msstatsptm_input = MaxQtoMSstatsTMTFormatHelper(evidence_sites,
                       annotation_ptm,
                       which.proteinid = which_proteinid_ptm,
                       rmPSM_withfewMea_withinRun = removeProtein_with1Peptide,
                       use_log_file=use_log_file,
                       append=append,
                       verbose=verbose,
                       log_file_path=log_file_path)
      
    } else if (labeling_type == "LF"){
      
      evidence_sites = MSstatsPTMSiteLocator(evidence, 
                                         protein_name_col= which_proteinid_ptm,
                                         unmod_pep_col = "Sequence",
                                         mod_pep_col = "Modified.sequence",
                                         fasta_file=fasta_path,
                                         fasta_protein_name="uniprot_ac", 
                                         mod_id=mod_id, 
                                         mod_id_is_numeric=FALSE, 
                                         terminus_included=FALSE)
      
      msstatsptm_input = MaxQtoMSstatsFormatHelper(evidence_sites,
                       annotation_ptm,
                       proteinID=which_proteinid_ptm,
                       removeMpeptides = removeMpeptides,
                       removeOxidationMpeptides = removeOxidationMpeptides,
                       removeProtein_with1Peptide = removeProtein_with1Peptide,
                       use_log_file=use_log_file,
                       append=append,
                       verbose=verbose,
                       log_file_path=log_file_path)

    
    }
  } else {
    setcolorder(annot.ptm, c("Run", "Channel", "Condition", "Mixture", 
                             "TechRepMixture", "Fraction", "BioReplicate"))
  
    annot.ptm = annot.ptm[, "Fraction" := NULL] ## TODO: should this be removed?
    annot.ptm = annot.ptm[, "Run" :=NULL] ## We recreate this column to match with PTM data
    annot.ptm = unique(annot.ptm)
    annot.ptm$Replicate = paste0("Reporter.intensity.corrected.",
                                 gsub('channel.', '', annot.ptm$Channel), 
                                 ".", TMT.keyword, 
                                 gsub("mixture", "", annot.ptm$Mixture), 
                                 ptm.keyword)
    
    msstatsptm_input = .convert.ptm.data(pho.data,
                                           annot.ptm,
                                           mod.num,
                                           ptm.keyword,
                                           which.proteinid.ptm,
                                           removeMpeptides)
    
    ## MaxQ phospho file has duplicate peptides per site 
    ## (if site has equal probability)
    ## Add protein into site to create unique identifier
    msstatsptm_input$PeptideSequence = paste(
      msstatsptm_input$PeptideSequence, msstatsptm_input$ProteinName, 
      sep = ':')
    
    setDT(msstatsptm_input)[, PeptideSequence := tstrsplit(PeptideSequence, ":",
                                                             keep = 1)]
  }
  
  MSstatsPTMformat = list('PTM' = msstatsptm_input)
  
  if (!is.null(evidence_prot)){
    annot = as.data.table(annotation.prot)
    
    ## Clean raw data
    #evidence = as.data.table(evidence)
    # evidence = evidence[!grepl("phos", evidence$Raw.file),]
    #proteinGroups = as.data.table(proteinGroups)
    
    if(labeling_type == "TMT"){
      msstats.abun = MaxQtoMSstatsTMTFormat(evidence = evidence_prot,
                                            proteinGroups = proteinGroups_prot,
                                            annotation = annotation_prot,
                                      which.proteinid = which_proteinid_protein)
    } else if (labeling_type == "LF"){
      msstats.abun = MaxQtoMSstatsFormat(evidence = evidence_prot,
                                            proteinGroups = proteinGroups_prot,
                                            annotation = annotation_prot,
                                      which.proteinid = which_proteinid_protein)
    }
    
    MSstatsPTMformat = list('PTM' = MSstatsPTMTMT.abun, 
                            "PROTEIN" = msstats.abun)
    
  }
  
  if (use_unmod_peptides){
    msstats.abun = msstatsptm_input[!grepl(mod_id, msstatsptm_input$PeptideSequence),]
    msstatsptm_input = msstatsptm_input[grepl(mod_id, msstatsptm_input$PeptideSequence),]
    
    MSstatsPTMformat = list(PTM = msstatsptm_input, 
                            PROTEIN = msstats.abun)
  }
  
  return(MSstatsPTMformat)
}

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
#' includes Condition, BioReplicate, Run, and Type (PTM or Protein) 
#' information. It will be matched with the column name of input for MS runs. 
#' Please note PTM and global Protein run names are often different, which is 
#' why an additional Type column indicating Protein or PTM is required.
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
#' @param removeFewMeasurements TRUE (default) will remove the features that 
#' have 1 or 2 measurements across runs.
#' @param removeOxidationMpeptides TRUE will remove the modified peptides
#' including 'Oxidation (M)' sequence. FALSE is default.
#' @param removeProtein_with1Peptide TRUE will remove the proteins which have
#' only 1 peptide and charge. FALSE is default.
#' @param mod.num For modified peptide dataset, must be one of `Single` or 
#' `Total`. The default is `Single`. The number modifications per peptide to be 
#' used. If "Single", only peptides with one modification will be 
#' used. Otherwise "Total" includes peptides with more than one modification.
#' Selecting "Total" may confound the effect of different modifications.
#' @return a list of two data.tables named 'PTM' and 'PROTEIN' in the format 
#' required by MSstatsPTM.
#' @examples
#' 
#' # Example annotation file
#' annotation = data.frame('Condition' = c('Control', 'Control', 'Control',
#'                          'Treatment', 'Treatment', 'Treatment'),
#'                          'BioReplicate' = c(1,2,3,4,5,6),
#'                          'Run' = c('prot_run_1', 'prot_run_2', 'prot_run_3',
#'                                   'phos_run_1', 'phos_run_2', 'phos_run_3'),
#'                          'Type' = c("Protein", "Protein", "Protein", "PTM", 
#'                                     "PTM", "PTM"))
#'                                     
#' # The output should be in the following format.
#' head(raw.input$PTM)
#' head(raw.input$PROTEIN)
ProgenesistoMSstatsPTMFormat = function(ptm_input,
                                        annotation,
                                        global_protein_input = FALSE,
                                        fasta_path = FALSE,
                                        useUniquePeptide=TRUE,
                                        summaryforMultipleRows=max, 
                                        removeFewMeasurements=TRUE,
                                        removeOxidationMpeptides=FALSE,
                                        removeProtein_with1Peptide=FALSE,
                                        mod.num = 'Single'){
  
  ptm_input = as.data.table(ptm_input)
  annotation = as.data.table(annotation)
  col_order = colnames(ptm_input)
  ptm_input$id = 1:nrow(ptm_input)
  
  annotation = as.data.table(annotation)
  ptm_annot = annotation[Type == "PTM"]
  protein_annot = annotation[Type == "Protein"]
  
  X.10 = X.9 = X.8 = NULL
  
  ## Format PTM data
  if (fasta_path == FALSE) {
    ptm_input[,X.10 := paste(X.10, X.8, X.9, sep = "_")]
    ptm_input[2, "X.10"] = "Accession"
  } else {
    .progensis.add.sites(ptm_input, fasta_path, col_order)
  }
  
  convert.ptm = ProgenesistoMSstatsFormat(ptm_input, ptm_annot,
                                          useUniquePeptide,
                                          summaryforMultipleRows,
                                          removeFewMeasurements,
                                          removeOxidationMpeptides,
                                          removeProtein_with1Peptide)
  
  if (global_protein_input[[1]][1] != FALSE){
    
    global_protein_input = as.data.table(global_protein_input)
    
    convert.prot = ProgenesistoMSstatsFormat(global_protein_input, 
                                             protein_annot,
                                             useUniquePeptide,
                                             summaryforMultipleRows,
                                             removeFewMeasurements,
                                             removeOxidationMpeptides,
                                             removeProtein_with1Peptide)
    MSstatsPTM.data = list("PTM" = convert.ptm,
                           "PROTEIN" = convert.prot)
  } else {
    MSstatsPTM.data = list("PTM" = convert.ptm,
                           "PROTEIN" = NULL)
  }
  
  return(MSstatsPTM.data)
  
}


#' Convert Spectronaut output into MSstatsPTM format
#' 
#' Currently only supports label-free quantification.
#' 
#' @param input name of Spectronaut PTM output, which is long-format. 
#' ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, 
#' IsotopeLabelType, Condition, BioReplicate, Run, Intensity, 
#' F.ExcludedFromQuantification are required. Rows with 
#' F.ExcludedFromQuantification=True will be removed.
#' @param annotation name of 'annotation.txt' data which includes Condition, 
#' BioReplicate, Run. If annotation is already complete in Spectronaut, 
#' use annotation=NULL (default). It will use the annotation information from 
#' input.
#' @param fasta_path string containing path to the corresponding fasta file for 
#' the modified peptide dataset.
#' @param protein_input name of Spectronaut global protein output, which is 
#' as in the same format as `input` parameter. 
#' @param protein_annotation name of annotation file for global protein data, in
#' the same format as above.
#' @param use_unmod_peptides If `protein_input` is not provided, 
#' unmodified peptides can be extracted from `input` to be used in place of a 
#' global profiling run. Default is `FALSE`.
#' @param intensity 'PeakArea'(default) uses not normalized peak area. 
#' 'NormalizedPeakArea' uses peak area normalized by Spectronaut. Default is 
#' NULL
#' @param mod_id Character that indicates the modification of interest. Default 
#' is `\\(Phospho\\)`. Note `\\` must be included before special characters.
#' @param fasta_protein_name Name of fasta column that matches with protein name
#' in evidence file. Default is `uniprot_iso`.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that 
#' have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will 
#' be replaced with zero and will be considered as censored missing values for 
#' imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. Default is 0.01.
#' @param useUniquePeptide TRUE (default) removes peptides that are assigned for
#'  more than one proteins. We assume to use unique peptide for each protein.
#' @param removeFewMeasurements TRUE (default) will remove the features that 
#' have 1 or 2 measurements across runs.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have 
#' only 1 feature, which is the combination of peptide, precursor charge, 
#' fragment and charge. FALSE is default.
#' @param summaryforMultipleRows max(default) or sum - when there are multiple 
#' measurements for certain feature and certain run, use highest or sum of 
#' multiple intensities.
#' @param use_log_file logical. If TRUE, information about data processing will 
#' be saved to a file.
#' @param append logical. If TRUE, information about data processing will be 
#' added to an existing log file.
#' @param verbose logical. If TRUE, information about data processing wil be 
#' printed to the console.
#' @param log_file_path character. Path to a file to which information about 
#' data processing will be saved. If not provided, such a file will be created 
#' automatically. If 'append = TRUE', has to be a valid path to a file.
#' @return a list of two data.tables named 'PTM' and 'PROTEIN' in the format 
#' required by MSstatsPTM.
#' @export
#' @importFrom MSstats SpectronauttoMSstatsFormat
#' @examples
#' # The output should be in the following format.
#' head(raw.input$PTM)
#' head(raw.input$PROTEIN)
SpectronauttoMSstatsPTMFormat = function(
    input,
    annotation = NULL,
    fasta_path = NULL,
    protein_input = NULL,
    protein_annotation = NULL,
    use_unmod_peptides=FALSE,
    intensity = "PeakArea",
    mod_id="\\[Phospho \\(STY\\)\\]",
    fasta_protein_name="uniprot_iso",
    filter_with_Qvalue = TRUE,
    qvalue_cutoff = 0.01,
    useUniquePeptide = TRUE,
    removeFewMeasurements = TRUE,
    removeProtein_with1Feature = FALSE,
    summaryforMultipleRows = max,
    use_log_file = TRUE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = NULL){
  
  ## TODO: Add checks on input and params
  input = as.data.table(input)
  ## Needed if input only call PSM
  if (!"EG.ModifiedSequence" %in% colnames(input)){
    input$PeptideSequence = .MSstatsPTMRemoveMods(input$EG.PrecursorId)
    mod_col = "EG.PrecursorId"
  } else {
    input$PeptideSequence = .MSstatsPTMRemoveMods(input$EG.ModifiedSequence)
    mod_col = "EG.ModifiedSequence"
  }
  
  
  input = MSstatsPTMSiteLocator(input, protein_name_col= "PG.ProteinGroups",
                                unmod_pep_col = "PeptideSequence",
                                mod_pep_col = mod_col,
                                # clean_mod=TRUE,
                                fasta_file=fasta_path, 
                                fasta_protein_name=fasta_protein_name,
                                mod_id=mod_id, 
                                mod_id_is_numeric=FALSE, 
                                terminus_included=FALSE, terminus_id="\\.")
  
  ptm_input = SpectronauttoMSstatsFormat(input, annotation, intensity, 
                                         filter_with_Qvalue,
                                         qvalue_cutoff, useUniquePeptide, 
                                         removeFewMeasurements, 
                                         removeProtein_with1Feature, 
                                         summaryforMultipleRows)
  
  msstats_input = list(PTM = ptm_input)
  if (!is.null(protein_input)) {
    protein_input = SpectronauttoMSstatsFormat(protein_input, 
                                               protein_annotation, intensity, 
                                               filter_with_Qvalue,
                                               qvalue_cutoff, useUniquePeptide, 
                                               removeFewMeasurements, 
                                               removeProtein_with1Feature, 
                                               summaryforMultipleRows)
    
    msstats_input = list(PTM = ptm_input, PROTEIN = protein_input)
  }
  
  if (use_unmod_peptides){
    protein_input = ptm_input[!grepl(mod_id, ptm_input$PeptideSequence),]
    ptm_input = ptm_input[grepl(mod_id, ptm_input$PeptideSequence),]
    
    msstats_input = list(PTM = ptm_input, PROTEIN = protein_input)
  }
  
  return(msstats_input)

}

## TODO Old Spect converter.. reconcile with new one and delete
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
#' @param fasta_path A string of path to a FASTA file, used to match PTM peptides.
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
#' @param modificationLabel String of modification name in `EG.ModifiedSequence`
#' column. Default is 'Phospho'. Note label could also include amino acids in 
#' abbreviations, for example 'Phospho (STY)'. If special characters are in the 
#' name, you must add // to the label, ex. 'Phospho \\(STY\\)'.
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
# SpectronauttoMSstatsPTMFormat = function(PTM.data,
#                                          fasta,
#                                          Protein.data = NULL,
#                                          annotation = NULL,
#                                          intensity = 'PeakArea',
#                                          filter_with_Qvalue = TRUE,
#                                          qvalue_cutoff = 0.01,
#                                          useUniquePeptide = TRUE,
#                                          removeFewMeasurements = TRUE,
#                                          removeProtein_with1Feature = FALSE,
#                                          removeNonUniqueProteins = TRUE,
#                                          modificationLabel = "Phospho",
#                                          removeiRT = TRUE,
#                                          summaryforMultipleRows=max,
#                                          which.Conditions = 'all'){
#   
#   ## TODO: Add logging
#   ## Check variable input
#   
#   ## Ensure format of input data
#   PTM.data = as.data.table(PTM.data)
#   if (!is.null(Protein.data)){
#     Protein.data = as.data.table(Protein.data)
#   }
#   
#   ## Check and filter for available conditions
#   if (which.Conditions != 'all') {
#     
#     PTM_conditions = unique(PTM.data$R.Condition)
#     if (!is.null(Protein.data)){
#       Protein_conditions = unique(Protein.data$R.Condition)
#     } else {
#       Protein_conditions = PTM_conditions
#     }
#     if((length(setdiff(PTM_conditions, which.Conditions)
#     ) == length(PTM_conditions)) |
#     (length(setdiff(Protein_conditions, which.Conditions)
#     ) == length(Protein_conditions))){
#       msg = (paste("None of the conditions specified in which.Conditions are",
#                    "available in one or both of LiP/TrP datasets. Please",
#                    "ensure the conditions listed in which.Conditions appear",
#                    "in the input datasets"))
#       stop(msg)
#     }
#     
#     PTM_conditions = PTM_conditions[(R.Condition %in% which.Conditions)]
#     if (!is.null(Protein_conditions)){
#       Protein_conditions = Protein_conditions[
#         (R.Condition %in% which.Conditions)]
#     }
#   }
#   
#   ## MSstats process
#   df.ptm = SpectronauttoMSstatsFormat(PTM.data, annotation, intensity,
#                                       filter_with_Qvalue, qvalue_cutoff,
#                                       useUniquePeptide, removeFewMeasurements,
#                                       removeProtein_with1Feature,
#                                       summaryforMultipleRows)
#   df.ptm = as.data.table(as.matrix(df.ptm))
#   if (!is.null(Protein.data)){
#     df.protein = SpectronauttoMSstatsFormat(as.data.frame(Protein.data), 
#                                             annotation, intensity,
#                                             filter_with_Qvalue, qvalue_cutoff,
#                                             useUniquePeptide, removeFewMeasurements,
#                                             removeProtein_with1Feature,
#                                             summaryforMultipleRows)
#     df.protein = as.data.table(as.matrix(df.protein))
#   }
#   
#   ## Remove non-unique proteins and modified peptides if requested
#   if (removeNonUniqueProteins){
#     df.ptm = df.ptm[!grepl(";", df.ptm$ProteinName),]
#   }
#   
#   if (removeiRT){
#     df.ptm = df.ptm[!grepl("iRT", df.ptm$ProteinName),]
#   }
#   
#   ## Format peptide data for locate peptide function
#   df.ptm$PeptideSequence = gsub("_", "", df.ptm$PeptideSequence)
#   df.ptm$PeptideSequence = gsub(paste0("\\[", modificationLabel, "\\]"),
#                                 "*", df.ptm$PeptideSequence)
#   
#   ## Remove modifications not in modificationLabel
#   df.ptm = df.ptm[!grepl("\\[", df.ptm$PeptideSequence),]
#   df.ptm$join_PeptideSequence = gsub("\\*", "", df.ptm$PeptideSequence)
#   
#   locate_mod_df = unique(df.ptm[, c("ProteinName", "PeptideSequence", 
#                                     "join_PeptideSequence")])
#   
#   locate_mod_df = locate_mod_df[grepl("\\*", locate_mod_df$PeptideSequence),]
#   
#   ## Load and format FASTA file
#   if (identical(typeof(fasta), "character")){
#     fasta = tidyFasta(fasta)
#   }
#   formated_fasta = as.data.table(fasta)
#   
#   min_len_peptide = 6
#   df.fasta.ptm = merge(locate_mod_df, 
#                        formated_fasta[, c("uniprot_iso", "sequence")], 
#                        by.x = "ProteinName", by.y = "uniprot_iso")
#   
#   df.fasta.ptm = df.fasta.ptm[
#     which(nchar(df.fasta.ptm$join_PeptideSequence) > min_len_peptide & str_count(
#       df.fasta.ptm$sequence, df.fasta.ptm$join_PeptideSequence) == 1),]
#   
#   start = sapply(seq_len(nrow(df.fasta.ptm)),
#                  function(i) gregexpr(df.fasta.ptm$PeptideSequence[i], 
#                                       df.fasta.ptm$sequence[i])[[1]])
#   mod_loc = sapply(df.fasta.ptm$PeptideSequence, 
#                    function(x) {gregexpr("\\*", x)})
#   
#   # mod_index = sapply(seq_along(mod_loc), function(i){mod_loc[i] = list(
#   #     as.integer(unlist(mod_loc[i])) + start[i][[1]][1])})
#   
#   peptide_mod = mapply(spectro_get_sites, mod_loc, start,
#                        df.fasta.ptm$join_PeptideSequence)
#   
#   df.fasta.ptm$Site = peptide_mod
#   
#   df.fasta.join = unique(df.fasta.ptm[, c("ProteinName", "PeptideSequence", 
#                                           "Site")])
#   
#   if (!removeNonUniqueProteins){
#     add_non_unique = unique(locate_mod_df[grepl(";", locate_mod_df$ProteinName), 
#                                           c("ProteinName", "PeptideSequence")])
#     add_non_unique$Site = add_non_unique$PeptideSequence
#     df.fasta.join = rbindlist(list(df.fasta.join, add_non_unique))
#   }
#   
#   #Data formatting for MSstatsPTM analysis
#   MSstats_PTM = merge(df.ptm, df.fasta.join,
#                       by = c("ProteinName", "PeptideSequence"))
#   MSstats_PTM$ProteinName = paste(MSstats_PTM$ProteinName,
#                                   MSstats_PTM$Site, sep = '_')
#   
#   MSstats_PTM$Intensity = ifelse(MSstats_PTM$Intensity <= 1, NA, 
#                                  MSstats_PTM$Intensity)
#   MSstats_PTM[, join_PeptideSequence := NULL]
#   
#   ## Check if ptm data is empty (indication that mod name is wrong)
#   if (nrow(MSstats_PTM) == 0){
#     stop("The PTM data.table is empty. Please check that the modificationLabel \
#          parameter is set correctly")
#   }
#   
#   if (!is.null(Protein.data)) {
#     if (removeNonUniqueProteins){
#       df.protein = df.protein[!grepl(";", df.protein$ProteinName),]
#     }
#     # df.protein$PeptideSequence = str_extract(df.protein$PeptideSequence,
#     #                                       "([ACDEFGHIKLMNPQRSTVWY]+)")
#     df.protein = df.protein[nchar(PeptideSequence) > min_len_peptide]
#     # 
#     # df.protein$Intensity = ifelse(df.protein$Intensity <= 1, NA, 
#     #                                df.protein$Intensity)
#     
#     MSstats_Protein = df.protein
#   } else {
#     MSstats_Protein = NULL
#   }
#   
#   PTMExpt = list(
#     PTM = MSstats_PTM,
#     PROTEIN = MSstats_Protein
#   )
#   
#   return(PTMExpt)
#   
# }