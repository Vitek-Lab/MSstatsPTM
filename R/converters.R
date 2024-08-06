#' Convert the output of DIA-NN PSM file into MSstatsPTM format
#' 
#' Takes as input the `report.tsv` file from DIA-NN and converts it into 
#' MSstatsPTM format. Requires PSM and an annotation file. Optionally an 
#' additional `report.tsv` file for a corresponding global profiling run can 
#' be included.
#' 
#' @importFrom data.table as.data.table
#' @importFrom MSstatsConvert DIANNtoMSstatsFormat
#' 
#' @param input data.frame of `report.tsv` file produced by Philosopher
#' @param annotation annotation with Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition columns or a path to file. Refer to the example 'annotation' for the meaning of each column.
#' @param input_protein same as `input` for global profiling run. Default is NULL.
#' @param annotation_protein same as `annotation` for global profiling run. Default is NULL.
#' @param fasta_path A string of path to a FASTA file, used to match PTM peptides.
#' @param use_unmod_peptides Boolean if the unmodified peptides in the input 
#' file should be used to construct the unmodified protein output. Only used if
#' `input_protein` is not provided. Default is `FALSE`.
#' @param protein_id_col Use 'Protein.Groups'(default) column for protein name. 
#' @param fasta_protein_name Name of column that matches with the protein names 
#' in `protein_id_col`. The protein names in these two columns must match in 
#' order to join the FASTA file with the DIA-NN output.
#' @param global_qvalue_cutoff The global qvalue cutoff. Default is 0.01.
#' @param qvalue_cutoff local qvalue cutoff for library. Default is 0.01.
#' @param pg_qvalue_cutoff local qvalue cutoff for protein groups Run should be 
#' the same as filename. Default is 0.01.
#' @param useUniquePeptide logical, if TRUE (default) removes peptides that are assigned for more than one proteins. 
#' We assume to use unique peptide for each protein.
#' @param removeFewMeasurements TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param removeOxidationMpeptides TRUE (default) will remove the peptides including oxidation (M) sequence.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param MBR If analaysis was done with match between runs or not. Default is TRUE.
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
#' @return `list` of one or two `data.frame` of class `MSstatsTMT`, named `PTM` and `PROTEIN`
#' 
#' @export
#' 
#' @examples
#' # ptm = read.csv("Phospho/report.tsv", sep="\t")
#' # protein = read.csv("Protein/report.tsv", sep="\t")
#' # annotation = read.csv("Phospho/annotation.csv")
#' # annotation_protein = read.csv("Protein/annotation.csv")
#' 
#' #DIANNtoMSstatsPTMFormat(ptm, annotation, 
#' #                        protein, annotation_protein,
#' #                        fasta_path="fasta_file.fasta")
#' 
DIANNtoMSstatsPTMFormat = function(input,
                                   annotation,
                                   input_protein=NULL,
                                   annotation_protein=NULL,
                                   fasta_path=NULL,
                                   use_unmod_peptides=FALSE,
                                   protein_id_col = "Protein.Group",
                                   fasta_protein_name="uniprot_ac",
                                   global_qvalue_cutoff = 0.01,
                                   qvalue_cutoff = 0.01,
                                   pg_qvalue_cutoff = 0.01,
                                   useUniquePeptide = TRUE,
                                   removeFewMeasurements = TRUE,
                                   removeOxidationMpeptides = TRUE,
                                   removeProtein_with1Feature = FALSE,
                                   MBR=TRUE,
                                   use_log_file = TRUE,
                                   append = FALSE,
                                   verbose = TRUE,
                                   log_file_path = NULL){
  
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path, 
                                      base = "MSstatsPTM_converter_log_")
  
  ## Check input parameters
  checkmate::assertTRUE(!is.null(input) & !is.null(annotation))
  .checkAnnotation(annotation, "LF")
  if (!is.null(annotation_protein)){
    .checkAnnotation(annotation_protein, "LF")
  }
  
  input = as.data.table(input)
  fasta = MSstatsPTM::tidyFasta(fasta_path)
  
  input = MSstatsPTMSiteLocator(input, 
                                protein_name_col=protein_id_col, 
                                unmod_pep_col="Stripped.Sequence",
                                mod_pep_col="Modified.Sequence",
                                clean_mod=TRUE,
                                fasta_file=fasta,
                                fasta_protein_name=fasta_protein_name,
                                mod_id="UniMod",
                                replace_text=TRUE)
  
  if (use_unmod_peptides){
    input_protein = input[input[,..protein_id_col][[1]]  == input$ProteinNameUnmod]  
    annotation_protein = annotation
  } else {
    input = input[input[,..protein_id_col][[1]] != input$ProteinNameUnmod]
  }
  
  ptm_input = DIANNtoMSstatsFormat(input, 
                                   annotation,
                                   global_qvalue_cutoff,
                                   qvalue_cutoff,
                                   pg_qvalue_cutoff,
                                   useUniquePeptide,
                                   removeFewMeasurements,
                                   removeOxidationMpeptides,
                                   removeProtein_with1Feature,
                                   use_log_file,
                                   append,
                                   verbose,
                                   log_file_path,
                                   MBR)
  
  msstats_format = list(PTM=ptm_input, PROTEIN=NULL)
  
  if (!is.null(input_protein)){
    checkmate::assertTRUE(!is.null(input_protein) & 
                            !is.null(annotation_protein))
    
    protein_input = DIANNtoMSstatsFormat(input_protein, 
                                         annotation_protein,
                                         global_qvalue_cutoff,
                                         qvalue_cutoff,
                                         pg_qvalue_cutoff,
                                         useUniquePeptide,
                                         removeFewMeasurements,
                                         removeOxidationMpeptides,
                                         removeProtein_with1Feature,
                                         use_log_file,
                                         append,
                                         verbose,
                                         log_file_path,
                                         MBR)
    
    msstats_format = list(PTM=ptm_input, PROTEIN=protein_input)
    
  }
  return(msstats_format)
}



#' Convert output of TMT labeled Fragpipe data into MSstatsPTM format.
#' 
#' Takes as input TMT experiments which are the output of Fragpipe and converts
#' into MSstatsPTM format. Requires `msstats.csv` file and an annotation file. 
#' Optionally an additional `msstats.csv` file can be uploaded if a 
#' corresponding global profiling run was performed. Site localization is 
#' performed and only high probability localizations are kept.
#' 
#' @export
#' @importFrom data.table as.data.table 
#' @importFrom MSstatsConvert FragPipetoMSstatsFormat
#' @importFrom MSstatsTMT PhilosophertoMSstatsTMTFormat
#' 
#' @param input data.frame of `msstats.csv` file produced by Philosopher
#' @param annotation annotation with Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition columns or a path to file. Refer to the example 'annotation' for the meaning of each column. Channel column should be 
#' consistent with the channel columns (Ignore the prefix "Channel ") in msstats.csv file. Run column should be consistent with the Spectrum.File columns in msstats.csv file.
#' @param input_protein same as `input` for global profiling run. Default is NULL.
#' @param annotation_protein same as `annotation` for global profiling run. Default is NULL.
#' @param use_unmod_peptides Boolean if the unmodified peptides in the input 
#' file should be used to construct the unmodified protein output. Only used if
#' `input_protein` is not provided. Default is `FALSE`.
#' @param label_type Type of labeling used for experiment. Must be one of "LF" 
#' or "TMT". Default is "TMT".
#' @param protein_id_col Use 'Protein'(default) column for TMT. This needs to be
#' changed to "ProteinName" for label free. For TMT, 'Master.Protein.Accessions'
#' can be used instead to get the protein ID with single protein.
#' @param peptide_id_col Use 'Peptide.Sequence'(default) column for TMT. Must be
#' changed to "PeptideSequence" for label free. "Modified.Peptide.Sequence" can 
#' be used instead to get the modified peptide sequence.
#' @param mod_id_col Column containing the modified Amino Acids. For example, a Phosphorylation experiment may pass `STY`. The corresponding column with `STY` combined with the mass (e.x. `STY.79.9663`) will be selected. Default is `STY`.
#' @param localization_cutoff Minimum localization score required to keep modification. Default is .75.
#' @param remove_unlocalized_peptides Boolean indicating if peptides without all sites localized should be kept. Default is TRUE (non-localized sites will be removed).
#' @param Purity_cutoff Cutoff for purity. Default is 0.6. Purity refers to how much of the detected ion signal 
#' within a specific inclusion window belongs to the target molecule or its closely related forms, 
#' compared to any other unwanted signals or noise. Higher values indicate greater purity.
#' @param PeptideProphet_prob_cutoff Cutoff for the peptide identification probability. Default is 0.7. 
#' The probability is confidence score determined by PeptideProphet and higher values indicate greater confidence.
#' @param useUniquePeptide logical, if TRUE (default) removes peptides that are assigned for more than one proteins. 
#' We assume to use unique peptide for each protein.
#' @param rmPSM_withfewMea_withinRun TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmPeptide_OxidationM TRUE (default) will remove the peptides including oxidation (M) sequence.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum (default) or max - when there are multiple measurements for certain feature in certain run, 
#' select the feature with the largest summation or maximal value.
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
#' @return `list` of one or two `data.frame` of class `MSstatsTMT`, named `PTM` and `PROTEIN`
#' 
#' @export
#' 
#' @examples 
#' # TMT Example (with global profiling run)
#' head(fragpipe_input)
#' head(fragpipe_annotation)
#' head(fragpipe_input_protein)
#' head(fragpipe_annotation_protein)
#' 
#' msstats_data = FragPipetoMSstatsPTMFormat(fragpipe_input,
#'                                           fragpipe_annotation,
#'                                           fragpipe_input_protein, 
#'                                           fragpipe_annotation_protein,
#'                                           label_type="TMT",
#'                                           mod_id_col = "STY",
#'                                           localization_cutoff=.75,
#'                                           remove_unlocalized_peptides=TRUE)
#' head(msstats_data$PTM)
#' head(msstats_data$PROTEIN)
#' 
#' # LFQ Example (w/out global profiling run)
#' input = system.file("tinytest/raw_data/Fragpipe/MSstats.csv", 
#'                                         package = "MSstatsPTM")
#' input = data.table::fread(input)
#' annot = system.file("tinytest/raw_data/Fragpipe/experiment_annotation.tsv", 
#'                                         package = "MSstatsPTM")
#' annot = data.table::fread(annot)                                        
#' 
#' msstats_data = FragPipetoMSstatsPTMFormat(input,
#'                                           annot,
#'                                           label_type="LF",
#'                                           mod_id_col = "STY",
#'                                           localization_cutoff=.75,
#'                                           protein_id_col = "ProteinName",
#'                                           peptide_id_col = "PeptideSequence")
#' head(msstats_data$PTM)
#' 
#' # Note that this is NULL because we did not include a global profiling run.
#' # Ideally, you should include an independent global profiling run.
#' head(msstats_data$PROTEIN)
#' 
FragPipetoMSstatsPTMFormat = function(input,
                                       annotation=NULL,
                                       input_protein=NULL,
                                       annotation_protein=NULL,
                                       use_unmod_peptides=FALSE,
                                       label_type="TMT",
                                       protein_id_col = "Protein",
                                       peptide_id_col = "Peptide.Sequence",
                                       mod_id_col = "STY",
                                       localization_cutoff=.75,
                                       remove_unlocalized_peptides=TRUE,
                                       Purity_cutoff = 0.6,
                                       PeptideProphet_prob_cutoff = 0.7,
                                       useUniquePeptide = TRUE,
                                       rmPSM_withfewMea_withinRun = FALSE,
                                       rmPeptide_OxidationM = TRUE,
                                       rmProtein_with1Feature = FALSE,
                                       summaryforMultipleRows = sum,
                                       use_log_file = TRUE,
                                       append = FALSE,
                                       verbose = TRUE,
                                       log_file_path = NULL){
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path, 
                                      base = "MSstatsPTM_converter_log_")
  
  ## Check input parameters
  checkmate::assertTRUE(!is.null(input))
  if (label_type == "TMT"){
    checkmate::assertTRUE(!is.null(annotation))
    .checkAnnotation(annotation, "TMT")
    if (!is.null(annotation_protein)){
      .checkAnnotation(annotation_protein, "TMT")
    }
  }
  
  input = as.data.table(input)
  
  mod_id_col = .getFullModID(input, mod_id_col)
  input$Start = input$Protein.Start
  if ("Is.Unique" %in% colnames(input)){
    input$Is.Unique = as.logical(input$Is.Unique)
  }
  
  input = MSstatsPTMSiteLocator(input, 
                       protein_name_col= protein_id_col,
                       unmod_pep_col = peptide_id_col,
                       mod_pep_col = mod_id_col,
                       clean_mod=FALSE,
                       fasta_file=NULL, 
                       fasta_protein_name="header",
                       mod_id="\\*", 
                       localization_scores=TRUE,
                       localization_cutoff=localization_cutoff,
                       remove_unlocalized_peptides=remove_unlocalized_peptides,
                       terminus_included=FALSE, 
                       terminus_id="\\.")
  
  if (use_unmod_peptides){
    input_protein = input[input[[protein_id_col]] == input$ProteinNameUnmod]
    input = input[input[[protein_id_col]] != input$ProteinNameUnmod]  
    annotation_protein = annotation
  } else {
    input = input[input[[protein_id_col]] != input$ProteinNameUnmod]  
  }
  
  input[[peptide_id_col]] = input[[colnames(input)[grepl("new_peptide_col", colnames(input))][[1]]]]
  if (label_type == "TMT"){
    
    ptm_input = PhilosophertoMSstatsTMTFormat(
      input, annotation, protein_id_col, peptide_id_col, Purity_cutoff, 
      PeptideProphet_prob_cutoff, useUniquePeptide, rmPSM_withfewMea_withinRun, 
      rmPeptide_OxidationM, rmProtein_with1Feature, summaryforMultipleRows, 
      use_log_file, append, verbose, log_file_path)
  } else if (label_type == "LF"){
    msstats_cols = c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                     "FragmentIon", "ProductCharge", "IsotopeLabelType", 
                     "Condition", "BioReplicate", "Run", "Intensity")
    ptm_input = FragPipetoMSstatsFormat(
      input[,..msstats_cols], useUniquePeptide, rmPSM_withfewMea_withinRun, 
      rmProtein_with1Feature, summaryforMultipleRows,
      use_log_file, append, verbose, log_file_path
    )
  }
  
  msstats_format = list(PTM=ptm_input, PROTEIN=NULL)
  
  if (!is.null(input_protein)){
    
    if (label_type == "TMT"){
      protein_input = PhilosophertoMSstatsTMTFormat(
        input_protein, annotation_protein, protein_id_col, peptide_id_col, 
        Purity_cutoff, PeptideProphet_prob_cutoff, useUniquePeptide, 
        rmPSM_withfewMea_withinRun, rmPeptide_OxidationM, rmProtein_with1Feature, 
        summaryforMultipleRows, use_log_file, append, verbose, log_file_path)
    } else if (label_type == "LF"){
      protein_input = FragPipetoMSstatsFormat(
        input_protein[,..msstats_cols], useUniquePeptide, rmPSM_withfewMea_withinRun, 
        rmProtein_with1Feature, summaryforMultipleRows,
        use_log_file, append, verbose, log_file_path)
    }
    msstats_format = list(PTM=ptm_input, PROTEIN=protein_input)
    
  }
  
  return(msstats_format)
}


#' Convert output of label-free or TMT MaxQuant experiments into MSstatsPTM format
#' 
#' Takes as input LF/TMT experiments from MaxQ and converts the data into the 
#' format needed for MSstatsPTM. Requires modified evidence.txt file from MaxQ 
#' and an annotation file for PTM data. To adjust modified peptides for changes 
#' in global protein level, unmodified TMT experimental data must also be 
#' returned. Optionally can use `Phospho(STY)Sites.txt` (or other PTM specific 
#' files) from MaxQuant, but this is not recommended. If PTM specific file 
#' provided, the raw intensities must be provided, not a ratio.
#'
#' @export
#' @importFrom stringr str_extract regex str_replace fixed str_split
#' @importFrom data.table melt as.data.table `:=` `%like%` setcolorder setDT tstrsplit setcolorder
#' @importFrom MSstatsConvert MaxQtoMSstatsFormat
#' @importFrom MSstatsTMT MaxQtoMSstatsTMTFormat
#' @importFrom checkmate assertChoice assertLogical
#' 
#' @param evidence name of 'evidence.txt' data, which includes feature-level 
#' data for enriched (PTM) data.
#' @param annotation data frame annotation file for the ptm level data.
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
#' @param annotation_protein data frame annotation file for the protein level data.
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
#' @param remove_other_mods Remove peptides which include modfications other 
#' than the one listed in `mod_id`. Default is `TRUE`. For example, in an 
#' experiment targeting Phosphorylation, setting this parameter to `TRUE` would 
#' remove peptides like 
#' (Acetyl (Protein N-term))AAAAPDSRVS(Phospho (STY))EEENLK. Set this parameter 
#' to `FALSE` to keep peptides with extraneous modifications.
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
#' # TMT experiment
#' head(maxq_tmt_evidence)
#' head(maxq_tmt_annotation)
#' 
#' msstats_format_tmt = MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
#'                         annotation=maxq_tmt_annotation,
#'                         fasta=system.file("extdata", "maxq_tmt_fasta.fasta", package="MSstatsPTM"),
#'                         fasta_protein_name="uniprot_ac",
#'                         mod_id="\\(Phospho \\(STY\\)\\)",
#'                         use_unmod_peptides=TRUE,
#'                         labeling_type = "TMT",
#'                         which_proteinid_ptm = "Proteins")
#' 
#' head(msstats_format_tmt$PTM)
#' head(msstats_format_tmt$PROTEIN)
#' 
#' # LF experiment
#' head(maxq_lf_evidence)
#' head(maxq_lf_annotation)
#' 
#' msstats_format_lf = MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
#'                         annotation=maxq_lf_annotation,
#'                         fasta=system.file("extdata", "maxq_lf_fasta.fasta", package="MSstatsPTM"),
#'                         fasta_protein_name="uniprot_ac",
#'                         mod_id="\\(Phospho \\(STY\\)\\)",
#'                         use_unmod_peptides=TRUE,
#'                         labeling_type = "LF",
#'                         which_proteinid_ptm = "Proteins")
#' head(msstats_format_lf$PTM)
#' head(msstats_format_lf$PROTEIN)
MaxQtoMSstatsPTMFormat = function(evidence=NULL,
                                  annotation=NULL,
                                  fasta_path,
                                  fasta_protein_name="uniprot_ac",
                                  mod_id="\\(Phospho \\(STY\\)\\)",
                                  sites_data=NULL,
                                  evidence_prot = NULL,
                                  proteinGroups = NULL,
                                  annotation_protein = NULL,
                                  use_unmod_peptides=FALSE,
                                  labeling_type = "LF",
                                  mod_num = 'Single',
                                  TMT_keyword = "TMT",
                                  ptm_keyword = "phos",
                                  which_proteinid_ptm = "Proteins",
                                  which_proteinid_protein = "Proteins",
                                  remove_other_mods=TRUE,
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
  
  .checkAnnotation(annotation, labeling_type)
  if (!is.null(annotation_protein)){
    .checkAnnotation(annotation_protein, labeling_type)
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
  
  annot.ptm = as.data.table(annotation)
  
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
                                             terminus_included = FALSE,
                                             remove_underscores=TRUE,
                                             remove_other_mods=remove_other_mods,
                                             bracket="(",
                                             replace_text=TRUE)
      
      msstatsptm_input = MaxQtoMSstatsTMTFormatHelper(evidence_sites,
                                                      annotation,
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
                                             terminus_included=FALSE,
                                             remove_underscores=TRUE,
                                             remove_other_mods=remove_other_mods,
                                             bracket="(",
                                             replace_text=TRUE)
      
      msstatsptm_input = MaxQtoMSstatsFormatHelper(evidence_sites,
                                                   annotation,
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
    annotation_protein = as.data.table(annotation_protein)
    
    ## Clean raw data
    #evidence = as.data.table(evidence)
    # evidence = evidence[!grepl("phos", evidence$Raw.file),]
    #proteinGroups = as.data.table(proteinGroups)
    
    if(labeling_type == "TMT"){
      msstats.abun = MaxQtoMSstatsTMTFormat(evidence = evidence_prot,
                                            proteinGroups = proteinGroups,
                                            annotation = annotation_protein,
                                            which.proteinid = which_proteinid_protein)
    } else if (labeling_type == "LF"){
      msstats.abun = MaxQtoMSstatsFormat(evidence = evidence_prot,
                                         proteinGroups = proteinGroups,
                                         annotation = annotation_protein,
                                         proteinID = which_proteinid_protein)
    }
    
    MSstatsPTMformat = list('PTM' = msstatsptm_input, 
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

#' Convert Proteome Discoverer output into MSstatsPTM format
#' 
#' Import Proteome Discoverer files, identify modification site location.
#' 
#' @param input PD report corresponding with enriched experimental data.
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which 
#' includes Condition, BioReplicate, Run information. 'Run' will be matched 
#' with 'Spectrum.File'
#' @param fasta_path string containing path to the corresponding fasta file for 
#' the modified peptide dataset.
#' @param protein_input PD report corresponding with unmodified experimental 
#' data.
#' @param annotation_protein Same format as `annotation` corresponding to 
#' unmodified data.
#' @param labeling_type type of experimental design, must be one of `LF` for 
#' label free or `TMT` for tandem mass tag.
#' @param mod_id Character that indicates the modification of interest. Default 
#' is `\\(Phospho\\)`. Note `\\` must be included before special characters.
#' @param keep_all_mods Boolean indicating whether to keep or remove peptides 
#' not in `mod_id`. Default is FALSE.
#' @param use_localization_cutoff Boolean indicating whether to use a custom 
#' localization cutoff or rely on PD's modifications column. `TRUE` is default 
#' and apply custom cutoff `localization_cutoff`.
#' @param use_unmod_peptides If `protein_input` is not provided, 
#' unmodified peptides can be extracted from `input` to be used in place of a 
#' global profiling run. Default is `FALSE`.
#' @param fasta_protein_name Name of fasta column that matches with protein name
#' in evidence file. Default is `uniprot_iso`.
#' @param localization_cutoff Minimum localization score required to keep modification. Default is .75.
#' @param remove_unlocalized_peptides Boolean indicating if peptides without all sites localized should be kept. Default is TRUE (non-localized sites will be removed).
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
#' @importFrom MSstatsConvert PDtoMSstatsFormat
#' @importFrom MSstatsTMT PDtoMSstatsTMTFormat
#' @importFrom stringr str_split str_trim str_split_i
#' @importFrom data.table setnames
#' @importFrom stringi stri_sub_replace_all
#' @return `list` of `data.table`
#' @export 
#' 
#' @examples
#' # Global profiling example
#' input = system.file("tinytest/raw_data/PD/pd-ptm-input.csv", 
#' package = "MSstatsPTM")
#' input = data.table::fread(input)
#' annot = system.file("tinytest/raw_data/PD/pd-ptm-annot.csv",
#'                     package = "MSstatsPTM")
#' annot = data.table::fread(annot)
#' input_protein = system.file("tinytest/raw_data/PD/protein-input.csv",
#'                             package = "MSstatsPTM")
#' input_protein = data.table::fread(input_protein)
#' annot_protein = system.file("tinytest/raw_data/PD/protein-annot.csv",
#'                             package = "MSstatsPTM")
#' annot_protein = data.table::fread(annot_protein)
#' fasta_path=system.file("extdata", "pd_with_proteome.fasta", 
#'                        package="MSstatsPTM")
#' pd_imported = PDtoMSstatsPTMFormat(
#'     input,
#'     annotation = annot,
#'     protein_input = input_protein,
#'     annotation_protein = annot_protein,
#'     fasta_path = fasta_path,
#'     mod_id = "\\(GG\\)",
#'     labeling_type = "TMT",
#'     use_localization_cutoff = FALSE,
#'     which_proteinid = "Master.Protein.Accessions")
#'    
#' head(pd_imported$PTM)
#' head(pd_imported$PROTEIN)
#' 
#' # No global profiling example
#' head(pd_psm_input)
#' head(pd_annotation)
#' 
#' msstats_format = PDtoMSstatsPTMFormat(pd_psm_input, 
#'                                       pd_annotation, 
#'                                       system.file("extdata", "pd_fasta.fasta", package="MSstatsPTM"),
#'                                       use_unmod_peptides=TRUE, 
#'                                       which_proteinid = "Master.Protein.Accessions")
#' 
#' head(msstats_format$PTM)
#' head(msstats_format$PROTEIN)
PDtoMSstatsPTMFormat = function(input,
                                annotation,
                                fasta_path,
                                protein_input=NULL,
                                annotation_protein=NULL,
                                labeling_type = "LF",
                                mod_id="\\(Phospho\\)",
                                use_localization_cutoff=FALSE,
                                keep_all_mods=FALSE,
                                use_unmod_peptides=FALSE,
                                fasta_protein_name="uniprot_iso",
                                localization_cutoff=75,
                                remove_unlocalized_peptides=TRUE,
                                useNumProteinsColumn = FALSE,
                                useUniquePeptide = TRUE,
                                summaryforMultipleRows = max,
                                removeFewMeasurements = TRUE,
                                removeOxidationMpeptides = FALSE,
                                removeProtein_with1Peptide = FALSE,
                                which_quantification = "Precursor.Area",
                                which_proteinid = "Protein.Group.Accessions",
                                use_log_file = TRUE,
                                append = FALSE,
                                verbose = TRUE,
                                log_file_path = NULL){
  
  input = as.data.table(input)
  
  if (is.null(fasta_path) & use_localization_cutoff == TRUE){
    stop("A FASTA file must be included if using a custom localization cutoff. Please pass a FASTA file to `fasta_path`.")
  }
  
  if (!is.null(protein_input) & use_unmod_peptides == TRUE){
    stop("Either pass protein_input data or set use_unmod_peptides = TRUE, not both")
  }
  
  .checkAnnotation(annotation, labeling_type)
  if (!is.null(annotation_protein)){
    .checkAnnotation(annotation_protein, labeling_type)
  }
  
  if (labeling_type == "LF"){
    sequence_col = "Sequence"
  } else if (labeling_type == "TMT"){
    sequence_col = "Annotated.Sequence"
    input[[sequence_col]] = toupper(str_split_i(input[[sequence_col]], 
                                                "\\.", i=2))
  }
  
  if (use_localization_cutoff){
    input = .getPDmods(input)
    
    input = MSstatsPTMSiteLocator(input, 
                                  protein_name_col= which_proteinid,
                                  unmod_pep_col = sequence_col,
                                  mod_pep_col = "ModSequence",
                                  clean_mod=FALSE,
                                  fasta_file=fasta_path, 
                                  fasta_protein_name=fasta_protein_name,
                                  mod_id="\\*", 
                                  localization_scores=TRUE,
                                  localization_cutoff=localization_cutoff,
                                  remove_unlocalized_peptides=remove_unlocalized_peptides,
                                  terminus_included=FALSE, 
                                  terminus_id="\\.")
  } else {
    input = .extract_pd_mods(input, mod_id, keep_all_mods, sequence_col)
    # input[,which_proteinid] = paste(input[,..which_proteinid][[1]], mods,sep="_")
    
    input = MSstatsPTMSiteLocator(input, 
                                  protein_name_col= which_proteinid,
                                  unmod_pep_col = sequence_col,
                                  mod_pep_col = "ModSequence",
                                  clean_mod=FALSE,
                                  fasta_file=fasta_path, 
                                  fasta_protein_name=fasta_protein_name,
                                  mod_id="\\*", 
                                  localization_scores=FALSE,
                                  localization_cutoff=localization_cutoff,
                                  remove_unlocalized_peptides=remove_unlocalized_peptides,
                                  terminus_included=FALSE, 
                                  terminus_id="\\.")
    
  }
  
  input[[sequence_col]] = input[, "ModSequence"]
  
  if (labeling_type == "LF"){
    ptm_input = PDtoMSstatsFormat(input, annotation, useNumProteinsColumn,
                                  useUniquePeptide, summaryforMultipleRows,
                                  removeFewMeasurements, removeOxidationMpeptides,
                                  removeProtein_with1Peptide, 
                                  which_quantification, which_proteinid, 
                                  sequence_col, use_log_file, append, verbose,
                                  log_file_path)
  } else if (labeling_type == "TMT"){
    ptm_input = PDtoMSstatsTMTFormat(input, annotation, which_proteinid, 
                                     useNumProteinsColumn,
                                     useUniquePeptide, TRUE,
                                     removeProtein_with1Peptide, 
                                     summaryforMultipleRows, 
                                     use_log_file, append, verbose,
                                     log_file_path)
  }
  
  if ("PeptideModifiedSequence" %in% colnames(ptm_input)){
    setnames(ptm_input, c("PeptideModifiedSequence"), c("PeptideSequence"))
  }
  
  if (!is.null(protein_input)) {
    if (labeling_type == "LF"){
      protein_input = PDtoMSstatsFormat(protein_input, annotation_protein, 
                                        useNumProteinsColumn,
                                        useUniquePeptide, summaryforMultipleRows,
                                        removeFewMeasurements, 
                                        removeOxidationMpeptides,
                                        removeProtein_with1Peptide, 
                                        which_quantification, which_proteinid, 
                                        "Sequence", use_log_file, append, 
                                        verbose, log_file_path)
    } else if (labeling_type == "TMT"){
      protein_input = PDtoMSstatsTMTFormat(protein_input, annotation_protein, 
                                       which_proteinid, useNumProteinsColumn,
                                       useUniquePeptide, TRUE,
                                       removeProtein_with1Peptide, 
                                       summaryforMultipleRows, 
                                       use_log_file, append, verbose,
                                       log_file_path)
    }
    
    if ("PeptideModifiedSequence" %in% colnames(protein_input)){
      setnames(protein_input, c("PeptideModifiedSequence"), 
               c("PeptideSequence"))
    }
      
    ptm_input = ptm_input[grepl("\\*", ptm_input[, "PeptideSequence"]), ]
    
    msstats_input = list(PTM = ptm_input, PROTEIN = protein_input)
  } else {
    if (use_unmod_peptides){
      ptm_input=as.data.frame(ptm_input)
      protein_input = ptm_input[!grepl("\\*", ptm_input$PeptideSequence),]
      ptm_input = ptm_input[grepl("\\*", ptm_input$PeptideSequence),]
      
      msstats_input = list(PTM = ptm_input, PROTEIN = protein_input)
    } else {
      ptm_input = ptm_input[grepl("\\*", ptm_input[, "PeptideSequence"]), ]
      
      msstats_input = list(PTM = ptm_input)
    }
  }
  
  return(msstats_input)
  
}

#' Converts non-TMT Progenesis output into the format needed for MSstatsPTM
#'
#' @export
#' @importFrom MSstatsConvert ProgenesistoMSstatsFormat
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
#' input = system.file("tinytest/raw_data/Progenesis/progenesis_peptide.csv", 
#'                                 package = "MSstatsPTM")
#' input = data.table::fread(input)
#' colnames(input) = unlist(input[1,])
#' input = input[-1,]
#' annot = system.file("tinytest/raw_data/Progenesis/phospho_annotation.csv", 
#'                                 package = "MSstatsPTM")
#' annot = data.table::fread(annot)
#' prog_imported = ProgenesistoMSstatsPTMFormat(
#'     input, 
#'     annot
#' )
#' head(prog_imported$PTM)
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
  
  .checkAnnotation(annotation, "LF")
  if (!is.null(annotation_protein)){
    .checkAnnotation(annotation_protein, "LF")
  }
  
  if (!"Source.File" %in% colnames(input)){
    message("Pivoting input data..")
    input = .pivotPS(input)
  }
  
  message("Merging with annotation..")
  input = merge(input, annotation, all.x = TRUE, by = "Raw.File")
  
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
#' @param fasta_protein_name Name of fasta column that matches with protein name
#' in evidence file. Default is `uniprot_iso`.
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
#' @importFrom MSstatsConvert SkylinetoMSstatsFormat
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
                                     fasta_protein_name="uniprot_iso",
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
    .checkAnnotation(annotation, "LF")
    if (!is.null(annotation_protein)){
      .checkAnnotation(annotation_prot, "LF")
    }
    
    input = merge(input, annotation, all.x = TRUE, by = "File Name")
  }
  
  ## Filter for required modifications
  message("Filtering modifications based on function arguements..")
  if (is.null(input_protein) & use_unmod_peptides){
    input_protein = input[!grepl("\\[+", input$`Peptide Modified Sequence`), ]
    
    if (nrow(input_protein)){
      stop("No unmodified peptides in PTM dataset. It cannot be used for \
           adjustment.")
    }
  }
  
  setnames(input, c("Protein.Name", "Peptide", "Peptide.Modified.Sequence"), 
           c("ProteinName","PeptideSequence", "PeptideModifiedSequence"))
  
  ## Ensure only modified peptides are used
  input = input[grepl("\\[+", input$`PeptideModifiedSequence`), ]
  
  ## Locate sites
  input = MSstatsPTMSiteLocator(input, 
                                fasta_file=fasta_path,
                                fasta_protein_name=fasta_protein_name,
                                terminus_included=FALSE,
                                mod_id_is_numeric=TRUE)

  
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


#' Convert Spectronaut output into MSstatsPTM format
#' 
#' Converters label-free Spectronaut data into MSstatsPTM format. Requires PSM 
#' output from Spectronaut and a custom made annotation file, mapping the run 
#' name to the condition and bioreplicate. Can optionally take a seperate PSM 
#' file for a global profiling run. If no global profiling run provided, the 
#' function can extract the unmodified peptides from the PTM PSM file and use 
#' them as a global profiling run (not recommended).
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
#' @param annotation_protein name of annotation file for global protein data, in
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
#' @param remove_other_mods Remove peptides which include modfications other 
#' than the one listed in `mod_id`. Default is `TRUE`. For example, in an 
#' experiment targeting Phosphorylation, setting this parameter to `TRUE` would 
#' remove peptides like 
#' (Acetyl (Protein N-term))AAAAPDSRVS(Phospho (STY))EEENLK. Set this parameter 
#' to `FALSE` to keep peptides with extraneous modifications.
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
#' @importFrom MSstatsConvert SpectronauttoMSstatsFormat
#' @examples
#' 
#' head(spectronaut_input)
#' head(spectronaut_annotation)
#' 
#' msstats_input = SpectronauttoMSstatsPTMFormat(spectronaut_input, 
#'                   annotation=spectronaut_annotation, 
#'                   fasta_path=system.file("extdata", "spectronaut_fasta.fasta", package="MSstatsPTM"),
#'                   use_unmod_peptides=TRUE,
#'                   mod_id = "\\[Phospho \\(STY\\)\\]",
#'                   fasta_protein_name = "uniprot_iso"
#'                   )
#' 
#' head(msstats_input$PTM)
#' head(msstats_input$PROTEIN)
SpectronauttoMSstatsPTMFormat = function(
    input,
    annotation = NULL,
    fasta_path = NULL,
    protein_input = NULL,
    annotation_protein = NULL,
    use_unmod_peptides=FALSE,
    intensity = "PeakArea",
    mod_id="\\[Phospho \\(STY\\)\\]",
    fasta_protein_name="uniprot_iso",
    remove_other_mods=TRUE,
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
                                fasta_file=fasta_path, 
                                fasta_protein_name=fasta_protein_name,
                                mod_id=mod_id, 
                                terminus_included=FALSE, terminus_id="\\.",
                                remove_underscores=TRUE,
                                remove_other_mods=remove_other_mods,
                                bracket="[",
                                replace_text=TRUE)
  
  ptm_input = SpectronauttoMSstatsFormat(input, annotation, intensity, 
                                         filter_with_Qvalue,
                                         qvalue_cutoff, useUniquePeptide, 
                                         removeFewMeasurements, 
                                         removeProtein_with1Feature, 
                                         summaryforMultipleRows)
  
  msstats_input = list(PTM = ptm_input)
  if (!is.null(protein_input)) {
    protein_input = SpectronauttoMSstatsFormat(protein_input, 
                                               annotation_protein, intensity, 
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

#' Import Metamorpheus files into PTM format
#' 
#' @author Anthony Wu
#' 
#' @param input name of Metamorpheus output file, which is tabular format. Use the AllQuantifiedPeaks.tsv file from the Metamorpheus output.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate.
#' @param fasta_path string containing path to the corresponding fasta file for 
#' the modified peptide dataset.
#' @param input_protein same as `input` for global profiling run. Default is NULL.
#' @param annotation_protein same as `annotation` for global profiling run. Default is NULL.
#' @param use_unmod_peptides If `protein_input` is not provided, 
#' unmodified peptides can be extracted from `input` to be used in place of a 
#' global profiling run. Default is `FALSE`.
#' @param mod_ids List of modifications of interest. Default 
#' is a list with only `Common Biological:Phosphorylation on S`. 
#' Please note that the 'mod_ids' parameter currently supports lists of size 1 only. 
#' Future updates aim to extend its functionality to accommodate lists of greater sizes.
#' Note `\\` must be included before special characters. 
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
#' @importFrom MSstatsConvert MetamorpheusToMSstatsFormat
#' @export
#' 
#' @examples 
#' input = system.file("tinytest/raw_data/Metamorpheus/AllQuantifiedPeaks.tsv", 
#'                                 package = "MSstatsPTM")
#' input = data.table::fread(input)
#' annot = system.file("tinytest/raw_data/Metamorpheus/ExperimentalDesign.tsv", 
#'                                 package = "MSstatsPTM")
#' annot = data.table::fread(annot)
#' input_protein = system.file("tinytest/raw_data/Metamorpheus/AllQuantifiedPeaksGlobalProteome.tsv",
#'                                 package = "MSstatsPTM")
#' input_protein = data.table::fread(input_protein)
#' annot_protein = system.file("tinytest/raw_data/Metamorpheus/ExperimentalDesignGlobalProteome.tsv", 
#'                                 package = "MSstatsPTM")
#' annot_protein = data.table::fread(annot_protein)
#' fasta_path=system.file("extdata", "metamorpheus_fasta.fasta", 
#'                                 package="MSstatsPTM")
#' metamorpheus_imported = MetamorpheusToMSstatsPTMFormat(
#'     input, 
#'     annot, 
#'     fasta_path=fasta_path,
#'     input_protein=input_protein,
#'     annotation_protein=annot_protein,
#'     use_unmod_peptides=FALSE,
#'     mod_ids = c("\\[Common Fixed:Carbamidomethyl on C\\]")
#' )
#' head(metamorpheus_imported$PTM)
#' head(metamorpheus_imported$PROTEIN)
MetamorpheusToMSstatsPTMFormat = function(input,
                                          annotation,
                                          fasta_path,
                                          input_protein = NULL,
                                          annotation_protein = NULL,
                                          use_unmod_peptides = FALSE,
                                          mod_ids = c("\\[Common Biological:Phosphorylation on S\\]"),
                                          useUniquePeptide = TRUE, 
                                          removeFewMeasurements = TRUE,
                                          removeProtein_with1Feature = FALSE, 
                                          summaryforMultipleRows = max,
                                          use_log_file = TRUE, 
                                          append = FALSE, 
                                          verbose = TRUE, 
                                          log_file_path = NULL) {
    
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path, 
                                        base = "MSstatsPTM_converter_log_")
    
    mod_id = mod_ids[[1]]
    input = as.data.table(input)
    
    ## Check input parameters
    checkmate::assertTRUE(!is.null(input))
    .checkAnnotation(annotation, "LF")
    if (!is.null(input_protein)) {
        checkmate::assertTRUE(!is.null(annotation_protein))
    }
    
    fasta = MSstatsPTM::tidyFasta(fasta_path)

    protein_id_col = "Protein Group"
    input = MSstatsPTMSiteLocator(input,
                                  protein_name_col=protein_id_col,
                                  unmod_pep_col="Base Sequence",
                                  mod_pep_col="Full Sequence",
                                  fasta_file=fasta,
                                  mod_id=mod_id,
                                  fasta_protein_name="uniprot_iso")

    if (use_unmod_peptides){
        input_protein = input[input[,..protein_id_col][[1]]  == input$ProteinNameUnmod]
        annotation_protein = annotation
    } else {
        input = input[input[,..protein_id_col][[1]] != input$ProteinNameUnmod]
    }
    
    ptm_input = MetamorpheusToMSstatsFormat(input,
                                            annotation,
                                            useUniquePeptide, 
                                            removeFewMeasurements,
                                            removeProtein_with1Feature, 
                                            summaryforMultipleRows,
                                            use_log_file, 
                                            append, 
                                            verbose, 
                                            log_file_path)
    
    msstats_format = list(PTM = ptm_input, PROTEIN = NULL)
    
    if (!is.null(input_protein)) {
        protein_input = MetamorpheusToMSstatsFormat(input_protein,
                                                    annotation_protein,
                                                    useUniquePeptide, 
                                                    removeFewMeasurements,
                                                    removeProtein_with1Feature, 
                                                    summaryforMultipleRows,
                                                    use_log_file, 
                                                    append, 
                                                    verbose, 
                                                    log_file_path)
        
        ptm_input = ptm_input[grepl(mod_id, ptm_input$PeptideSequence),]
        
        msstats_format = list(PTM = ptm_input, PROTEIN = protein_input)
    }
    
    return(msstats_format)
    
}