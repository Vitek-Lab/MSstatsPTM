#' Convert Peaks Studio output into MSstatsPTM format
#' 
#' Currently only supports label-free quantification.
#' 
#' @param input name of Peaks Studio PTM output
#' @param annotation name of annotation file which includes Raw.file, Condition,
#' BioReplicate, Run. For example annoation see example below.
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
#' @importFrom data.table as.data.table melt
#' @importFrom MSstatsConvert MSstatsBalancedDesign
#' @return `list` of `data.table`
#' @export 
#' 
#' @examples
#' #TODO: Add examples
PStoMSstatsPTMFormat = function(
    input, annotation, input_protein = NULL, annotation_protein = NULL,
    use_unmod_peptides = FALSE, target_modification = NULL, 
    remove_oxidation_peptides = FALSE, remove_multi_mod_types = FALSE,
    summaryforMultipleRows = max, use_log_file = TRUE, append = FALSE, 
    verbose = TRUE, log_file_path = NULL
){
  
  input = as.data.table(input)
  input_prot = NULL
  message("Pivoting input data..")
  input = .pivotPS(input)
  message("Merging with annotation..")
  input = merge(input, annotation, all.x = TRUE, by = "Raw.File")
  
  if (!is.null(input_protein)){
    message("Converting unmodified protein data..")
    ## TODO: Add unmod protein converter
    stop("Converter currently does not support seperate unmodified protein run.\
         Please contact package maintainer for more information.")
  }
  
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
  
  if (!is.null(input_prot)){
    input_prot[, c("Raw.File", "Mod", "End", "Start"):=NULL]
    input_prot = MSstatsBalancedDesign(input_prot, 
                                       c('PeptideSequence', 'ProductCharge'))
    input_prot = as.data.table(input_prot)[, list(
      Intensity = summaryforMultipleRows(Intensity, na.rm = TRUE)), 
      by = feature_cols]
    
    input_prot = as.data.frame(input_prot)
  }
  
  return(list("PTM" = as.data.frame(input),
              "PROTEIN" = input_prot))
}

#' Convert Skyline output into MSstatsPTM format
#' 
#' Currently only supports label-free quantification.
#' 
#' @param input name of Skyline PTM output
#' @param fasta A string of path to a FASTA file, used to match PTM peptides.
#' @param annotation name of 'annotation.txt' data which includes Condition, 
#' BioReplicate, Run. If annotation is already complete in Skyline, use 
#' annotation=NULL (default). It will use the annotation information from input.
#' @param input_protein name of Skyline unmodified protein output (optional)
#' @param annotation_protein name of 'annotation.txt' data which includes Condition, 
#' BioReplicate, Run for unmodified protein output. This can be the same as 
#' `annotation`.
#' @param use_unmod_peptides Boolean if the unmodified peptides in the input 
#' file should be used to construct the unmodified protein output. Only used if
#' `input_protein` is not provided. Default is `FALSE`
#' @param removeiRT TRUE (default) will remove the proteins or peptides which 
#' are labeld 'iRT' in 'StandardType' column. FALSE will keep them.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that 
#' have greater than qvalue_cutoff in DetectionQValue column. Those intensities 
#' will be replaced with zero and will be considered as censored missing values 
#' for imputation purpose.
#' @param qvalue_cutoff Cutoff for DetectionQValue. default is 0.01.
#' @param removeFewMeasurements TRUE will remove the features that 
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
#' 
#' @return `list` of `data.table`
#' @export
#' 
#' @examples
#' #TODO: Add examples
SkylinetoMSstatsPTMFormat = function(input,
                                     fasta,
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
  input = MSstatsPTMSiteLocator(input, fasta_file=fasta, mod_id_is_numeric=TRUE, 
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