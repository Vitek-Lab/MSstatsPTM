#' Check annotation file to ensure it is formatted correctly
#' @noRd
.checkAnnotation = function(annotation, label_type){
  columns = colnames(annotation)
  
  if (!("Run" %in% columns | "Raw.file" %in% columns)){
    msg = "Run column missing in annotation file. Annotation must include either
    `Run` or `Raw.file` column to match with input data."
    stop(msg)
  }
  
  if (label_type == "LF"){
    max_columns = c("Run", "Raw.file", "Condition", 
                    "BioReplicate", "IsotopeLabelType", "Fraction")
  } else if (label_type == "TMT"){
    max_columns = c("Run", "Raw.file", "Fraction", "TechRepMixture", "Channel", 
                    "Condition", "Mixture", "BioReplicate")
  } else {
    msg = paste("Labeling type must be either LF or TMT")
    stop(msg)
  }
  
  if (sum(columns %in% max_columns) != length(columns)){
    msg = paste0(
    "Extra columns included in the annotation file that are not required. 
    This can cause problems in the converter dependencies. Please 
    only include the following columns in the annotation file: 
                 ", paste0(max_columns, collapse=", "))
    stop(msg)
  }
}



#' Check validity of parameters to MaxQ Converter function
#' @noRd
.checkMaxQconverterParams = function(mod.num = 'Total',
                                     keyword = "phos",
                                     which.proteinid.ptm = "Protein",
                                     which.proteinid.protein = "Leading.razor.protein",
                                     removeMpeptides = FALSE) {
  
  assertChoice(toupper(mod.num), c('SINGLE', 'TOTAL'), 
                          .var.name = 'NumberOfModifications')
  if (is.null(keyword)){
    msg = ('Keyword must not be null - stop')
    # getOption("MSstatsLog")("ERROR", msg)
    stop(msg)
  }
  assertChoice(which.proteinid.ptm, c("Proteins", "Leading.proteins",
                                                 "Protein"),
                          .var.name = "ProteinColumnNamePTM") 
  ## TODO: Use log to track these var choices
  # getOption("MSstatsLog")("INFO", paste("Summary method:", 
  #                                       summarization$method))
  assertChoice(which.proteinid.protein, c("Leading.proteins", 
                                                     "Leading.razor.protein", 
                                                     "Gene.names"),
                          .var.name = "ProteinColumnNameProtein")
  # getOption("MSstatsLog")("INFO", paste("cutOffCensored:", imputation$cutoff))
  assertLogical(removeMpeptides, .var.name = "removeMpeptides")
}

#' Function to check for global protein data
#' @noRd
.check.global.protein = function(evidence = NULL, 
                                  proteinGroups = NULL){
  
  if ((is.null(evidence) & !is.null(proteinGroups))|
      (!is.null(evidence) & is.null(proteinGroups))){
    msg = paste("To convert global protein data",
                "both the evidence and proteinGroups",
                "files must be provided - stop")
    #getOption("MSstatsLog")("ERROR", msg)
    stop(msg)
  
  } else if (is.null(evidence) & is.null(proteinGroups)){
    convert.global = FALSE
  } else {
    convert.global = TRUE
  }
  
  return(convert.global)
}

#' Check PTM and protein datasets
#' @noRd
.summarizeCheck = function(data, type) {
  # Check the PTM data
  if (is.null(data[["PTM"]])) {
    msg = paste('PTM peak list is missing. Input data must be of type list with',
                  'an element named \"PTM\" - stop')
    # getOption("MSstatsLog")("ERROR", msg)
    stop(msg)
  }
  ##TODO: Add a check to make sure Site was added into protein name?
  if (type == 'LF') {
    
    LF.cols = c("BioReplicate", "Condition", "FragmentIon", "Intensity", 
                 "IsotopeLabelType", "PeptideSequence", "PrecursorCharge",
                 "ProductCharge", "ProteinName", "Run")
    provided_cols = intersect(LF.cols, colnames(data[["PTM"]]))
    
    if (length(provided_cols) < 10){
      msg = paste("Missing columns in the PTM input:",
                   paste(setdiff(LF.cols, colnames(data[["PTM"]])), 
                         collapse = " "))
      #getOption("MSstatsLog")("ERROR", msg)
      stop(msg)
    }
  } else if (type == 'TMT'){
    ## Add peptide sequence if not available
    if (!'PeptideSequence' %in% colnames(data)){
      data[["PTM"]]$PeptideSequence = data[["PTM"]]$ProteinName
    }
    
    tmt.columns = c('ProteinName', 'PeptideSequence', 'Charge', 'PSM', 
                     'Mixture', 'TechRepMixture', 'Run', 'Channel', 'Condition',
                     'BioReplicate', 'Intensity')
    provided_cols = intersect(tmt.columns, colnames(data[["PTM"]]))
    
    if (length(provided_cols) < 11){
      msg = paste("Missing columns in the PTM input:",
                   paste(setdiff(tmt.columns, colnames(data[["PTM"]])),
                         collapse = " "))
      #getOption("MSstatsLog")("ERROR", msg)
      stop(msg)
    }
  }
  # Check the PROTEIN data
  if (is.null(data[["PROTEIN"]])) {
    adjust.protein = FALSE
  } else {
    adjust.protein = TRUE
    
    ## Check Protein data.table columns
    if (type == 'LF') {
      LF.cols = c("BioReplicate", "Condition", "FragmentIon", "Intensity", 
                   "IsotopeLabelType", "PeptideSequence", "PrecursorCharge",
                   "ProductCharge", "ProteinName", "Run")
      provided_cols = intersect(LF.cols, colnames(data[["PROTEIN"]]))
      
      if (length(provided_cols) < 10){
        msg = paste("Missing columns in the PROTEIN input:",
                     paste(setdiff(LF.cols, colnames(data[["PROTEIN"]])), 
                           collapse = " "))
        #getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
      }
    } else if (type == 'TMT'){
      tmt.columns = c('ProteinName', 'PeptideSequence', 'Charge', 'PSM', 
                       'Mixture', 'TechRepMixture', 'Run', 'Channel', 
                       'Condition', 'BioReplicate', 'Intensity')
      provided_cols = intersect(tmt.columns, colnames(data[["PROTEIN"]]))
      
      if (length(provided_cols) < 11){
        msg = paste("Missing columns in the PROTEIN input:",
                     paste(setdiff(tmt.columns, colnames(data[["PROTEIN"]])),
                           collapse = " "))
        #getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
      }
    }
  }
  ## Conditional for protein adjustment
  return(adjust.protein)
}

#' Check validity of parameters to dataProcess function.
#' @param log_base of logarithmic transformation
#' @param normalization_method string: "quantile", "equalizemedians", "FALSE",
#' "NONE" or "globalStandards"
#' @param normalization_method.PTM string: "quantile", "equalizemedians", 
#' "FALSE", "NONE" or "globalStandards"
#' @param address string
#' @param fill_rows logical, if TRUE, missing run observations for each feature
#' will be added with INTENSITY = NA
#' @param feature_selection list with elements: remove_uninformative
#' @param feature_selection.PTM list with elements: remove_uninformative
#' @param summarization list with elements: method.
#' @param imputation list with elements: cutoff, symbol.
#' @param n_clusters integer
#' @noRd
.checkDataProcessParams = function(log_base, normalization_method,
                                   normalization_method.PTM,
                                   standards_names, standards_names.PTM, 
                                   address,
                                   feature_selection, feature_selection.PTM,
                                   summarization, imputation) {
  
  checkmate::assertChoice(log_base, c(2, 10), .var.name = "logTrans")
  checkmate::assertChoice(summarization$method, c("linear", "TMP"),
                          .var.name = "summaryMethod")
  #getOption("MSstatsLog")("INFO", paste("Summary method:", 
  #                                      summarization$method))
  #getOption("MSstatsLog")("INFO", paste("cutOffCensored:", imputation$cutoff))
  checkmate::assertChoice(imputation$symbol, c("0", "NA"), 
                          null.ok = TRUE, .var.name = "censoredInt")
  #getOption("MSstatsLog")("INFO", paste("censoredInt:", imputation$symbol))
  if (summarization$method == "TMP" & imputation$MB & is.null(imputation$symbol)
      ) {
    msg = paste("The combination of required input",
                "MBimpute == TRUE and censoredInt = NULL",
                "has no censore missing values.",
                "Imputation will not be performed.- stop")
    getOption("MSstatsLog")("ERROR", msg)
    stop(msg)
  }
  checkmate::assertChoice(toupper(as.character(normalization_method)),
                          c("FALSE", "EQUALIZEMEDIANS", "QUANTILE", 
                            "GLOBALSTANDARDS"), .var.name = "normalization")
  checkmate::assertChoice(toupper(as.character(normalization_method.PTM)),
                          c("FALSE", "EQUALIZEMEDIANS", "QUANTILE", 
                            "GLOBALSTANDARDS"), .var.name = "PTM normalization")
  if (toupper(as.character(normalization_method)) == "GLOBALSTANDARDS" &
      is.null(standards_names)) {
    msg = paste("For normalization with global standards,",
                "the names of global standards are needed.",
                "Please add 'nameStandards' input.")
    #getOption("MSstatsLog")("ERROR", msg)
    stop(msg)
  }
  ## Remove this if global standards is not needed for PTM data
  if (toupper(as.character(normalization_method.PTM)) == "GLOBALSTANDARDS" &
      is.null(standards_names.PTM)) {
    msg = paste("For normalization with global standards,",
                "the names of global standards are needed.",
                "Please add 'nameStandards' input.")
    #getOption("MSstatsLog")("ERROR", msg)
    stop(msg)
  }
}

#' Check validity of parameters to TMT proteinSummarization function.
#' @param method string: "msstats", "MedianPolish", "Median", "LogSum"
#' @param global_norm logical: if True will normalize to protein-level
#' @param global_norm.PTM logical: if True will normalize to PTM-level
#' @param reference_norm logical: Indicates if reference channel was included
#' @param reference_norm,PTM logical: Indicates if reference channel was included
#' @param remove_norm_channel logical: remove Norm channel
#' @param remove_empty_channel logical: remove empty channels
#' @param MBimpute logical: only for method "msstats"
#' @param maxQuantileforCensored: Null or double
#' @noRd
.checkDataProcessParams.TMT = function(method, global_norm, global_norm.PTM,
                            reference_norm, reference_norm.PTM, 
                            remove_norm_channel, remove_empty_channel,
                            MBimpute, maxQuantileforCensored) {
  
  checkmate::assertChoice(method, c("msstats", "MedianPolish", "Median", 
                                    "LogSum"), .var.name = "Method")
  checkmate::assertLogical(global_norm, .var.name = "global_norm")
  checkmate::assertLogical(global_norm.PTM, .var.name = "global_norm.PTM")
  checkmate::assertLogical(reference_norm, .var.name = "reference_norm")
  checkmate::assertLogical(reference_norm.PTM, .var.name = "reference_norm.PTM")
  checkmate::assertLogical(remove_norm_channel, 
                           .var.name = "remove_norm_channel")
  checkmate::assertLogical(remove_empty_channel, 
                           .var.name = "remove_empty_channel")
  checkmate::assertLogical(MBimpute, .var.name = "MBimpute")
  
  if(!(is.double(maxQuantileforCensored)|is.null(maxQuantileforCensored))){
    msg = "maxQuantileforCensored must be either NULL or double - stop"
    stop(msg)
  }
}