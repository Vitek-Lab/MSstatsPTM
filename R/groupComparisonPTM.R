#' Perform differential analysis on MS-based proteomics experiments targeting PTMs
#'
#' Takes summarized PTM and protein data from `dataSummarizationPTM` or 
#' `dataSummarizationPTM_TMT` functions and performs differential analysis. 
#' Leverages unmodified protein data to perform adjustment and deconvolute the 
#' effect of the PTM and unmodified protein. If protein data is unavailable, 
#' PTM data can still be passed into the function, however adjustment can not be
#' performed. All model results are returned for completeness.
#' 
#' @export
#' @importFrom data.table data.table as.data.table `:=`
#' @importFrom MSstats groupComparison MSstatsContrastMatrix
#' @importFrom MSstatsTMT groupComparisonTMT
#' @importFrom MSstatsConvert MSstatsLogsSettings MSstatsSaveSessionInfo
#' @importFrom stats p.adjust xtabs
#' @importFrom stringr str_match
#' 
#' @param data list of summarized datasets. Output of MSstatsPTM summarization 
#' function \code{\link[MSstatsPTM]{dataSummarizationPTM}}  or 
#' \code{\link[MSstatsPTM]{dataSummarizationPTM_TMT}} depending on acquisition
#' type.
#' @param data.type Type of data. Must be one of `LF` or `TMT`. Will be deprecated 
#' in favor of ptm_label_type and protein_label_type.
#' @param contrast.matrix comparison between conditions of interests. Default 
#' models full pairwise comparison between all conditions
#' @param moderated For TMT experiments only. TRUE will moderate t statistic; 
#' FALSE (default) uses ordinary t statistic. Default is FALSE.
#' @param adj.method For TMT experiments only. Adjusted method for multiple 
#' comparison. "BH" is default. "BH" is used for all other experiment types
#' @param log_base For non-TMT experiments only. The base of the logarithm used 
#' in summarization.
#' @param save_fitted_models logical, if TRUE, fitted models will be added to 
#' the output.
#' @param use_log_file logical. If TRUE, information about data processing
#' will be saved to a file.
#' @param append logical. If TRUE, information about data processing will be 
#' added to an existing log file.
#' @param verbose logical. If TRUE, information about data processing will be 
#' printed to the console.
#' @param log_file_path character. Path to a file to which information about 
#' data processing will be saved. 
#' If not provided, such a file will be created automatically.
#' If `append = TRUE`, has to be a valid path to a file.
#' @param base start of the file name.
#' @param ptm_label_type Indicator of labeling type for PTM dataset. Must be one
#' of `LF` or `TMT`
#' @param protein_label_type Indicator of labeling type for PROTEIN dataset. 
#' Must be one of `LF` or `TMT`
#' @return list of modeling results. Includes PTM, PROTEIN, and ADJUSTED
#'         data.tables with their corresponding model results.
#'         
#' @examples 
#' 
#' model.lf.msstatsptm = groupComparisonPTM(summary.data, 
#'                                      ptm_label_type="LF",
#'                                      protein_label_type="LF",
#'                                      verbose = FALSE)
groupComparisonPTM = function(data, 
                              data.type = NULL,
                              contrast.matrix = "pairwise",
                              moderated = FALSE, 
                              adj.method = "BH",
                              log_base = 2,
                              save_fitted_models=TRUE,
                              use_log_file = TRUE, 
                              append = FALSE,
                              verbose = TRUE, 
                              log_file_path = NULL,
                              base = "MSstatsPTM_log_",
                              ptm_label_type = "LF",
                              protein_label_type = "LF") {
  
  ## Start log  
  # if (is.null(log_file_path) & use_log_file == TRUE){
  #   time_now = Sys.time()
  #   path = paste0(base, gsub("[ :\\-]", "_", time_now), 
  #                 ".log")
  #   file.create(path)
  # } else {path = log_file_path}
  # 
  # if (data.type == 'TMT'){
  #   pkg = "MSstatsTMT"
  #   option_log = "MSstatsTMTLog"
  # } else {
  #   pkg = "MSstats"
  #   option_log = "MSstatsLog"
  # }
  # 
  # MSstatsLogsSettings(use_log_file, append,
  #                     verbose, log_file_path = path, 
  #                     pkg_name = pkg)
  
  # getOption(option_log)("INFO", "Starting parameter and data checks..")
  if (!is.null(data.type) && (data.type == "TMT" || data.type == "LF")) {
      ptm_label_type = data.type
      protein_label_type = data.type
  }
  
  Label = Site = NULL
  
  data.ptm = data[["PTM"]]
  data.protein = data[["PROTEIN"]]
  
  ## Check for missing variables in PTM
  ## Determine if PTM should be adjusted for protein level.
  if (!is.null(data.protein)){
    adj.protein = TRUE
  } else{
    adj.protein = FALSE
  }
  
  ## Create pairwise matrix for label free
  if (contrast.matrix[1] == "pairwise"){
    # getOption(option_log)("INFO", "Building pairwise matrix.")
    if ("GROUP" %in% colnames(data.ptm$ProteinLevelData)) {
      labels <- unique(data.ptm$ProteinLevelData$GROUP)
    } else if ("Condition" %in% colnames(data.ptm$ProteinLevelData)) {
      labels <- unique(data.ptm$ProteinLevelData$Condition)
    } else {
      stop("Error: Neither 'GROUP' nor 'Condition' column exists in data.ptm$ProteinLevelData.")
    }
    contrast.matrix = MSstatsContrastMatrix('pairwise', labels)
  }
  
  ## PTM Modeling
  message("Starting PTM modeling...")
  if (ptm_label_type == "TMT"){
    # getOption(option_log)("INFO", "Starting TMT PTM Model")
    ptm_model_full = groupComparisonTMT(data.ptm, 
                                        contrast.matrix = contrast.matrix,
                                        moderated = moderated, 
                                        adj.method = adj.method,
                                        save_fitted_models=save_fitted_models)#,
                                        # save_fitted_models = TRUE,
                                        # use_log_file = use_log_file, 
                                        # append = append, verbose = verbose, 
                                        # log_file_path = path)
    ptm_model = ptm_model_full$ComparisonResult
    ptm_model_site_sep = ptm_model_full$ComparisonResult
    ptm_model_details = ptm_model_full$FittedModel
  } else if (ptm_label_type == "LF") {
    # getOption(option_log)("INFO", "Starting non-TMT PTM Model")
    ptm_model_full = groupComparison(contrast.matrix,
                                      data.ptm, save_fitted_models, log_base)#, 
                                      # use_log_file, append, verbose, 
                                      # log_file_path = path)
    ptm_model = ptm_model_full$ComparisonResult
    ptm_model_site_sep = ptm_model_full$ComparisonResult
    ptm_model_details = ptm_model_full$FittedModel
  }
  
  models = list('PTM.Model'=ptm_model, 'Model.Details'=ptm_model_details)
  
  if (adj.protein) {
    
    ## Protein Modeling
    message("Starting Protein modeling...")
    if (protein_label_type == "TMT"){
      # getOption(option_log)("INFO", "Starting TMT Protein Model")
      protein_model_full = groupComparisonTMT(data.protein, 
                                          contrast.matrix = contrast.matrix,
                                          moderated = moderated, 
                                          adj.method = adj.method,
                                          save_fitted_models = save_fitted_models)
                                          # use_log_file = use_log_file,
                                          # append = append,
                                          # verbose = verbose, 
                                          # log_file_path = path)
      protein_model = protein_model_full$ComparisonResult
      protein_model_details = protein_model_full$FittedModel
    } else if (protein_label_type == "LF") {
      # getOption(option_log)("INFO", "Starting non-TMT Protein Model")
      protein_model_full = groupComparison(contrast.matrix, 
                                           data.protein, save_fitted_models, 
                                           log_base, use_log_file)#, 
                                            # append, verbose, 
                                            # log_file_path = path)
      protein_model = protein_model_full$ComparisonResult
      protein_model_details = protein_model_full$FittedModel
    }
    
    ptm_model = as.data.table(ptm_model)
    protein_model = as.data.table(protein_model)
    
    message("Starting adjustment...")
    # getOption(option_log)("INFO", "Starting Protein Adjustment")
    ptm_model_site_sep = copy(ptm_model)
    
    ## extract global protein name
    ptm_model_site_sep = .extractProtein(ptm_model_site_sep, protein_model)
    # getOption(option_log)("INFO", "Rcpp function extracted protein info")
    
    ## adjustProteinLevel function can only compare one label at a time
    comparisons = unique(ptm_model_site_sep[, Label])
    
    adjusted_model_list = list()
    for (i in seq_len(length(comparisons))) {
      # getOption(option_log)("INFO", paste0("Adjusting for Comparison - ", 
                                             # as.character(i)))
      temp_adjusted_model = .applyPtmAdjustment(comparisons[[i]],
                                                   ptm_model_site_sep,
                                                   protein_model)
      adjusted_model_list[[i]] = temp_adjusted_model
    }
    
    adjusted_models = rbindlist(adjusted_model_list)
    
    adjusted_models$GlobalProtein = adjusted_models$Protein
    adjusted_models$Protein = adjusted_models$Site
    adjusted_models[, Site := NULL]
    
    ## Add unadjustable ptms into final dataframe
    adjusted_models$temp_check = paste(adjusted_models$Protein, 
                                       adjusted_models$Label, sep = "_")
    ptm_model$temp_check = paste(ptm_model$Protein, 
                                 ptm_model$Label, sep = "_")
    missing_ptms = setdiff(ptm_model$temp_check, adjusted_models$temp_check)
    
    missing_rows = ptm_model[ptm_model$temp_check %in% missing_ptms]
    missing_rows$GlobalProtein = "missing"
    missing_rows$Adjusted = FALSE
    missing_rows[, c('issue', 'MissingPercentage', 
                     'ImputationPercentage', 'temp_check'):=NULL]
    
    adjusted_models$Adjusted = TRUE
    
    if (ptm_label_type == 'TMT'){
      missing_rows$Tvalue = missing_rows$log2FC / missing_rows$SE
    }
    
    ptm_model$temp_check = NULL
    adjusted_models$temp_check = NULL
    
    adjusted_models = rbindlist(list(adjusted_models, missing_rows), 
                                use.names=TRUE)
    adjusted_models = adjusted_models[!is.na(adjusted_models$Protein)]
    
    # getOption(option_log)("INFO", "Adjustment complete, returning models.")
    models = list('PTM.Model'=ptm_model, 
                  'PROTEIN.Model'=protein_model,
                  'ADJUSTED.Model'=adjusted_models, 
                  'Model.Details'=list('PTM'=ptm_model_details,
                                       'PROTEIN'=protein_model_details))
  }
  
  if (!is.null(data.type) && (data.type == "TMT" || data.type == "LF")) {
      warning("DEPRECATION NOTICE: The `data.type` argument is being deprecated. 
              Please use `ptm_label_type` and `protein_label_type` instead ahead
              of Release 3.22")
  }

  return(models)

}
