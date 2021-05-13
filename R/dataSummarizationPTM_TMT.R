#' Process MS PTM and global protein data produced via tandem mass tag labeling
#'
#' Utilizes functionality from MSstatsTMT to clean, summarize, and
#' normalize PTM and protein level data. Imputes missing values, protein and PTM
#' level summarization from peptide level quantification. Applies global median
#' normalization on peptide level data and normalizes between runs.
#'
#' @export
#' @importFrom MSstatsTMT proteinSummarization
#' @importFrom MSstatsConvert MSstatsLogsSettings MSstatsSaveSessionInfo
#' @importFrom data.table as.data.table is.data.table
#' @param data Name of the output of MSstatsPTM converter function or 
#' peptide-level quantified data from other tools. It should be a list 
#' containing one or two data tables, named PTM and PROTEIN for modified and 
#' unmodified datasets. The list must at least contain the PTM dataset. The data
#' should have columns ProteinName, PeptideSequence, Charge, PSM, Mixture, 
#' TechRepMixture, Run, Channel, Condition, BioReplicate, Intensity
#' @param method Four different summarization methods to protein-level can be 
#' performed : "msstats"(default), "MedianPolish", "Median", "LogSum".
#' @param global_norm Global median normalization on for unmodified peptide 
#' level data (equalizing the medians across all the channels and MS runs). 
#' Default is TRUE. It will be performed before protein-level summarization.
#' @param global_norm.PTM Same as above for modified peptide level data. Default
#' is TRUE.
#' @param reference_norm Reference channel based normalization between MS runs 
#' on unmodified protein level data. TRUE(default) needs at least one reference 
#' channel in each MS run, annotated by 'Norm' in Condtion column. It will be 
#' performed after protein-level summarization. FALSE will not perform this 
#' normalization step. If data only has one run, then reference_norm=FALSE.
#' @param reference_norm.PTM Same as above for modified peptide level data. 
#' Default is TRUE.
#' @param remove_norm_channel TRUE(default) removes 'Norm' channels from 
#' protein level data.
#' @param remove_empty_channel TRUE(default) removes 'Empty' channels from 
#' protein level data.
#' @param MBimpute only for method="msstats". TRUE (default) imputes missing 
#' values by Accelated failure model. FALSE uses minimum value to impute the 
#' missing value for each peptide precursor ion.
#' @param MBimpute.PTM Same as above for modified peptide level data. Default is
#' TRUE
#' @param maxQuantileforCensored We assume missing values are censored. 
#' maxQuantileforCensored is Maximum quantile for deciding censored missing 
#' value, for instance, 0.999. Default is Null.
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
#' @return list of two data.tables
#' @examples
#' head(raw.input.tmt$PTM)
#' head(raw.input.tmt$PROTEIN)
#' 
#' quant.tmt.msstatsptm <- dataSummarizationPTM_TMT(raw.input.tmt,
#'                                                  method = "msstats", 
#'                                                  verbose = FALSE)
#' head(quant.tmt.msstatsptm$PTM)
dataSummarizationPTM_TMT <- function(
  data,
  method = "msstats",
  global_norm = TRUE,
  global_norm.PTM = TRUE,
  reference_norm = TRUE,
  reference_norm.PTM = TRUE,
  remove_norm_channel = TRUE,
  remove_empty_channel = TRUE,
  MBimpute = TRUE,
  MBimpute.PTM = TRUE,
  maxQuantileforCensored = NULL,
  use_log_file = TRUE, 
  append = FALSE,
  verbose = TRUE, 
  log_file_path = NULL
  ) {
  
  ## Start log  
  if (is.null(log_file_path) & use_log_file == TRUE){
    time_now <- Sys.time()
    path <- paste0("MSstatsPTM_log_", gsub("[ :\\-]", "_", time_now), 
                  ".log")
    file.create(path)
  } else {path <- log_file_path}
  
  MSstatsLogsSettings(use_log_file, append,
                      verbose, log_file_path = path,
                      base = "MSstatsPTM_log_", pkg_name = "MSstatsTMT")
  
  getOption("MSstatsTMTLog")("INFO", "Starting parameter and data checks..")
  
  .checkDataProcessParams.TMT(method, global_norm, global_norm.PTM,
                          reference_norm, reference_norm.PTM, 
                          remove_norm_channel, remove_empty_channel,
                          MBimpute, maxQuantileforCensored)

  adj.protein <- FALSE
  
  PTM.dataset <- data[["PTM"]]
  
  protein.dataset <- data[["PROTEIN"]]
  
  # Check PTM and PROTEIN data for correct format
  adj.protein <- .summarizeCheck(data, 'TMT')
  PTM.dataset <- as.data.table(PTM.dataset)
  
  getOption("MSstatsTMTLog")("INFO", "Parameter and data checks complete.")
  
  ## Determine if protein level should also be summarized
  if (adj.protein) {
    getOption("MSstatsTMTLog")("INFO", "Protein dataset was included.")
    protein.dataset <- as.data.table(protein.dataset)
  }
  
  ## Run summarization function from MSstatsTMT
  getOption("MSstatsTMTLog")("INFO", "Starting PTM summarization..")
  ptm.summarized <- proteinSummarization(PTM.dataset,
                                         method, global_norm.PTM, 
                                         reference_norm.PTM,
                                         remove_norm_channel,
                                         remove_empty_channel, 
                                         MBimpute.PTM,
                                         maxQuantileforCensored, use_log_file, 
                                         append, verbose, log_file_path = path,
                                         msstats_log_path = path)
  
  if (adj.protein) {
    getOption("MSstatsTMTLog")("INFO", "Starting Protein summarization..")
    protein.summarized <- proteinSummarization(protein.dataset,
                                               method, global_norm, 
                                               reference_norm,
                                               remove_norm_channel,
                                               remove_empty_channel, 
                                               MBimpute,
                                               maxQuantileforCensored, 
                                               use_log_file, append, verbose, 
                                               log_file_path = path)
  }
  
  ## Compile and return summarized results
  getOption("MSstatsTMTLog")("INFO", "Summarization complete. Returning output")
  msstatsptm.summarized <- list("PTM" = ptm.summarized)
  if (adj.protein) {
    msstatsptm.summarized <- c(msstatsptm.summarized,
                               "PROTEIN" = list(protein.summarized))
  }
  
  return(msstatsptm.summarized)
}