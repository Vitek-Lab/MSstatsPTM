#' Data summarization function for label-free MS experiments targeting PTMs.
#'
#' Utilizes functionality from MSstats to clean, summarize, and
#' normalize PTM and protein level data. Imputes missing values, performs 
#' normalization, and summarizes data. PTM data is summarized up to the 
#' modification and protein data is summarized up to the protein level. Takes 
#' as input the output of the included converters (see included `raw.input` 
#' data object for required input format).
#'
#' @export
#' @importFrom MSstats dataProcess
#' @importFrom MSstatsConvert MSstatsLogsSettings MSstatsSaveSessionInfo
#' @importFrom data.table as.data.table is.data.table
#' 
#' @param data name of the list with PTM and (optionally) unmodified protein 
#' data.tables, which can be the output of the MSstatsPTM converter functions
#' @param logTrans logarithm transformation with base 2(default) or 10
#' @param normalization normalization for the protein level dataset, to remove 
#' systematic bias between MS runs. There are three different normalizations 
#' supported. 'equalizeMedians'(default) represents constant normalization 
#' (equalizing the medians) based on reference signals is performed. 'quantile' 
#' represents quantile normalization based on reference signals is performed. 
#' 'globalStandards' represents normalization with global standards proteins. 
#' FALSE represents no normalization is performed
#' @param normalization.PTM normalization for PTM level dataset. Default is 
#' "equalizeMedians" Can be adjusted to any of the options described above.
#' @param nameStandards vector of global standard peptide names for protein 
#' dataset. only for normalization with global standard peptides.
#' @param nameStandards.PTM Same as above for PTM dataset.
#' @param featureSubset "all" (default) uses all features that the data set has. 
#' "top3" uses top 3 features which have highest average of log-intensity across runs. 
#' "topN" uses top N features which has highest average of log-intensity across runs. 
#' It needs the input for n_top_feature option. 
#' "highQuality" flags uninformative feature and outliers.
#' @param featureSubset.PTM For PTM dataset only. Options same as above.
#' @param remove_uninformative_feature_outlier For protein dataset only. It only
#'  works after users used featureSubset="highQuality" in dataProcess. TRUE 
#'  allows to remove 1) the features are flagged in the column, 
#'  feature_quality="Uninformative" which are features with bad quality, 2) 
#'  outliers that are flagged in the column, is_outlier=TRUE, for run-level 
#'  summarization. FALSE (default) uses all features and intensities for 
#'  run-level summarization.
#' @param remove_uninformative_feature_outlier.PTM For PTM dataset only. Options
#' same as above.
#' @param min_feature_count optional. Only required if featureSubset = "highQuality".
#' Defines a minimum number of informative features a protein needs to be considered
#' in the feature selection algorithm.
#' @param min_feature_count.PTM For PTM dataset only. Options the same as above.
#' Default is 1 due to low average feature count for PTMs. 
#' @param n_top_feature For protein dataset only. The number of top features for
#'  featureSubset='topN'. Default is 3, which means to use top 3 features.
#' @param n_top_feature.PTM For PTM dataset only. Options same as above.
#' @param summaryMethod "TMP"(default) means Tukey's median polish, which is 
#' robust estimation method. "linear" uses linear mixed model. 
#' @param equalFeatureVar only for summaryMethod="linear". default is TRUE. 
#' Logical variable for whether the model should account for heterogeneous 
#' variation among intensities from different features. Default is TRUE, which 
#' assume equal variance among intensities from features. FALSE means that we 
#' cannot assume equal variance among intensities from features, then we will 
#' account for heterogeneous variation from different features.
#' @param censoredInt Missing values are censored or at random. 'NA' (default) 
#' assumes that all 'NA's in 'Intensity' column are censored. '0' uses zero 
#' intensities as censored intensity. In this case, NA intensities are missing 
#' at random. The output from Skyline should use '0'. Null assumes that all NA 
#' intensites are randomly missing.
#' @param MBimpute For protein dataset only. only for summaryMethod="TMP" and 
#' censoredInt='NA' or '0'. TRUE (default) imputes 'NA' or '0' (depending on 
#' censoredInt option) by Accelated failure model. FALSE uses the values 
#' assigned by cutoffCensored.
#' @param MBimpute.PTM For PTM dataset only. Options same as above.
#' @param remove50missing only for summaryMethod="TMP". TRUE removes the runs 
#' which have more than 50% missing values. FALSE is default.
#' @param maxQuantileforCensored Maximum quantile for deciding censored missing 
#' values. default is 0.999
#' @param fix_missing Default is Null. Optional, same as the 'fix_missing' 
#' parameter in MSstatsConvert::MSstatsBalancedDesign function
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
#' @return list of summarized PTM and Protein results. These results contain 
#' the reformatted input to the summarization function, as well as run-level
#' summarization results.
#' @examples
#' head(raw.input$PTM)
#' head(raw.input$PROTEIN)
#' 
#' quant.lf.msstatsptm = dataSummarizationPTM(raw.input, verbose = FALSE)
#' head(quant.lf.msstatsptm$PTM$ProteinLevelData)
dataSummarizationPTM = function(
  data,
  logTrans = 2, 
  normalization = "equalizeMedians", 
  normalization.PTM = "equalizeMedians", 
  nameStandards = NULL,
  nameStandards.PTM = NULL,
  featureSubset = "all",
  featureSubset.PTM = "all",
  remove_uninformative_feature_outlier = FALSE, 
  remove_uninformative_feature_outlier.PTM = FALSE,
  min_feature_count = 2,
  min_feature_count.PTM = 1,
  n_top_feature = 3, 
  n_top_feature.PTM = 3,
  summaryMethod = "TMP", 
  equalFeatureVar = TRUE, 
  censoredInt = "NA",
  MBimpute = TRUE, 
  MBimpute.PTM = TRUE, 
  remove50missing = FALSE,
  fix_missing = NULL,
  maxQuantileforCensored = 0.999,
  use_log_file = TRUE, 
  append = TRUE,
  verbose = TRUE, 
  log_file_path = NULL,
  base = "MSstatsPTM_log_") {

  ## Start log  
  if (is.null(log_file_path) & use_log_file == TRUE){
    time_now = Sys.time()
    path = paste0(base, gsub("[ :\\-]", "_", time_now), 
                  ".log")
    file.create(path)
  } else {path = log_file_path}
  
  MSstatsLogsSettings(use_log_file, append,
                      verbose, log_file_path = path)
  
  getOption("MSstatsLog")("INFO", "Starting parameter and data checks..")
  
  .checkDataProcessParams(logTrans, normalization, normalization.PTM, 
                          nameStandards, nameStandards.PTM, 
                          address,
                          list(method = featureSubset, n_top = n_top_feature,
                               remove_uninformative = 
                                 remove_uninformative_feature_outlier),
                          list(method = featureSubset.PTM, 
                               n_top = n_top_feature.PTM,
                               remove_uninformative = 
                                 remove_uninformative_feature_outlier.PTM),
                          list(method = summaryMethod, 
                               equal_var = equalFeatureVar),
                          list(symbol = censoredInt,
                               MB = MBimpute))
  
  adj.protein = FALSE

  PTM.dataset = data[["PTM"]]
  protein.dataset = data[["PROTEIN"]]
  
  # Check PTM and PROTEIN data for correct format
  adj.protein = .summarizeCheck(data, 'LF')
  
  getOption("MSstatsLog")("INFO", "Parameter and data checks complete.")
  
  message("Starting PTM summarization...")
  ## Run summarization function from MSstats
  getOption("MSstatsLog")("INFO", "Starting PTM summarization..")
  ptm.summarized = dataProcess(as.data.frame(PTM.dataset),
                                logTrans, normalization.PTM, nameStandards.PTM,
                                featureSubset.PTM, 
                                remove_uninformative_feature_outlier.PTM,
                                min_feature_count.PTM, n_top_feature.PTM,
                                summaryMethod, equalFeatureVar, censoredInt,
                                MBimpute.PTM, remove50missing, fix_missing,
                                maxQuantileforCensored, use_log_file, append,
                                verbose, log_file_path = path)
  if (adj.protein) {
    message("Starting Protein summarization...")
    getOption("MSstatsLog")("INFO", "Starting Protein summarization..")
    protein.summarized = dataProcess(as.data.frame(protein.dataset),
                                  logTrans, normalization, nameStandards,
                                  featureSubset, 
                                  remove_uninformative_feature_outlier,
                                  min_feature_count, n_top_feature,
                                  summaryMethod, equalFeatureVar, censoredInt,
                                  MBimpute, remove50missing, fix_missing,
                                  maxQuantileforCensored, use_log_file, append,
                                  verbose, log_file_path = path)
  }

  ## Compile and return summarized results
  getOption("MSstatsLog")("INFO", "Summarization complete. Returning output")
  msstatsptm.summarized = list("PTM" = ptm.summarized)
  if (adj.protein) {
    msstatsptm.summarized = c(msstatsptm.summarized,
                               "PROTEIN" = list(protein.summarized))
  }

  return(msstatsptm.summarized)
}

