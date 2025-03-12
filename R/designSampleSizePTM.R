#' Planning future experimental designs of PTM experiments in sample size calculation
#'
#' @description Calculate sample size for future experiments of a PTM 
#' experiment based on intensity-based linear model. Calculation is only 
#' available for group comparison experimental designs (not including time 
#' series).
#' Two options of the calculation: 
#' (1) number of biological replicates per condition, 
#' (2) power.
#' 
#' @param data output of the groupComparisonPTM function.
#' @param desiredFC the range of a desired fold change which includes the lower 
#' and upper values of the desired fold change.
#' @param FDR a pre-specified false discovery ratio (FDR) to control the overall 
#' false positive rate. Default is 0.05
#' @param numSample minimal number of biological replicates per condition. 
#' TRUE represents you require to calculate the sample size for this category, 
#' else you should input the exact number of biological replicates.
#' @param power a pre-specified statistical power which defined as the probability 
#' of detecting a true fold change. TRUE represent you require to calculate the power 
#' for this category, else you should input the average of power you expect. Default is 0.9
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
#' 
#' @details The function fits the model and uses variance components to calculate 
#' sample size. The underlying model fitting with intensity-based linear model with 
#' technical MS run replication. Estimated sample size is rounded to 0 decimal.
#' The function can only obtain either one of the categories of the sample size 
#' calculation (numSample, numPep, numTran, power) at the same time.
#' 
#' @return data.frame - sample size calculation results including varibles:
#' desiredFC, numSample, FDR,  and power.
#' 
#' @importFrom MSstats designSampleSize
#' @importFrom stats median
#' 
#' @export
#' 
#' @examples
#' model.lf.msstatsptm = groupComparisonPTM(summary.data, 
#'                                      ptm_label_type="LF",
#'                                      protein_label_type="LF",
#'                                      verbose = FALSE)
#'                                      
#' #(1) Minimal number of biological replicates per condition
#' designSampleSizePTM(data=model.lf.msstatsptm, numSample=TRUE,
#'                  desiredFC=c(2.0,2.75), FDR=0.05, power=0.8)
#' #(2) Power calculation
#' designSampleSizePTM(data=model.lf.msstatsptm, numSample=5,
#'                  desiredFC=c(2.0,2.75), FDR=0.05, power=TRUE)                                  
#' 
#' 
designSampleSizePTM = function(
  data, desiredFC, FDR = 0.05, numSample = TRUE, power = 0.8, 
  use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
  base = "MSstatsPTM_log_"
){
  
  ## Start log  
  if (is.null(log_file_path) & use_log_file == TRUE){
    time_now = Sys.time()
    path = paste0(base, gsub("[ :\\-]", "_", time_now), 
                  ".log")
    file.create(path)
  } else {path = log_file_path}
  
  # if (data.type == 'TMT'){
  #   pkg = "MSstatsTMT"
  #   option_log = "MSstatsTMTLog"
  # } else {
  #   pkg = "MSstats"
  #   option_log = "MSstatsLog"
  # }
  pkg = "MSstats"
  option_log = "MSstatsLog"
  
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path, 
                                      "MSstatsPTM_sampleSize_log_")
  getOption(option_log)("INFO", "** MSstatsPTM - designSampleSizePTM function")
  getOption(option_log)("INFO", paste0("Desired fold change = ", 
                                         paste(desiredFC, collapse=" - ")))
  getOption(option_log)("INFO", paste0("FDR = ", FDR))
  getOption(option_log)("INFO", paste0("Power = ", power))
  
  fitted_models = data$Model.Details
  
  if (!is.null(fitted_models$PROTEIN)){
    
    ptm_var_component = MSstats:::.getVarComponent(fitted_models$PTM)
    protein_var_component = MSstats:::.getVarComponent(fitted_models$PROTEIN)
    
    ptm_median_sigma_error = median(ptm_var_component[["Error"]], na.rm = TRUE)
    protein_median_sigma_error = median(protein_var_component[["Error"]], 
                                        na.rm = TRUE)
    
    ptm_median_sigma_subject = MSstats:::.getMedianSigmaSubject(
      ptm_var_component)
    protein_median_sigma_subject = MSstats:::.getMedianSigmaSubject(
      protein_var_component)
    
    ## power calculation
    if (isTRUE(power)) {
      delta = log2(seq(desiredFC[1], desiredFC[2], 0.025))
      desiredFC = 2 ^ delta
      power_output = .calculatePowerPTM(desiredFC, FDR, delta, 
                                        ptm_median_sigma_error, 
                                        protein_median_sigma_error, 
                                        ptm_median_sigma_subject,
                                        protein_median_sigma_subject, numSample)
      getOption("MSstatsLog")("INFO", "Power is calculated. - okay")
      sample_size = data.frame(desiredFC, numSample, FDR, 
                               power = power_output)
    }	
    
    if (is.numeric(power)) {
      delta = log2(seq(desiredFC[1], desiredFC[2], 0.025))
      desiredFC = 2 ^ delta
      ## Large portion of proteins are not changing
      m0_m1 = 99 ## it means m0/m1=99, m0/(m0+m1)=0.99
      alpha = power * FDR / (1 + (1 - FDR) * m0_m1)
      if (isTRUE(numSample)) {
        numSample = .getNumSamplePTM(desiredFC, power, alpha, delta,
                                     ptm_median_sigma_error, 
                                     protein_median_sigma_error, 
                                     ptm_median_sigma_subject,
                                     protein_median_sigma_subject)
        
        getOption(option_log)("INFO", 
                              "The number of sample is calculated. - okay")
        sample_size = data.frame(desiredFC, numSample, FDR, power)
      }
    } 
    
    
  } else {
    sample_size = designSampleSize(fitted_models$PTM, desiredFC, FDR, numSample,
                                   power, use_log_file, append, verbose, 
                                   log_file_path)
  }
  
  return(sample_size)
  
}

#' Power calculation for PTM experiment
#' @inheritParams designSampleSizePTM
#' @importFrom stats qnorm
#' @return `float` of power
#' @keywords internal
.calculatePowerPTM = function(desiredFC, FDR, delta, ptm_median_sigma_error, 
                           protein_median_sigma_error, ptm_median_sigma_subject,
                           protein_median_sigma_subject, numSample) {
  m0_m1 = 99
  t = delta / sqrt(2 * (ptm_median_sigma_error/numSample + 
                          protein_median_sigma_error/numSample + 
                          ptm_median_sigma_subject/numSample + 
                          protein_median_sigma_subject/numSample))
  powerTemp = seq(0, 1, 0.01)
  power = numeric(length(t))
  for (i in seq_along(t)) {
    diff = qnorm(powerTemp) + qnorm(1 - powerTemp * 
                                      FDR / (1 + (1 - FDR) * m0_m1) / 2) - t[i]
    min(abs(diff), na.rm = TRUE)
    power[i] = powerTemp[order(abs(diff))][1]
  }
  return(power)
}

#' Get sample size for PTM experiment
#' @inheritParams designSampleSizePTM
#' @inheritParams .calculatePowerPTM
#' @importFrom stats qnorm
#' @return `int` of samples
#' @keywords internal
.getNumSamplePTM = function(desiredFC, power, alpha, delta, 
                            ptm_median_sigma_error, protein_median_sigma_error, 
                            ptm_median_sigma_subject, 
                            protein_median_sigma_subject){
  z_alpha = qnorm(1 - alpha / 2)
  z_beta = qnorm(power)
  aa = (delta / (z_alpha + z_beta)) ^ 2
  numSample = round(2 * (ptm_median_sigma_error + protein_median_sigma_error + 
                           ptm_median_sigma_subject + 
                           protein_median_sigma_subject) / aa, 0)
  return(numSample)
}
