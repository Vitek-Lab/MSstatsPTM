#' Convert Peaks Studio output into MSstatsPTM format
#' 
#' Currently only supports label-free quantification.
#' 
#' @param input name of Peaks Studio PTM output
#' @param annotation name of annotation file which includes Raw.file, Condition,
#' BioReplicate, Run. For example annoation see example below.
#' @param input_protein name of Peaks Studio unmodified protein output 
#' (optional)
#' @param annotation name of annotation file which includes Raw.file, Condition,
#' BioReplicate, Run for unmodified protein output.
#' @param use_unmod_peptides
#' @importFrom data.table as.data.table melt
#' @importFrom dataframe as.data.frame
#' @return list of data.table
#' @export 
#' 
#' @examples
#' #TODO: Add examples
PStoMSstatsPTMFormat = function(
    input, annotation, input_protein = NULL, annotation_protein = NULL,
    use_unmod_peptides = FALSE, target_modification = NULL, 
    remove_oxidation_peptides = FALSE, remove_multi_mod_types = FALSE,
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL
    ){
  
  input = as.data.table(input)
  input_prot = NULL
  message("Pivoting input data..")
  input = pivotPS(input)
  input = add_missing_cols(input)
  message("Merging with annotation..")
  input = merge(input, annotation, all.x = TRUE, by = "Raw.File")

  if (!is.null(input_protein)){
    message("Converting unmodified protein data..")
    ## TODO: Add unmod protein converter
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
  input = add_mod(input)
  
  input[, c("Raw.File", "Mod", "End", "Start"):=NULL]
  
  feature_cols = c("ProteinName", "PeptideSequence", "FragmentIon", 
                   "ProductCharge", "PrecursorCharge", "IsotopeLabelType",
                   "BioReplicate", "Condition", "Run" )
  input = input[, list(Intensity = max(Intensity, na.rm = TRUE)), 
                by = feature_cols]
  
  if (!is.null(input_prot)){
    input_prot[, c("Raw.File", "Mod", "End", "Start"):=NULL]
    input_prot = input_prot[, list(Intensity = max(Intensity, 
                                                          na.rm = TRUE)), 
                            by = feature_cols]
    input_prot = as.data.frame(input_prot)
    
  }
  
  return(list("PTM" = as.data.frame(input),
              "PROTEIN" = input_prot))
}

#' Pivot PS data into long format
#' @noRd
#' @keywords internal
pivotPS = function(input){
  
  columns = colnames(input)
  len_cols = length(columns)
  end_runs = which(grepl("os", columns))[[1]] - 1
  keep_cols = c(3, 4, len_cols, len_cols-1, len_cols-2, 11:end_runs)
  
  
  filter_input = input[, ..keep_cols]
  
  
  pivot_input = melt(filter_input, 
                     measure.vars = 6:length(colnames(filter_input)), 
                     variable.name = "Raw.File", value.name = "Intensity")
  colnames(pivot_input) = c("ProteinName", "PeptideSequence", 
                            "Mod", "End", "Start", "Raw.File", "Intensity")
  
  return(pivot_input)
}

#' Add missing colums to PS
#' @noRd
#' @keywords internal
add_missing_cols = function(input){
  input$FragmentIon <- NA
  input$ProductCharge <- NA
  input$PrecursorCharge <- NA
  input$IsotopeLabelType <- 'L'
  return(input)
}

#' Add modification to protein name
#' @noRd
#' @keywords internal
add_mod = function(input){
  ## Replace mod locations with *
  input$PeptideSequence = lapply(input$PeptideSequence, 
                                 function(x){gsub("[0-9.()]", "", x)})
  input$PeptideSequence = lapply(input$PeptideSequence, 
                                 function(x){gsub("[+]", "*", x)})
  ## Fix instances of double mod aa (caused by Acetylation)
  input$PeptideSequence = lapply(input$PeptideSequence, 
                                 function(x){gsub("\\*\\*", "*", x)})
  
  string_num = lapply(input$PeptideSequence, function(x){unlist(gregexpr("\\*", x))})
  string_num = lapply(string_num, function(x){x - 1:length(x)})
  site_num = mapply(function(x,y){x + y - 2}, string_num, input$Start)
  
  site_aa = lapply(input$PeptideSequence, function(x){
    unlist(lapply(unlist(gregexpr("\\*", x)) - 1, function(y){substr(x, y, y)}))})
  
  
  full_site = lapply(mapply(function(x,y){paste(x, y, sep="")}, 
                            site_aa, site_num),
                     function(z){paste(z, collapse='_')})
  input$ProteinName = paste(input$ProteinName, full_site, sep="_")
  
  input$PeptideSequence = unlist(input$PeptideSequence)
  input$ProteinName = unlist(input$ProteinName)
  
  return(input)
}