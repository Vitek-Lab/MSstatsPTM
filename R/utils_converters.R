#' @noRd
.convert.ptm.data = function(pho.data,
                              annot,
                              mod.num,
                              keyword,
                              which.proteinid,
                              removeMpeptides){
  
  ## Remove Contaminant --------------------------------------------------------
  if (is.element("Contaminant", colnames(pho.data))) {
    pho.data = pho.data["Contaminant" == "", "Contaminant":=NA]
    pho.data = pho.data["Contaminant" != "+"]
  }
  
  if (is.element("Potential.contaminant", colnames(pho.data))) {
    pho.data = pho.data["Potential.contaminant" == "", 
                         "Potential.contaminant" :=NA]
    pho.data = pho.data["Potential.contaminant" != "+"]
  }
  
  ## Remove Reverse ------------------------------------------------------------
  if (is.element("Reverse", colnames(pho.data))) {
    pho.data = pho.data["Reverse" == "", "Reverse" := NA]
    pho.data = pho.data["Reverse" != "+"]
  }
  
  message('** + Contaminant, + Reverse, + Only.identified.by.site, PTMs are removed.')
  
  ## column for protein  -------------------------------------------------------
  ## default : Proteins
  which.pro = NULL
  
  if (which.proteinid == 'Proteins') {
    which.pro = 'Proteins'
  } else if (which.proteinid == 'Leading.proteins') {
    which.pro = 'Leading.proteins'
  } else if (which.proteinid == 'Protein') {
    which.pro = 'Protein'
  }
  
  if (which.proteinid == 'Protein' & !is.element('Protein', colnames(pho.data))) {
    which.pro = 'Leading.proteins'
    message('** Use Leading.proteins instead of Protein.')
  }
  
  if (which.proteinid == 'Leading.proteins' &
      !is.element('Leading.proteins', colnames(pho.data))) {
    which.pro = 'Proteins'
    message('** Use Proteins instead of Leading.proteins.')
  }
  
  ## at least 'Proteins' should be in the evidence.txt
  if (!is.element(which.pro, colnames(pho.data))) {
    stop('** Please select which columns should be used for protein ids,
         among three options (Proteins, Leading.proteins, Protein).')
  }
  
  if (mod.num == 'Single'){
    pep = "___1"
    pho.abun = .convert.peps(pho.data,
                              annot,
                              pep,
                              keyword,
                              which.pro,
                              removeMpeptides)
  }
  else if (mod.num == 'Total'){
    ## Extract and convert peps
    pho.pep.list = unique(unlist(lapply(
      colnames(pho.data), function(x) {str_extract(x, '___\\w+')})))
    pho.pep.list = pho.pep.list[!is.na(pho.pep.list)]
    
    temp.data.table.list = list()
    for (pep in pho.pep.list) {
      temp_pep = .convert.peps(pho.data,
                                annot,
                                pep,
                                keyword,
                                which.pro,
                                removeMpeptides)
      temp.data.table.list = append(temp.data.table.list, list(temp_pep))
    }
    names(temp.data.table.list) = pho.pep.list
    
    ## Combine into one datatable
    pho.abun = .pho.combine.pep(temp.data.table.list)
  }
  else{
    stop("mod.num must be one of Single or Total.")
  }
  
  #pho.abun = pho.abun[Abundance != 0 & !is.na(Abundance)]
  pho.abun$Run = paste(pho.abun$Mixture, pho.abun$TechRepMixture, sep = "_")
  
  colnames(pho.abun) = c('Replicate', 'PSM', 'Intensity', 'Channel', 
                          'Condition', 'Mixture', 'TechRepMixture', 
                          'BioReplicate', 'PeptideSequence', 'Charge', 
                          'ProteinName', 'Run')
  return(pho.abun)
  }


#' @noRd
.convert.peps = function(data,
                          annot,
                          pho_pep,
                          keyword,
                          which.pro,
                          removeMpeptides) {
  
  
  ## extract required columns --------------------------------------------------
  if (keyword == ""){
    required_data_columns = colnames(data)[
      grepl("Reporter.intensity.corrected.", colnames(data)) &
        lengths(regmatches(colnames(data), gregexpr("\\.", colnames(data)))
                ) >= 4
      & grepl(pho_pep, colnames(data))]
  } else {
    required_data_columns = colnames(data)[
    grepl("Reporter.intensity.corrected.", colnames(data)) &
      grepl(keyword, colnames(data)) & grepl(pho_pep, colnames(data))]
  }
  
  data = data[, c(which.pro, "Positions.within.proteins", "Amino.acid",
                   "Phospho..STY..Probabilities","Position.in.peptide",
                   "Charge", required_data_columns), with=FALSE]
  
  if (is.element('Proteins', colnames(data))){
    colnames(data)[colnames(data) == 'Proteins'] = 'ProteinName'
  }
  if (is.element('Leading.proteins', colnames(data))){
    colnames(data)[colnames(data) == 'Leading.proteins'] = 'ProteinName'
    ## Remove sites that are not in lead protein
    data$Positions.within.proteins = vapply(
      paste0(data$Positions.within.proteins, ';'), function(x){
        str_match(x, regex("^.*?(?=;)"))[1]})
  }
  if (is.element('Protein', colnames(data))){
    colnames(data)[colnames(data) == 'Protein'] = 'ProteinName'
    ## Remove sites that are not in lead protein
    data$Positions.within.proteins = sapply(
      paste0(data$Positions.within.proteins, ';'), function(x){
        str_match(x, regex("^.*?(?=;)"))[1]})
  }
  
  ## remove M site -------------------------------------------------------------
  ## only when removeMpeptides=TRUE
  if(removeMpeptides){
    remove_m_sequence = unique(data[
      grep("M", data$Phospho..STY..Probabilities),
      "Phospho..STY..Probabilities"])
    
    if(length(remove_m_sequence) > 0){
      data = data[-grep("M", data$Phospho..STY..Probabilities),]
    }
  }
  
  ## remove rows with all zero values ------------------------------------------
  ### some rows have all zero values across all MS runs. They should be removed.
  tmp = data[, c(required_data_columns), with=FALSE]
  count = apply(tmp, 1, function(x) sum(x==0))
  ## remove
  data = data[count != ncol(tmp), ]
  
  ## check whether there are multiple rows or not.------------------------------
  data = data[data$ProteinName!="", ]
  
  ## unique identifier for each row : phospho size + position + charge
  data$fea = paste(data$Phospho..STY..Probabilities, data$Charge,
                    data$ProteinName,
                    paste0(data$Amino.acid, data$Positions.within.proteins),
                    sep="_")
  
  count = xtabs(~fea, data)
  if(length(count[count > 1]) > 0){
    ## keep the features that have no issue.
    pho.single = data[which(data$fea %in% names(count[count == 1])), ]

    for(i in seq_len(length(names(count[count > 1])))){

      ## keep the feature that have issue.
      pho.mul = data[which(data$fea %in% names(count[count > 1])[i]), ]

      ## keep the feature with maximal summation of intensities
      pho.new = rbind(pho.single[, c(required_data_columns, 'fea'), with=FALSE],
                       pho.mul[
                         which.max(rowSums(pho.mul[, c(required_data_columns),
                                                   with=FALSE])),
                         c(required_data_columns, 'fea'), with=FALSE])
    }

  } else{
    pho.new = data[, c(required_data_columns, 'fea'), with=FALSE]
  }
  
  ## Assign the column names ---------------------------------------------------
  colnames(pho.new)[colnames(pho.new) == 'fea'] = 'PSM'
  
  ## make long format
  pho.l = melt(pho.new, id=c("PSM"))
  colnames(pho.l)[which(
    colnames(pho.l) %in% c('variable','value'))] = c("Replicate","Abundance")
  pho.l$Replicate = gsub(pho_pep, "", pho.l$Replicate)
  
  ## merge with annotation
  pho.l = merge(pho.l, annot, by = 'Replicate', all.x = TRUE)
  
  ## replace zero with NA for MQ output
  pho.l[pho.l$Abundance == 0, 'Abundance'] = NA
  pho.l$Abundance = as.numeric(as.character(pho.l$Abundance))
  #pho.l$Abundance = log2(pho.l$Abundance)
  pho.l$PeptideSequence = str_match(pho.l$PSM, regex("^(?:[^_]*){1}"))
  pho.l$Charge = str_match(pho.l$PSM, regex("(?<=_)[0-9](?=_)")[1])
  pho.l$Protein = str_replace(pho.l$PSM, fixed(paste0(pho.l$PeptideSequence, 
                                                       '_', pho.l$Charge, '_')),
                               "")
  pep.pho.abun = pho.l
  
  return(pep.pho.abun)
}

#' @noRd
.pho.combine.pep = function(pep.list){
  
  ## Columns to join on
  join.cols = colnames(pep.list[[1]])[
    colnames(pep.list[[1]]) != 'Abundance']
  
  ## Combine pep data.tables
  pho.abun = Reduce(function(x,y) {
    merge(x, y, by = join.cols, all = TRUE)}, pep.list)
  
  ## sum pep abundance
  pho.abun$Total = rowSums(
    pho.abun[, c(colnames(pho.abun)[
      colnames(pho.abun) %like% 'Abundance']),
      with=FALSE], na.rm=TRUE)
  
  ## Clean combined data.table
  pho.abun = pho.abun[, c(colnames(pho.abun)[
    !colnames(pho.abun) %like% 'Abundance']), with = FALSE]
  colnames(pho.abun)[colnames(pho.abun) == 'Total'] = 'Abundance'
  setcolorder(pho.abun, c("PSM", "Replicate", "Abundance", "Channel", 
                          "Condition", "Mixture", "TechRepMixture", 
                          "BioReplicate", "Protein", "Charge"))
  
  combined.pho.abun = pho.abun
  return(combined.pho.abun)
}

#' Add modified data into protein for Progenesis input
#' @noRd
.progensis.add.sites = function(ptm_input, fasta_path, col_order){
  
  X.8 = X.9 = X.10 = NULL
  
  ## Load fasta file
  format_fasta = tidyFasta(fasta_path)
  
  mod_df = ptm_input[X.9 != "", ]
  
  ## Extract peptide and protein info
  peptide = mod_df[3:nrow(mod_df), X.8]
  protein = mod_df[3:nrow(mod_df), X.10]
  mods = mod_df[3:nrow(mod_df), X.9]
  
  ## Extract modifications
  loc_df = data.table("uniprot_iso" = protein, "peptide" = peptide, 
                       "mods" = mods)
  loc_df = merge(loc_df, format_fasta, by = c("uniprot_iso"), all.x = TRUE)
  loc_df = loc_df[!is.na(sequence),]
  
  mod_locations = regmatches(loc_df[, mods], gregexpr("[[:digit:]]+", 
                                                       loc_df[, mods]))
  
  start = sapply(seq_len(nrow(loc_df)),
                  function(i) gregexpr(loc_df$peptide[i], 
                                       loc_df$sequence[i])[[1]])
  mod_index = sapply(seq_len(length(mod_locations)), function(i){
    mod_locations[i] = list(
      as.integer(unlist(mod_locations[i])) + start[i])})
  
  peptide_mod = mapply(get_sites, mod_locations, mod_index,
                        loc_df$peptide)
  
  loc_df$site = peptide_mod
  loc_df$site = paste(loc_df$uniprot_iso, loc_df$site, sep = "_")
  
  ptm_input = merge(ptm_input, 
                     unique(loc_df[, c("uniprot_iso", "peptide", 
                                       "mods", "site")]), 
                     by.x = c("X.10" , "X.8", "X.9"), 
                     by.y = c("uniprot_iso", "peptide", "mods"), all.x = TRUE)
  
  ## Put data back into Progenesis format
  ptm_input = ptm_input[order(ptm_input$id), ]
  ptm_input$X.10 = ptm_input$site
  ptm_input[2, "X.10"] = "Accession"
  ptm_input[1, "X.10"] = ""
  ptm_input = ptm_input[, ..col_order]
  
  ## Remove proteins missing from fasta file
  ptm_input = ptm_input[!is.na(X.10), ]
  
  return(ptm_input)
}

#' Helper function for progenesis site identification
#' @noRd
get_sites = function(str_idx, pep_idx, pep) {
  
  if (!is.na(pep_idx[1])){
    mod_site = ""
    for (i in seq(str_idx)){
      site = substr(pep, as.numeric(str_idx[[i]]), as.numeric(str_idx[[i]]))
      temp = paste0(site, as.character(pep_idx[i]), "_")
      mod_site = paste0(mod_site, temp)
    }
    mod_site = substr(mod_site,1,nchar(mod_site)-1)
  } else {mod_site = NA}
  return(mod_site)
}

#' Helper function for spectronaut site identification
#' combine with prot one at some point
#' @noRd
spectro_get_sites = function(str_idx, start, pep) {
  
  mod_site = ""
  for (i in seq_along(str_idx)){
    site = substr(pep, as.numeric(str_idx[[i]])-i, as.numeric(str_idx[[i]])-i)
    pep_idx = start[[1]] + as.numeric(str_idx[[i]])-i-1
    temp = paste0(site, as.character(pep_idx), "_")
    mod_site = paste0(mod_site, temp)
  }
  mod_site = substr(mod_site,1,nchar(mod_site)-1)
  
  return(mod_site)
}


#' Locate modification site number and amino acid
#' 
#' @param data `data.table` of enriched experimental run. Must include 
#' `ProteinName`, `PeptideSequence`, `PeptideModifiedSequence`, and (optionally)
#' `Start` columns.
#' @param protein_name_col Name of column indicating protein. Default is 
#' `ProteinName`.
#' @param unmod_pep_col Name of column indicating unmodified peptide sequence. Default
#' is `PeptideSequence`.
#' @param mod_pep_col Name of column indicating modified peptide sequence. Default
#' is `PeptideModifiedSequence`.
#' @param clean_mod Remove special characters and numbers around modification 
#' name. Default is `FALSE`
#' @param fasta_file File path to FASTA file that matches with proteins in 
#' `data`. Can be either string or `data.table` processed with `tidyFasta()` 
#' function. Default to NULL if peptide number included in `data`.
#' @param fasta_protein_name Name of fasta file column that matches with 
#' `protein_name_col`. Default is `header`.
#' @param mod_id String that indicates what amino acid was modified in 
#' `PeptideSequence`.
#' @param mod_id_is_numeric Boolean indicating if mod id is numeric. If TRUE, 
#' `mod_id` will be ignored. Default is FALSE.
#' @param terminus_included Boolean indicating if the `PeptideSequence` includes
#' the terminus amino acid.
#' @param terminus_id String that indicates what the terminus amino acid is. 
#' Default is '.'.
#' 
#' @importFrom data.table as.data.table
#' 
#' @return `data.table` with site location added into `Protein` column.
#' @export 
#' 
#' @examples
#' ##TODO
#' 
MSstatsPTMSiteLocator = function(data, 
                                 protein_name_col= "ProteinName",
                                 unmod_pep_col = "PeptideSequence",
                                 mod_pep_col = "PeptideModifiedSequence",
                                 clean_mod=FALSE,
                                 fasta_file=NULL, 
                                 fasta_protein_name="header",
                                 mod_id="\\*", 
                                 mod_id_is_numeric=FALSE, 
                                 terminus_included=FALSE, 
                                 terminus_id="\\."){
  
  ## Check if peptide number included in data
  if ("Start" %in% colnames(data)){
    id_data = data[,c(protein_name_col, unmod_pep_col, 
                      mod_pep_col, "Start")]
  } else {
    if (is.null(fasta_file)){
      stop("FASTA file not provided and `Start` column missing from data. \
           MSstatsPTMSiteLocator requires one of these to be provided to \
           identify modification site number.")
    }
    
    id_data = data[,c(protein_name_col, unmod_pep_col, mod_pep_col), with=FALSE]
    id_data = .joinFasta(id_data, fasta_file, fasta_protein_name,
                         protein_name_col, unmod_pep_col, mod_pep_col)
  }
  
  if (terminus_included){
    id_data = .fixTerminus(id_data, terminus_id, unmod_pep_col)
  }
  
  if (clean_mod){
    id_data[,mod_pep_col] = lapply(id_data[,mod_pep_col], 
                                   function(x){gsub("[0-9]+|[[:punct:]]", "", x)})
  }
  
  id_data = .locateSites(id_data, mod_id, mod_id_is_numeric,
                         protein_name_col, unmod_pep_col, mod_pep_col)
  id_data = id_data[!duplicated(id_data),]
  
  data = merge(data, id_data, by=c(protein_name_col, unmod_pep_col, 
                                   mod_pep_col))
  
  data[, "ProteinNameUnmod"] = data[, protein_name_col, with=FALSE]
  data[, protein_name_col] = data$ProteinName_mod
  data=as.data.table(data)
  
  if ("Start.y" %in% colnames(data)){
    data$Start = data$Start.y
  }
  remove = intersect(c("PeptideModifiedSequence_adj", "ProteinName_mod", 
                       "start_fix", "Start.x", "Start.y"),
                     colnames(data))
  data[, (remove):=NULL]
  return(data)
}
  
#' Add FASTA data into dataframe
#' @param data data.table
#' @param fasta_file string or data.table
#' @return data.table
#' @keywords internal
.joinFasta = function(data, fasta_file, fasta_protein_name,
                      protein_name_col, unmod_pep_col, mod_pep_col){
  
  ## Load FASTA if not already loaded
  if (is.character(fasta_file)){
    fasta_file = tidyFasta(fasta_file)
  }
  
  ## Join data and FASTA
  data = merge(data, fasta_file, by.x=protein_name_col, 
                  by.y=fasta_protein_name, all.x=TRUE)

  missing_prot = unique(data[is.na(data$sequence), 
                             protein_name_col, with=FALSE])
  if (nrow(missing_prot) > 0){
    print(paste0("FASTA file missing ", as.character(length(missing_prot)),
                 " Proteins. These will be removed. This may be due to non-unique identifications."))
  }
  
  data = data[!is.na(data$sequence),]
              
  ## Locate peptide sequence start
  data$Start = mapply(function(x,y){unlist(gregexpr(pattern=x, y))[[1]]}, 
                      data[, unmod_pep_col, with=FALSE][[1]], 
                      data[, "sequence"][[1]])
  data$Start = unlist(data$Start)
  
  data = data[, c(protein_name_col, unmod_pep_col, mod_pep_col, "Start"), 
              with=FALSE]
  
  return(data)
}
#' Add site location and aa
#' @param data data.table
#' @param fasta_file string or data.table
#' @return data.table
#' @keywords internal
.fixTerminus = function(data, terminus_id, unmod_pep_col){
  
  ## Terminus makes start location off
  data$start_fix = lapply(data[unmod_pep_col], function(x){
    as.integer(unlist(gregexpr(terminus_id, 
                               substr(x, start=1, stop=2))) == 2)})
  data$start_fix = unlist(data$start_fix)
  data$Start = data$Start - unlist(data$start_fix)
  
  return(data)
}

#' Add site location and aa
#' @param data data.table
#' @param mod_id string
#' @return data.table
#' @keywords internal
.locateSites = function(data, mod_id, mod_id_is_numeric,
                        protein_name_col, unmod_pep_col, mod_pep_col){
  
  data[, "PeptideModifiedSequence_adj"] = data[,mod_pep_col, with=FALSE]
  
  if (mod_id_is_numeric){
    mod_id="\\*"
    ## Replace mod locations with *
    data$PeptideModifiedSequence_adj = lapply(data$PeptideModifiedSequence_adj, 
                                           function(x){gsub("[0-9.()]", "", x)})
    data$PeptideModifiedSequence_adj = lapply(data$PeptideModifiedSequence_adj, 
                                           function(x){gsub("\\[", "", x)})
    data$PeptideModifiedSequence_adj = lapply(data$PeptideModifiedSequence_adj, 
                                              function(x){gsub("\\]", "", x)})
    data$PeptideModifiedSequence_adj = lapply(data$PeptideModifiedSequence_adj, 
                                              function(x){gsub("\\+", 
                                                               mod_id, x)})
    ## Fix instances of double mod aa (caused by Acetylation)
    data$PeptideModifiedSequence_adj = lapply(data$PeptideModifiedSequence_adj, 
                                           function(x){
                                             gsub("\\*\\*", mod_id, x)})
    data$PeptideModifiedSequence_adj = unlist(data$PeptideModifiedSequence_adj)
  }
  
  unmod_data = data[!grepl(mod_id, data[,"PeptideModifiedSequence_adj"][[1]]), ]
  mod_data = data[grepl(mod_id, data[,"PeptideModifiedSequence_adj"][[1]]),]
  
  ## Locate number and aa  
  string_num = lapply(mod_data$PeptideModifiedSequence_adj, 
                      function(x){unlist(gregexpr(mod_id, x))})
  
  string_num = lapply(string_num, function(x){x - 1:length(x)})
  site_num = mapply(function(x,y){x + y - 1}, string_num, mod_data$Start)
  
  site_aa = lapply(mod_data$PeptideModifiedSequence_adj, function(x){
    unlist(lapply(unlist(gregexpr(mod_id, x)) - 1, 
                  function(y){substr(x, y, y)}))})
  
  ## Combine number and aa
  full_site = lapply(mapply(function(x,y){paste(x, y, sep="")}, 
                            site_aa, site_num),
                     function(z){paste(z, collapse='_')})
  mod_data$ProteinName_mod = paste(mod_data[, protein_name_col, with=FALSE][[1]], 
                                   full_site, sep="_")
  
  ## Ensure columns are not lists
  mod_data[,unmod_pep_col] = mod_data[, unmod_pep_col, with=FALSE][[1]]
  mod_data$PeptideModifiedSequence_adj = unlist(
    mod_data$PeptideModifiedSequence_adj)
  mod_data$ProteinName_mod = unlist(mod_data$ProteinName_mod)
  
  unmod_data$ProteinName_mod = unmod_data[, protein_name_col, with=FALSE][[1]]
  data = rbindlist(list(mod_data, unmod_data))
  
  return(data)
}

#' Pivot Peak Studio data into long format
#' @noRd
#' @keywords internal
.pivotPS = function(input){
  
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
  
  pivot_input$FragmentIon <- NA
  pivot_input$ProductCharge <- NA
  pivot_input$PrecursorCharge <- NA
  pivot_input$IsotopeLabelType <- 'L'
  
  return(pivot_input)
}

#' Remove modifications from pep sequence that are not required
#' @noRd
#' @keywords internal
.MSstatsPTMRemoveMods = function(sequence, remove_special_char=TRUE){
  if (remove_special_char){
    sequence = lapply(sequence, function(x){gsub("\\[(.*?)\\]|[0-9]+|[[:punct:]]", "", x)})
  } else {
    sequence = lapply(sequence, function(x){gsub("\\[(.*?)\\]|[0-9]+", "", x)})
  }
  
  return(unlist(sequence))
}

#' Pull out modifications from PD PTM data for input into protein name
#' @noRd
#' @keywords internal
.extract_pd_mods = function(modifications, mod_id){
  split_mods = str_split(modifications, ";")
  
  quick_filter = function(mod_list){
    mod_list = str_trim(mod_list)
    mask = ifelse(grepl(mod_id, mod_list), TRUE, FALSE)
    return(mod_list[mask])
  }
  
  target_mods = lapply(split_mods, quick_filter)
  join_mods = lapply(target_mods, function(x){paste(x, collapse="_")})
  
  return(unlist(join_mods))
}


#' Remake MaxQ MSstatsConvert::.cleanRawMaxQuant function to not require proteinGroups file
#' 
#' @noRd
#' @keywords internal
#' @import MSstatsConvert
.cleanRawMaxQuantPTM = function(msstats_object, protein_id_col, data.type,
                                remove_by_site = FALSE,
                                channel_columns = "Reporterintensitycorrected"){
  
  ProteinIDs = id = PSM = PeptideSequence = NULL
  ProteingroupIDs = PrecursorCharge = Intensity = NULL
  channels = character(0)
  
  mq_input = MSstatsConvert::getInputFile(msstats_object, "evidence")
  
  filter_cols = c("Contaminant", "Potentialcontaminant", "Reverse")
  msg = paste("** + Contaminant, + Reverse, + Potential.contaminant",
              "proteins are removed.")
  if (remove_by_site) {
    filter_cols = c(filter_cols, "Onlyidentifiedbysite")
    msg = paste("** + Contaminant, + Reverse, + Potential.contaminant,",
                "+ Only.identified.by.site proteins are removed.")
  }
  
  mq_input = MSstatsConvert:::.filterManyColumns(mq_input, filter_cols, "+")
  
  
  if (data.type == "TMT") {
    channels = MSstatsConvert:::.getChannelColumns(colnames(mq_input), 
                                                   channel_columns)
  }
  
  data.table::setnames(
    mq_input, 
    c(protein_id_col, "Modifiedsequence", "Charge", "Rawfile"), 
    c("ProteinName", "PeptideSequence", "PrecursorCharge", "Run"),
    skip_absent = TRUE)
  mq_input[["PeptideSequence"]] = gsub("_", "", mq_input[["PeptideSequence"]])
  mq_cols = c("ProteinName", "PeptideSequence", "Modifications", 
              "PrecursorCharge", "Run", "Intensity", 
              "Fraction", "TechReplicate", "Run", "BioReplicate",
              "PSM", "Score")
  mq_cols = intersect(c(mq_cols, channels),
                      colnames(mq_input))
  mq_input = mq_input[, mq_cols, with = FALSE]
  mq_input = unique(mq_input)
  
  if (data.type == "TMT") {
    mq_input[, PSM := paste(PeptideSequence, PrecursorCharge, 
                            1:nrow(mq_input), sep = "_")]
    mq_input = melt(mq_input, measure.vars = channels,
                    id.vars = c("ProteinName", "PeptideSequence", 
                                "PrecursorCharge", "PSM", "Run", "Score"),
                    variable.name = "Channel", value.name = "Intensity")
    mq_input$Channel = gsub(channel_columns, "channel", mq_input$Channel)
    mq_input$Channel = MSstatsConvert:::.standardizeColnames(mq_input$Channel)
    suppressWarnings({
      mq_input$Intensity = ifelse(mq_input$Intensity == 0, NA,
                                  mq_input$Intensity)
    })
    mq_input = MSstatsConvert:::.filterFewMeasurements(mq_input, 0, FALSE, 
                                      c("PeptideSequence", 
                                        "PrecursorCharge", "Run"))
  }
  
  mq_input = mq_input[!is.na(Intensity), ]
  
  return(mq_input)
}

#' Remake MSstatstoMaxQuant converter to work without proteinGroups file
#' 
#' @noRd
#' @keywords internal
#' @import MSstatsConvert
MaxQtoMSstatsFormatHelper = function(
    evidence, annotation, data.type = "LF",   
    proteinID = "Proteins",
    useUniquePeptide = TRUE,
    summaryforMultipleRows = max,
    removeFewMeasurements = TRUE,
    removeMpeptides = FALSE,
    removeOxidationMpeptides = FALSE,
    removeProtein_with1Peptide = FALSE,
    use_log_file = TRUE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = NULL) {
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path)
  
  input = MSstatsConvert::MSstatsImport(list(evidence = evidence), 
                                        type = "MSstats",
                                        tool = "MaxQuant")
  input = .cleanRawMaxQuantPTM(input, proteinID, data.type)
  annotation = MSstatsConvert::MSstatsMakeAnnotation(input, 
                                                     annotation, 
                                                     "Run" = "Rawfile")
  
  m_filter = list(col_name = "PeptideSequence", 
                  pattern = "M", 
                  filter = removeMpeptides, 
                  drop_column = FALSE)
  
  oxidation_filter = list(col_name = "Modifications", 
                          pattern = "Oxidation", 
                          filter = removeOxidationMpeptides, 
                          drop_column = TRUE)
  
  feature_columns = c("PeptideSequence", "PrecursorCharge")
  input = MSstatsConvert::MSstatsPreprocess(
    input, 
    annotation,
    feature_columns,
    remove_shared_peptides = useUniquePeptide, 
    remove_single_feature_proteins = removeProtein_with1Peptide,
    pattern_filtering = list(oxidation = oxidation_filter,
                             m = m_filter),
    feature_cleaning = list(
      remove_features_with_few_measurements = removeFewMeasurements,
      summarize_multiple_psms = summaryforMultipleRows),
    columns_to_fill = list("FragmentIon" = NA,
                           "ProductCharge" = NA,
                           "IsotopeLabelType" = "L"))
  input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns)
  
  msg_final = paste("** Finished preprocessing. The dataset is ready",
                    "to be processed by the dataProcess function.")
  getOption("MSstatsLog")("INFO", msg_final)
  getOption("MSstatsMsg")("INFO", msg_final)
  getOption("MSstatsLog")("INFO", "\n")
  
  return(input)
}

#' Remake MSstatstoMaxQuantTMT converter to work without proteinGroups file
#' 
#' @noRd
#' @keywords internal
#' @import MSstatsConvert
MaxQtoMSstatsTMTFormatHelper = function(
    evidence, annotation, data.type="TMT", which.proteinid = 'Proteins',
    rmProt_Only.identified.by.site = FALSE, useUniquePeptide = TRUE,
    rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE, 
    summaryforMultipleRows = sum, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL
) {
  MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                      log_file_path, 
                                      base = "MSstatsTMT_converter_log_")
  
  input = MSstatsConvert::MSstatsImport(list(evidence = evidence), 
                                        "MSstatsTMT", "MaxQuant")
  input = .cleanRawMaxQuantPTM(input, which.proteinid, data.type,
                               remove_by_site = rmProt_Only.identified.by.site,
                               channel_columns = "Reporterintensitycorrected")
  annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
  
  feature_columns = c("PeptideSequence", "PrecursorCharge")
  input = MSstatsConvert::MSstatsPreprocess(
    input, 
    annotation, 
    feature_columns,
    remove_shared_peptides = useUniquePeptide,
    remove_single_feature_proteins = rmProtein_with1Feature,
    feature_cleaning = list(remove_features_with_few_measurements = rmPSM_withfewMea_withinRun,
                            summarize_multiple_psms = summaryforMultipleRows))
  input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                fix_missing = "zero_to_na")
  data.table::setnames(input, "PrecursorCharge", "Charge", skip_absent = TRUE)
  
  msg_final = paste("** Finished preprocessing. The dataset is ready",
                    "to be processed by the proteinSummarization function.")
  getOption("MSstatsLog")("INFO", msg_final)
  getOption("MSstatsMsg")("INFO", msg_final)
  getOption("MSstatsLog")("INFO", "\n")
  input
}