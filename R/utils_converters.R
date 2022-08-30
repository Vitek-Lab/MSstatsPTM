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
