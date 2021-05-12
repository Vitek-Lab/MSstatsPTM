#' Extract global protein from combined protein and site column
#' @noRd
.extractProtein <- function(ptm_model, protein_model){
  
  ## All proteins
  available_proteins <- unique(as.character(protein_model$Protein))
  available_proteins <- available_proteins[order(nchar(available_proteins),
                                                 available_proteins,
                                                 decreasing = TRUE)]
  available_ptms <- unique(as.character(ptm_model$Protein))
  
  ## Set site
  ptm_model$Site <- ptm_model$Protein
  ## Call Rcpp function
  ptm_proteins <- extract_protein_name(available_ptms,
                                       available_proteins)
  global_protein_lookup <- data.table(Site = available_ptms, 
                                      Protein = ptm_proteins)
  
  ## Add extracted protein name into model
  ptm_model[, "Protein" := NULL]
  ptm_model <- merge(ptm_model, global_protein_lookup,
                     all.x = TRUE, by = 'Site')

  return(ptm_model)
}

#' Apply global protein adjustment to each condition
#' @noRd
.applyPtmAdjustment <- function(label, ptm_model, protein_model){
  
  Label <- NULL
  
  temp_ptm_model <- ptm_model[Label == label]
  temp_protein_model <- protein_model[Label == label]
  
  ## Function from MSstatsPTM Compare
  temp_adjusted_model <- .adjustProteinLevel(temp_ptm_model, temp_protein_model)
  temp_adjusted_model$adj.pvalue <- p.adjust(temp_adjusted_model$pvalue,
                                             method = 'BH')
  temp_adjusted_model
}

#' Calculate parameter changes
#' @noRd
.adjustProteinLevel <- function(diff_site, diff_protein) {
  diff_ref <- diff_protein[, c("Protein", "Label", "log2FC", "SE", "DF")]
  names(diff_ref)[names(diff_ref) == "log2FC"] <- "log2FC_ref"
  names(diff_ref)[names(diff_ref) == "SE"] <- "SE_ref"
  names(diff_ref)[names(diff_ref) == "DF"] <- "DF_ref"
  joined <- merge(diff_site, diff_ref, by = c("Protein", "Label"))
  
  missing_ctrl <- joined[joined$log2FC == Inf, ]  # PTM missing in control
  missing_case <- joined[joined$log2FC == -Inf, ] # PTM missing in case
  res_mctrl <- data.table(Protein = missing_ctrl$Protein,
                          Site = missing_ctrl$Site, Label = missing_ctrl$Label,
                          log2FC = Inf, SE = NA, Tvalue = NA, DF = NA, 
                          pvalue = NA)
  res_mcase <- data.table(Protein = missing_case$Protein,
                          Site = missing_case$Site, Label = missing_case$Label,
                          log2FC = -Inf, SE = NA, Tvalue = NA, DF = NA, 
                          pvalue = NA)
  
  idx_full <- abs(joined$log2FC) != Inf & abs(joined$log2FC_ref) != Inf
  full <- joined[idx_full, ]
  log2fc <- full$log2FC - full$log2FC_ref
  s2 <- full$SE ^ 2
  s2_ref <- full$SE_ref ^ 2
  stderr <- sqrt(s2 + s2_ref)
  numer <- (s2 + s2_ref) ^ 2
  denom <- (s2 ^ 2 / full$DF + s2_ref ^ 2 / full$DF_ref)
  df <- numer / denom
  tval <- log2fc / stderr
  pval <- 2 * stats::pt(abs(tval), df, lower.tail = FALSE)

  res_full <- data.table(Protein = full$Protein,
                         Site = full$Site,
                         Label = full$Label,
                         log2FC = log2fc,
                         SE = stderr,
                         Tvalue = tval,
                         DF = df,
                         pvalue = pval)
  
  rbindlist(list(res_full, res_mctrl, res_mcase))
}

#' make contrast matrix for pairwise comparisons
#' @noRd
.makeContrast <- function(groups) {
  
  ncomp <- length(groups) * (length(groups) - 1) / 2 # Number of comparison
  contrast_matrix <- matrix(rep(0, length(groups) * ncomp), 
                            ncol = length(groups))
  colnames(contrast_matrix) <- groups
  
  count <- 0
  contrast_matrix_rownames <- NULL
  for(j in seq_len(length(groups)-1)){
    for(k in (j+1):length(groups)){
      
      count <- count + 1
      # save row name
      contrast_matrix_rownames <- c(contrast_matrix_rownames, 
                                    paste(groups[j], groups[k], sep = "-"))
      # set constrast value
      contrast_matrix[count, groups[j]] <- 1
      contrast_matrix[count, groups[k]] <- -1
    }
  }
  rownames(contrast_matrix) <- contrast_matrix_rownames
  
  return(contrast_matrix)
}
