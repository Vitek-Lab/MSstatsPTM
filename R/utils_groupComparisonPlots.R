#' Check input data and parameters
#' @noRd
.check.plotting.data = function(data, type, sig, FCcutoff, logBase.pvalue, 
                                 ylimUp, ylimDown, xlimUp, x.axis.size, 
                                 y.axis.size, dot.size, text.size, text.angle, 
                                 legend.size, ProteinName,colorkey, numProtein, 
                                 width, height, which.Comparison, which.Protein,
                                 address){
  
  ## Column check
  assertChoice(toupper(type), c("HEATMAP", "VOLCANOPLOT"),
               .var.name = "Type") 
  assertNumeric(sig, .var.name = "Sig")
  if (!is.logical(FCcutoff)){
    assertNumeric(FCcutoff, .var.name = "FCcutoff")
  }
  
  assertChoice(logBase.pvalue, c(2, 10),
               .var.name = "Logbase Pvalue") 
  
  if (!is.logical(ylimUp)){
    assertNumeric(ylimUp, .var.name = "ylimUp")
  }
  if (!is.logical(ylimDown)){
    assertNumeric(ylimDown, .var.name = "ylimDown")
  }
  if (!is.logical(xlimUp)){
    assertNumeric(xlimUp, .var.name = "xlimUp")
  }
  
  ## Check plotting size vars
  assertNumeric(x.axis.size, .var.name = ("x.axis.size"))
  assertNumeric(y.axis.size, .var.name = ("y.axis.size"))
  assertNumeric(dot.size, .var.name = ("dot.size"))
  assertNumeric(text.size, .var.name = ("text.size"))
  assertNumeric(text.angle, .var.name = ("text.angle"))
  assertNumeric(legend.size, .var.name = ("legend.size"))
  
  assertLogical(colorkey, .var.name = ("colorkey"))
  
  assertNumeric(width, .var.name = ("width"))
  assertNumeric(height, .var.name = ("height"))
  
  ## Check data columns
  min.cols = c("Protein", "Label", "log2FC", "SE", "DF", "pvalue","adj.pvalue")
  
  if (!is.null(data[["PTM.Model"]])){
    missing = setdiff(min.cols, colnames(data[["PTM.Model"]]))
    if (length(missing) > 0){
      msg = paste("Missing columns in the PTM input:",
                  paste(missing, collapse = " "))
      stop(msg)
    }
  }
  if (!is.null(data[["PROTEIN.Model"]])){
    missing = setdiff(min.cols, colnames(data[["PROTEIN.Model"]]))
    if (length(missing) > 0){
      msg = paste("Missing columns in the Protein input:",
                   paste(missing, collapse = " "))
      stop(msg)
    }
  }
  if (!is.null(data[["ADJUSTED.Model"]])){
    missing = setdiff(min.cols, colnames(data[["ADJUSTED.Model"]]))
    if (length(missing) > 0){
      msg = paste("Missing columns in the Adjusted PTM input:",
                   paste(missing, collapse = " "))
      stop(msg)
    }
  }
                
}

#' Format models into plotting format
#' @noRd
.format.model.plots = function(data, which.Comparison, which.Protein){
  
  ptm.model = data[['PTM.Model']]
  ptm.model = as.data.table(ptm.model)
  protein.model = data[['PROTEIN.Model']]
  protein.model = as.data.table(protein.model)
  adjusted.model = data[['ADJUSTED.Model']]
  adjusted.model = as.data.table(adjusted.model)
  
  ptm.model$Protein = factor(ptm.model$Protein)
  ptm.model$Label = factor(ptm.model$Label)
  
  if (nrow(protein.model) > 0){
    protein.model$Protein = factor(protein.model$Protein)
    protein.model$Label = factor(protein.model$Label)
    adjusted.model$Protein = factor(adjusted.model$Protein)
    adjusted.model$Label = factor(adjusted.model$Label)
  }
  
  if (which.Comparison[[1]] != "all") {
    ## check which.comparison is name of comparison
    if (is.character(which.Comparison)) {
      
      temp.name = which.Comparison
      ## message if name of comparison is wrong.
      if (length(setdiff(temp.name, unique(ptm.model$Label))) > 0) {
        
        ## TODO: Logging        
        # processout = rbind(processout, 
        #                     paste0("Please check labels of comparions. ", 
        #                            "Result does not have this comparison. - ", 
        #                            toString(temp.name)))
        # write.table(processout, file=finalfile, row.names=FALSE)
        
        stop(paste0("Please check labels of comparisons. ", 
                    "Result does not have this comparison. - ", 
                    toString(temp.name)))
      }
    }
    ## check which.comparison is order number of comparison
    if (is.numeric(which.Comparison)) {
      
      temp.name = levels(ptm.model$Label)[which.Comparison]
      
      ## message if name of comparison is wrong.
      if (length(levels(ptm.model$Label))<max(which.Comparison)) {
        stop(paste0("Please check your selection of comparisons. There are ",
                    length(levels(ptm.model$Label)), 
                    " comparisons in this result."))
      }
    }  
    
    ## use only assigned comparisons
    ptm.model = ptm.model[which(ptm.model$Label %in% temp.name), ]
    if (nrow(protein.model) > 0){
      protein.model = protein.model[which(
        protein.model$Label %in% temp.name), ]
      adjusted.model = adjusted.model[which(
        adjusted.model$Label %in% temp.name),]
    }

  }
  
  if (which.Protein[[1]] != "all") {
    ## check which.comparison is name of comparison
    if (is.character(which.Protein)) {
      
      temp.name = which.Protein
      ## message if name of comparison is wrong.
      if (length(setdiff(temp.name, unique(ptm.model$Protein))) > 0) {
        
        ## TODO: Logging        
        # processout = rbind(processout, 
        #                     paste0("Please check labels of comparions. ", 
        #                            "Result does not have this comparison. - ", 
        #                            toString(temp.name)))
        # write.table(processout, file=finalfile, row.names=FALSE)
        
        stop(paste0("Please check labels of proteins ", 
                    "Result does not have all these proteins. - ", 
                    toString(temp.name)))
      }
    }
    ## check which.comparison is order number of comparison
    if (is.numeric(which.Protein)) {
      
      temp.name = levels(ptm.model$Protein)[which.Protein]
      
      ## message if name of comparison is wrong.
      if (length(levels(ptm.model$Protein)) < max(which.Protein)) {
        stop(paste0("Please check your selection of proteins. There are ",
                    length(levels(ptm.model$Protein)), 
                    " comparisons in this result."))
      }
    }
    
    ## use only assigned comparisons
    ptm.model = ptm.model[which(ptm.model$Protein %in% temp.name), ]
    if (nrow(protein.model) > 0){
      adjusted.model = adjusted.model[
        which(adjusted.model$Protein %in% temp.name), ]
      protein.model = protein.model[which(protein.model$Protein %in% unique(
        adjusted.model$GlobalProtein)), ]
    }
  }
  
  ptm.model$Protein = factor(ptm.model$Protein)
  ptm.model$Label = factor(ptm.model$Label)
  
  if (nrow(protein.model) > 0){
    protein.model$Protein = factor(protein.model$Protein)
    protein.model$Label = factor(protein.model$Label)
    adjusted.model$Protein = factor(adjusted.model$Protein)
    adjusted.model$Label = factor(adjusted.model$Label)
  }
  
  return(list('PTM.Model' = ptm.model, 'PROTEIN.Model' = protein.model, 
              'ADJUSTED.Model' = adjusted.model))
}

#' Wrapper for plotting heatmap for every model provided
#' @noRd
.plotHeatmap = function(data, sig, FCcutoff, logBase.pvalue, ylimUp, ylimDown, 
                     text.angle, x.axis.size, y.axis.size, dot.size, colorkey, 
                     numProtein, width, height, address, isPlotly){
  ## If there are the file with the same name, 
  ## add next numbering at the end of file name
  if (!isPlotly && address != FALSE) {
    allfiles = list.files()
    
    num = 0
    filenaming = paste0(address, "Heatmap")
    finalfile = paste0(address, "Heatmap.pdf")
    
    while (is.element(finalfile, allfiles)) {
      num = num + 1
      finalfile = paste0(paste(filenaming, num, sep="-"), ".pdf")
    }	
    
    pdf(finalfile, width=width, height=height)
  }
  final_plots = list()
  ptm_heatmap_plot <- .plot.model.heatmap(data[['PTM.Model']], sig, FCcutoff, logBase.pvalue, 
                      ylimUp, ylimDown, x.axis.size, y.axis.size, text.angle,
                      colorkey, numProtein,
                      "Unadjusted PTM", isPlotly)
  final_plots[[1]] <- ptm_heatmap_plot
  if (nrow(data[['PROTEIN.Model']]) > 0){
    protein_heatmap_plot <- .plot.model.heatmap(data[['PROTEIN.Model']], sig, FCcutoff, logBase.pvalue,
                        ylimUp, ylimDown, x.axis.size, y.axis.size, text.angle,
                        colorkey, numProtein,
                        "Protein", isPlotly)
    final_plots[[2]] <- protein_heatmap_plot
    adjusted_heatmap_plot <- .plot.model.heatmap(data[['ADJUSTED.Model']], sig, FCcutoff, logBase.pvalue,
                        ylimUp, ylimDown, x.axis.size, y.axis.size, text.angle,
                        colorkey, numProtein,
                        "Adjusted PTM", isPlotly)
    final_plots[[3]] <- adjusted_heatmap_plot
  }
  
  if (address != FALSE && !isPlotly) {
    dev.off()
  }
  if (isPlotly) {
    final_plots
  }
}



#' Plots individual model heatmap
#' @noRd
.plot.model.heatmap = function(data, sig, FCcutoff, logBase.pvalue, ylimUp,
                                 ylimDown, x.axis.size, y.axis.size, text.angle,
                                 colorkey, numProtein, model, isPlotly){

  Label = Protein = sign_adj_pval = NULL

  if (logBase.pvalue == 2) {
    y.limUp = 30
  } else if (logBase.pvalue == 10) {
    y.limUp = 10
  }

  if (is.numeric(ylimUp)) {
    y.limUp = ylimUp
  }


  ## if FCcutoff is assigned, make p-value insignificant.
  if (is.numeric(FCcutoff)) {
    if (colnames(data)[3] == "log2FC") {
      data$adj.pvalue[data[, 3] < log2(FCcutoff) & data[, 3] > (
        -log2(FCcutoff))] = 1
    }
    if (colnames(data)[3] == "log10FC") {
      data$adj.pvalue[data[, 3] < log10(FCcutoff) & data[, 3] > (
        -log10(FCcutoff))] = 1
    }
  }

  ## based on p-value
  if (logBase.pvalue == 2) {
      temp =  -log2(data$adj.pvalue) * sign(data[, 3])
  } else if (logBase.pvalue == 10) {
      temp =  -log10(data$adj.pvalue) * sign(data[, 3])
  }

  data$sign_adj_pval = temp
  
  obj = data[, c("Protein", "Label", "sign_adj_pval")]
  
  if(!isPlotly) {
    ## maximum number of proteins per heatmap
    unique_prot = unique(obj$Protein)
    totalpro = length(unique_prot)
    if(is.null(numProtein)) {
      numProtein = 50
    }
    numheatmap = totalpro %/% numProtein + 1
    ## draw heatmap
    ## loop for numProtein
    for (j in seq_len(numheatmap)) {
      
      plot_prot = unique_prot[((j - 1) * numProtein + 1):(j * numProtein)]
      temp_obj = obj[obj$Protein %in% plot_prot]
      
      limits = c(-1,1)*max(abs(temp_obj$sign_adj_pval[
        is.finite(temp_obj$sign_adj_pval)]))
      temp_heatmap = ggplot(temp_obj, aes(.data$Label, .data$Protein, fill = .data$sign_adj_pval)
      ) + geom_tile() + scale_fill_distiller(
        palette = "RdBu", name = "(sign) Adj pvalue",
        limits = limits) + labs(
          title = paste0(model, " - Page #", as.character(j)), y = model,
          x = "Comparison") +
        theme_msstats(type = "COMPARISONPLOT", x.axis.size, y.axis.size, 13,
                      element_rect(fill = "gray95"),
                      element_text(colour = c("#00B0F6"), size = 14),
                      "bottom", text_angle = text.angle)
      print(temp_heatmap)
    }
  } else {
    limits = c(-1, 1) * max(abs(obj$sign_adj_pval[is.finite(obj$sign_adj_pval)]))
    data$sign_adj_pval[is.infinite(data$sign_adj_pval)] <- NaN
    obj = data[, c("Protein", "Label", "sign_adj_pval")]
    if (!is.null(numProtein)) {
      # subset numProtein unique proteins
      obj$Protein <- as.character(obj$Protein)
      unique_proteins <- unique(obj$Protein)
      top_10_proteins <- unique_proteins[1:numProtein]
      obj <- obj[obj$Protein %in% top_10_proteins, ]
      limits = c(-1, 1) * max(abs(obj$sign_adj_pval[is.finite(obj$sign_adj_pval)]))
    }
    
    colors <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
    heatmap_plot <- plot_ly(
      data = obj,
      x = as.character(obj$Label),
      y = obj$Protein,
      z = ~sign_adj_pval,
      width = 800,
      height = 600,
      zmin = limits[1],
      zmax = limits[2],
      xgap = 0,
      ygap = 0,
      colorbar = list(title = "(sign) Adj pvalue", tickfont = list(size = 13)),
      colors = colors,
      zauto = FALSE,
      type = "heatmap"
    )
    
    # Customize the layout
    heatmap_plot <- plotly::layout(
      heatmap_plot,
      title = paste0(model, " : All Pages"),
      xaxis = list(title = "Comparison", tickangle = 45, tickfont = list(size = x.axis.size)),
      yaxis = list(title = model, tickfont = list(size = y.axis.size)),
      plot_bgcolor = "grey"
    )
    return(heatmap_plot)
  }

  
}

#' Wrapper for plotting volcano plot for every model provided
#' @noRd
.plotVolcano = function(data, sig, FCcutoff, logBase.pvalue, ylimUp, ylimDown, 
                     xlimUp, x.axis.size, y.axis.size, dot.size, text.size, 
                     legend.size, ProteinName, colorkey, numProtein,
                     width, height, address, plot_name_list, isPlotly){
  ## If there are the file with the same name
  ## add next numbering at the end of file name
  if (!isPlotly && address != FALSE) {
    allfiles = list.files()
    
    num = 0
    filenaming = paste0(address, "VolcanoPlot")
    finalfile = paste0(address, "VolcanoPlot.pdf")
    
    while (is.element(finalfile, allfiles)) {
      num = num + 1
      finalfile = paste0(paste(filenaming, num, sep="-"), ".pdf")
    }	
    
    pdf(finalfile, width=width, height=height)
  }

  ptm_plots <- .plot.model.volcano(data[['PTM.Model']], sig, FCcutoff, logBase.pvalue, 
                      ylimUp, ylimDown, xlimUp, x.axis.size, y.axis.size, 
                      dot.size, text.size, legend.size, ProteinName,
                      plot_name_list[[1]], isPlotly)
  final_plots = ptm_plots
  if (nrow(data[['PROTEIN.Model']])) {
    protein_plots <- .plot.model.volcano(data[['PROTEIN.Model']], sig, FCcutoff, logBase.pvalue, 
                        ylimUp, ylimDown, xlimUp, x.axis.size, y.axis.size, 
                        dot.size, text.size, legend.size, ProteinName,
                        plot_name_list[[2]], isPlotly)
    adjusted_plots <- .plot.model.volcano(data[['ADJUSTED.Model']], sig, FCcutoff, logBase.pvalue,
                        ylimUp, ylimDown, xlimUp, x.axis.size, y.axis.size, 
                        dot.size, text.size, legend.size, ProteinName,
                        plot_name_list[[3]], isPlotly)
    final_plots <- c(final_plots, protein_plots, adjusted_plots)
  }
  if (address != FALSE) {
    dev.off()
  }
  if (isPlotly) {
    final_plots
  }
}

#' Plots individual model volcano plot
#' @noRd
.plot.model.volcano = function(data, sig, FCcutoff, logBase.pvalue, ylimUp,
                                ylimDown, xlimUp, x.axis.size, y.axis.size, 
                                dot.size, text.size, legend.size, ProteinName,
                                model, isPlotly){
  log2FC = Protein = NULL
  
  if (logBase.pvalue == 2) {
    y.limUp  = 30
  } else if (logBase.pvalue == 10) {
    y.limUp  = 10
  }
  
  if (is.numeric(ylimUp)) {
    y.limUp = ylimUp 
  }
  
  ## remove the result, NA and Inf
  data = data[!is.na(data$adj.pvalue), ]
  data = data[is.finite(data$log2FC)]
  
  ## group for coloring dots
  if (!FCcutoff) {  
    data[data$adj.pvalue >= sig, "colgroup"] = "black"
    data[which(data$adj.pvalue < sig & data[, 3] > 0), "colgroup"] = "red"
    data[which(data$adj.pvalue < sig & data[, 3] < 0), "colgroup"] = "blue" 
  }
  
  if (is.numeric(FCcutoff)) {
    data$colgroup = "black"
    
    if (colnames(data)[3] == "log2FC") {
      data[which(data$adj.pvalue < sig & data[, 3] > log2(FCcutoff)), 
           "colgroup"] = "red"
      data[which(data$adj.pvalue < sig & data[, 3] < (-log2(FCcutoff))), 
           "colgroup"] = "blue"
    }
    
    if (colnames(data)[3] == "log10FC") {
      data = data[is.finite(data$log10FC)]
      data[which(data$adj.pvalue < sig & data[, 3] > log10(FCcutoff)), 
           "colgroup"] = "red"
      data[which(data$adj.pvalue < sig & data[, 3] < (-log10(FCcutoff))), 
           "colgroup"] = "blue"
    }
  }
  
  data$colgroup = factor(data$colgroup, levels=c("black", "blue", "red"))
  plots <- vector("list", nlevels(data$Label))
  ## for multiple volcano plots, 
  for (i in seq_len(nlevels(data$Label))) {
    
    sub = data[data$Label == levels(data$Label)[i], ]
    
    if (logBase.pvalue == 2) {
      sub$adj.pvalue[sub$adj.pvalue < 2^(-y.limUp)] = 2^(-y.limUp)
    } else if (logBase.pvalue == 10) {
      sub$adj.pvalue[sub$adj.pvalue < 10^(-y.limUp)] = 10^(-y.limUp)
    }
    
    sub = as.data.table(sub)
    
    if (logBase.pvalue == 2) {
      y.limup = ceiling(max(-log2(sub[!is.na(sub$adj.pvalue), "adj.pvalue"])))
      if (y.limup < (-log2(sig))) {
        y.limup = (-log2(sig) + 1) ## for too small y.lim
      }
    } else if (logBase.pvalue == 10) {
      y.limup = ceiling(max(-log10(sub[!is.na(sub$adj.pvalue), "adj.pvalue"])))
      if (y.limup < (-log10(sig))) {
        y.limup = (-log10(sig) + 1) ## for too small y.lim
      }
    }
    
    ## ylimDown
    y.limdown = 0 ## default is zero
    if (is.numeric(ylimDown)) {
      y.limdown = ylimDown
    }
    
    ## x.lim
    x.lim = ceiling(max(abs(sub[is.finite(log2FC), log2FC])))
    
    if (x.lim < 3) {
      x.lim = 3
    }
    if (is.numeric(xlimUp)) {
      x.lim = xlimUp
    }
    ## for assigning x in ggplot2
    subtemp = sub
    colnames(subtemp)[3] = "logFC"
    
    if (logBase.pvalue == 2) {
      subtemp$logadjp = (-log2(subtemp$adj.pvalue))
    } else if (logBase.pvalue == 10) {
      subtemp$logadjp = (-log10(subtemp$adj.pvalue))
    }
    
    ## for x limit for inf or -inf
    subtemp$newlogFC = subtemp$logFC
    subtemp[!is.na(subtemp$issue) &
              subtemp$issue == "oneConditionMissing" & 
              subtemp$logFC == Inf, "newlogFC"] = (x.lim - 0.2)
    subtemp[!is.na(subtemp$issue) & 
              subtemp$issue == "oneConditionMissing" & 
              subtemp$logFC == (-Inf), "newlogFC"] = (x.lim - 0.2) * (-1)
    
    ## add (*) in Protein name for Inf or -Inf
    subtemp$Protein = as.character(subtemp$Protein)
    subtemp[!is.na(subtemp$issue) & 
              subtemp$issue == "oneConditionMissing", "Protein"] = 
      paste0("*", subtemp[!is.na(subtemp$issue) & 
                            subtemp$issue == "oneConditionMissing", "Protein"])
    
    if (logBase.pvalue == 2) {
      log_title = 'Log2'
    } else if (logBase.pvalue == 10) {
      log_title = 'Log10'
    }
    
    ptemp = ggplot(aes(x = .data$logFC, y = .data$logadjp,
                               color = .data$colgroup,
                               label = .data$Protein),
                    data=subtemp) +
      geom_point(size=dot.size) +
      scale_colour_manual(values=c("gray65", "blue", "red"), 
                          limits=c("black", "blue", "red"), 
                          breaks=c("black", "blue", "red"), 
                          labels=c("No regulation", "Down-regulated", 
                                   "Up-regulated")) +
      scale_y_continuous(paste0('-', log_title, ' (adjusted p-value)'), 
                         limits=c(y.limdown, y.limup)) + 
      scale_x_continuous(paste0(log_title, ' fold change'), 
                         limits=c(-x.lim, x.lim)) +
      labs(title= paste0(model, ": ", unique(sub$Label)))
    
    ## add protein name
    if (ProteinName) {
      if (length(unique(subtemp$colgroup)) == 1 & 
          any(unique(subtemp$colgroup) == 'black')) {
        message(paste0("The volcano plot for ", unique(subtemp$Label), 
      " does not show the protein names because none of them are significant."))
      } else {
        ptemp = ptemp +
          geom_text_repel(data=subtemp[subtemp$colgroup != "black", ],
                          aes(label = .data$Protein),
                          size=text.size,
                          col='black')
      }
    }
    
    ## For legend of linetype for cutoffs
    ## first assign line type
    ltypes = c("type1"="twodash", "type2"="dotted")
    
    ## cutoff lines, FDR only
    if (!FCcutoff) {
      if (logBase.pvalue == 2){
        sigcut = data.table("Protein"='sigline', 
                             "logFC"=seq(-x.lim, x.lim, length.out=20),
                             "logadjp"=(-log2(sig)),
                             "line"='twodash')
      } else if (logBase.pvalue == 10) {
        sigcut = data.table("Protein"='sigline', 
                             "logFC"=seq(-x.lim, x.lim, length.out=20),
                             "logadjp"=(-log10(sig)),
                             "line"='twodash')
      }
      
      pfinal = ptemp +
        geom_line(data=sigcut,
                  aes(x = .data$logFC, y = .data$logadjp, linetype = .data$line),
                  colour="darkgrey",
                  linewidth=0.6,
                  show.legend=TRUE) +
        scale_linetype_manual(values=c('twodash'=6),
                              labels=c(paste0("Adj p-value cutoff (", sig, ")"))
                              ) +
        guides(colour=guide_legend(override.aes=list(linetype=0)),
               linetype=guide_legend())
    }
    if (is.numeric(FCcutoff)) {
      if (logBase.pvalue == 2) {
        ## three different lines
        sigcut = data.table("Protein"='sigline', 
                             "logFC"=seq(-x.lim, x.lim, length.out=10), 
                             "logadjp"=(-log2(sig)),
                             "line"='twodash')
        FCcutpos = data.table("Protein"='sigline', 
                               "logFC"=log2(FCcutoff), 
                               "logadjp"=seq(y.limdown, y.limup, length.out=10), 
                               "line"='dotted')
        FCcutneg = data.table("Protein"='sigline', 
                               "logFC"=(-log2(FCcutoff)), 
                               "logadjp"=seq(y.limdown, y.limup, length.out=10), 
                               "line"='dotted')
      } else if (logBase.pvalue == 10) {
        ## three different lines
        sigcut = data.table("Protein"='sigline', 
                             "logFC"=seq(-x.lim, x.lim, length.out=10), 
                             "logadjp"=(-log10(sig)),
                             "line"='twodash')
        FCcutpos = data.table("Protein"='sigline', 
                               "logFC"=log10(FCcutoff), 
                               "logadjp"=seq(y.limdown, y.limup, length.out=10), 
                               "line"='dotted')
        FCcutneg = data.table("Protein"='sigline', 
                               "logFC"=(-log10(FCcutoff)), 
                               "logadjp"=seq(y.limdown, y.limup, length.out=10), 
                               "line"='dotted')
      }
      ## three lines, with order color first and then assign linetype manual
      pfinal = ptemp +
        geom_line(data=sigcut, 
                  aes(x = .data$logFC, y = .data$logadjp, linetype = .data$line),
                  colour="darkgrey",
                  linewidth=0.6,
                  show.legend=TRUE) +
        geom_line(data=FCcutpos,
                  aes(x = .data$logFC, y = .data$logadjp, linetype = .data$line),
                  colour="darkgrey",
                  linewidth=0.6,
                  show.legend=TRUE) +
        geom_line(data=FCcutneg,
                  aes(x = .data$logFC, y = .data$logadjp, linetype = .data$line),
                  colour="darkgrey",
                  linewidth=0.6) +
        scale_linetype_manual(values=c('dotted'=3, 'twodash'=6),
                              labels=c(paste0("Fold change cutoff (", FCcutoff, 
                                              ")"),
                                       paste0("Adj p-value cutoff (", sig, ")"))
                              ) +
        guides(colour=guide_legend(override.aes=list(linetype=0)),
               linetype=guide_legend())
    }
    
    pfinal = pfinal +
      theme_msstats(type = "VOLCANO", x.axis.size, y.axis.size, legend.size, 
                    element_rect(fill = "gray95"),
                    element_text(colour = c("#00B0F6"), size = 14),
                    "bottom", text_angle = text.angle, 
                    legend.title = element_blank())
      # theme(
      #   panel.background = element_rect(fill='white', colour="black"),
      #   panel.grid.minor = element_blank(),
      #   axis.text.x = element_text(size=x.axis.size, colour="black"),
      #   axis.text.y = element_text(size=y.axis.size, colour="black"),
      #   axis.ticks = element_line(colour="black"),
      #   axis.title.x = element_text(size=x.axis.size+5, vjust=-0.4),
      #   axis.title.y = element_text(size=y.axis.size+5, vjust=0.3),
      #   title = element_text(size=x.axis.size+8, vjust=1.5),
      #   legend.position="bottom",
      #   legend.key = element_rect(fill='white', colour='white'),
      #   legend.text = element_text(size=legend.size),
      #   legend.title = element_blank()
      # )
    
    print(pfinal)
    plots[[i]] = pfinal
  }
  if (isPlotly) {
    plots
  }
}