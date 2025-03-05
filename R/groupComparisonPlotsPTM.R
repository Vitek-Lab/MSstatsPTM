#' Visualization for model-based analysis and summarization
#'
#' To analyze the results of modeling changes in abundance of modified peptides 
#' and overall protein, groupComparisonPlotsPTM takes as input the results of 
#' the groupComparisonPTM function. It asses the results of three models: 
#' unadjusted PTM, adjusted PTM, and overall protein. To asses the results of
#' the model, the following visualizations can be created:
#' (1) VolcanoPlot (specify "VolcanoPlot" in option type), to plot peptides or
#' proteins and their significance for each model.
#' (2) Heatmap (specify "Heatmap" in option type), to evaluate the fold change
#' between conditions and peptides/proteins
#' 
#' @export
#' @import ggplot2
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom data.table data.table as.data.table 
#' @importFrom ggrepel geom_text_repel
#' @importFrom checkmate assertNumeric assertChoice assertLogical
#' @importFrom MSstats theme_msstats
#' @importFrom plotly ggplotly style add_trace plot_ly subplot layout
#' 
#' @param data name of the list with models, which can be the output of the 
#' MSstatsPTM \code{\link[MSstatsPTM]{groupComparisonPTM}} function
#' @param type choice of visualization, one of VolcanoPlot or Heatmap
#' @param sig FDR cutoff for the adjusted p-values in heatmap and volcano plot. 
#' level of significance for comparison plot. 100(1-sig)% confidence interval 
#' will be drawn. sig=0.05 is default.
#' @param FCcutoff For volcano plot or heatmap, whether involve fold change 
#' cutoff or not. FALSE (default) means no fold change cutoff is applied for 
#' significance analysis. FCcutoff = specific value means specific fold change 
#' cutoff is applied.
#' @param logBase.pvalue for volcano plot or heatmap, (-) logarithm 
#' transformation of adjusted p-value with base 2 or 10(default).
#' @param ylimUp for all three plots, upper limit for y-axis. FALSE (default) 
#' for volcano plot/heatmap use maximum of -log2 (adjusted p-value) or -log10 
#' (adjusted p-value). FALSE (default) for comparison plot uses maximum of 
#' log-fold change + CI.
#' @param ylimDown for all three plots, lower limit for y-axis. FALSE (default) 
#' for volcano plot/heatmap use minimum of -log2 (adjusted p-value) or -log10 
#' (adjusted p-value). FALSE (default) for comparison plot uses minimum of 
#' log-fold change - CI.
#' @param xlimUp for Volcano plot, the limit for x-axis. FALSE (default) for 
#' use maximum for absolute value of log-fold change or 3 as default if maximum 
#' for absolute value of log-fold change is less than 3.
#' @param x.axis.size size of axes labels, e.g. name of the comparisons in 
#' heatmap, and in comparison plot. Default is 10.
#' @param y.axis.size size of axes labels, e.g. name of targeted proteins in 
#' heatmap. Default is 10.
#' @param dot.size size of dots in volcano plot and comparison plot. Default is 
#' 3.
#' @param text.size size of ProteinName label in the graph for Volcano Plot. 
#' Default is 4.
#' @param text.angle angle of x-axis labels represented each comparison at the 
#' bottom of graph in comparison plot. Default is 0.
#' @param legend.size size of legend for color at the bottom of volcano plot. 
#' Default is 7.
#' @param ProteinName for volcano plot only, whether display protein names or 
#' not. TRUE (default) means protein names, which are significant, are 
#' displayed next to the points. FALSE means no protein names are displayed.
#' @param colorkey TRUE(default) shows colorkey.
#' @param numProtein The number of proteins which will be presented in each 
#' heatmap. Default is 50.
#' @param width width of the saved file. Default is 10.
#' @param height height of the saved file. Default is 10.
#' @param which.Comparison list of comparisons to draw plots. List can be 
#' labels of comparisons or order numbers of comparisons from levels(data$Label)
#' , such as levels(testResultMultiComparisons$ComparisonResult$Label). 
#' Default is "all", which generates all plots for each protein.
#' @param which.PTM Protein list to draw comparison plots. List can be 
#' names of Proteins or order numbers of Proteins from 
#' levels(testResultMultiComparisons$ComparisonResult$Protein). Default is 
#' "all", which generates all comparison plots for each protein.
#' @param address the name of folder that will store the results. Default 
#' folder is the current working directory. The other assigned folder has to 
#' be existed under the current working directory. An output pdf file is 
#' automatically created with the default name of "VolcanoPlot.pdf" or 
#' "Heatmap.pdf". The command address can help to specify where to store the 
#' file as well as how to modify the beginning of the file name. If 
#' address=FALSE, plot will be not saved as pdf file but showed in window
#' @param isPlotly Parameter to use Plotly or ggplot2. If set to TRUE, MSstats 
#' will save Plotly plots as HTML files. If set to FALSE MSstats will save ggplot2 plots
#' as PDF files
#' @return plot or pdf
#' @examples 
#' 
#' model.lf.msstatsptm = groupComparisonPTM(summary.data, 
#'                                      data.type = "LabelFree")
#' groupComparisonPlotsPTM(data = model.lf.msstatsptm,
#'                         type = "VolcanoPlot",
#'                         FCcutoff= 2,
#'                         logBase.pvalue = 2,
#'                         address=FALSE)
groupComparisonPlotsPTM = function(data = data,
                                    type,
                                    sig=0.05,
                                    FCcutoff=FALSE,
                                    logBase.pvalue=10,
                                    ylimUp=FALSE,
                                    ylimDown=FALSE,
                                    xlimUp=FALSE,
                                    x.axis.size=10,
                                    y.axis.size=10,
                                    dot.size=3,
                                    text.size=4,
                                    text.angle=0,
                                    legend.size=13,
                                    ProteinName=TRUE,
                                    colorkey=TRUE,
                                    numProtein=NULL, 
                                    width=10, 
                                    height=10, 
                                    which.Comparison = "all", 
                                    which.PTM = "all",
                                    address="",
                                    isPlotly = FALSE) {
  
  .check.plotting.data(data, type, sig, FCcutoff, logBase.pvalue, ylimUp,
                       ylimDown, xlimUp, x.axis.size, y.axis.size, dot.size,
                       text.size, text.angle, legend.size, ProteinName,
                       colorkey, numProtein, width, height, which.Comparison, 
                       which.PTM, address)
  
  data = .format.model.plots(data, which.Comparison, which.PTM)
  type = toupper(type)
  
  plot_name_list = c("Unadjusted Peptides", "Global Protein", 
                     "Adjusted Peptide")
  
  if(isPlotly & address != FALSE) {
    print("Plots will be saved as .HTML file as plotly is selected, set isPlotly = FALSE, if 
            you want to generate PDF using ggplot2")
  }
  
  if (type == 'HEATMAP') {
    plotly_plots <- .plotHeatmap(data, sig, FCcutoff, logBase.pvalue, ylimUp, ylimDown, text.angle,
                  x.axis.size, y.axis.size, dot.size, colorkey, numProtein, 
                  width, height, address, isPlotly)
    
    if(isPlotly) {
      if(address != FALSE) {
        .savePlotlyPlotHTML(plotly_plots,address,"Heatmap" ,width, height)
      }
      plotly_plots
    }
    
  } else if (type == 'VOLCANOPLOT') {
    plots <- .plotVolcano(data, sig, FCcutoff, logBase.pvalue, ylimUp, ylimDown, xlimUp, 
                  x.axis.size, y.axis.size, dot.size, text.size, legend.size, 
                  ProteinName, colorkey, numProtein, width, height, 
                  address, plot_name_list, isPlotly)
    
    plotly_plots <- vector("list", length(plots))
    if(isPlotly) {
      for(i in seq_along(plots)) {
        plot <- plots[[i]]
        plotly_plot <- .convertGgplot2Plotly(plot, width = 1000)
        plotly_plot <- .fixLegendPlotlyPlotsVolcano(plotly_plot)
        plotly_plots[[i]] = list(plotly_plot)
      }
      if(address != FALSE) {
        .savePlotlyPlotHTML(plotly_plots,address,"VolcanoPlot" ,width, height)
      }
      plotly_plots <- unlist(plotly_plots, recursive = FALSE)
      plotly_plots
    }
  }
}

.fixLegendPlotlyPlotsVolcano = function(plot) {
  df <- data.frame(id = seq_along(plot$x$data), legend_entries = unlist(lapply(plot$x$data, `[[`, "name")))
  # Create a mapping
  color_mapping <- c("black" = "No regulation", "red" = "Up-regulated", "blue" = "Down-regulated")
  # Update the legend_entries column
  df$legend_group <- sapply(df$legend_entries, function(entry) {
    for (color in names(color_mapping)) {
      if (grepl(color, entry)) {
        entry <- gsub(color, color_mapping[color], entry)
        break
      }
    }
    entry <- gsub(",.+", "", entry)
    entry <- gsub("\\(|\\)", "", entry)  # Remove any remaining parentheses
    entry
  })
  for (i in df$id) {
    if(length(grep(df$legend_group[[i]], color_mapping)) == 0) { # keep only 3 legends
      plot$x$data[[i]]$showlegend <- FALSE
    }
    plot$x$data[[i]]$name <- df$legend_group[[i]]
    plot$x$data[[i]]$legendgroup <- plot$x$data[[i]]$name
  }
  plot <- plotly::layout(plot,legend=list(title=list(text="")))
  plot
}
  
  