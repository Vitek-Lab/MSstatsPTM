#' Visualization for explanatory data analysis
#'
#' To illustrate the quantitative data and quality control of MS runs,
#' dataProcessPlotsPTM takes the quantitative data from dataSummarizationPTM or
#' dataSummarizationPTM_TMT to plot the following :
#' (1) profile plot (specify "ProfilePlot" in option type), to identify the
#' potential sources of variation for each protein;
#' (2) quality control plot (specify "QCPlot" in option type), to evaluate the
#' systematic bias between MS runs.
#'
#' @export
#' @import ggplot2
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom gridExtra grid.arrange
#' @importFrom stringr str_match
#' @importFrom data.table data.table as.data.table melt dcast `:=` setnames copy rbindlist
#' @importFrom checkmate assertCharacter assertNumeric assertChoice
#' @importFrom MSstats theme_msstats
#' @importFrom plotly ggplotly style add_trace plot_ly subplot layout
#' 
#' @param data name of the list with PTM and (optionally) Protein data, which
#' can be the output of the MSstatsPTM 
#' \code{\link[MSstatsPTM]{dataSummarizationPTM}} or 
#' \code{\link[MSstatsPTM]{dataSummarizationPTM_TMT}} functions.
#' @param type choice of visualization. "ProfilePlot" represents profile plot of
#'  log intensities across MS runs.
#' "QCPlot" represents box plots of log intensities across channels and MS runs.
#' @param ylimUp upper limit for y-axis in the log scale.
#' FALSE(Default) for Profile Plot and QC Plot uses the upper limit as rounded
#' off maximum of log2(intensities) after normalization + 3..
#' @param ylimDown lower limit for y-axis in the log scale. FALSE(Default) for
#' Profile Plot and QC Plot uses 0..
#' @param x.axis.size size of x-axis labeling for "Run" and "channel in Profile
#' Plot and QC Plot.
#' @param y.axis.size size of y-axis labels. Default is 10.
#' @param text.size size of labels represented each condition at the top of
#' Profile plot and QC plot. Default is 4.
#' @param text.angle angle of labels represented each condition at the top of
#' Profile plot and QC plot. Default is 0.
#' @param legend.size size of legend above Profile plot. Default is 7.
#' @param dot.size.profile size of dots in Profile plot. Default is 2.
#' @param ncol.guide number of columns for legends at the top of plot. Default
#' is 5.
#' @param width width of the saved pdf file. Default is 10.
#' @param height height of the saved pdf file. Default is 10.
#' @param ptm.title title of overall PTM QC plot
#' @param protein.title title of overall Protein QC plot
#' @param which.PTM PTM list to draw plots. List can be names of
#' PTMs or order numbers of PTMs.
#' Default is "all", which generates all plots for each protein. For QC plot,
#' "allonly" will generate one QC plot with all proteins.
#' @param which.Protein List of proteins to plot. Will plot all PTMs associated 
#' with listed Proteins. Default is NULL which will default to which.PTM.
#' @param originalPlot TRUE(default) draws original profile plots, without
#' normalization.
#' @param summaryPlot TRUE(default) draws profile plots with protein
#' summarization for each channel and MS run.
#' @param address the name of folder that will store the results. Default folder
#'  is the current working directory.
#' @param isPlotly Parameter to use Plotly or ggplot2. If set to TRUE, MSstats 
#' will save Plotly plots as HTML files. If set to FALSE MSstats will save ggplot2 plots
#' as PDF files
#' The other assigned folder has to be existed under the current working
#' directory.
#' An output pdf file is automatically created with the default name of
#' "ProfilePlot.pdf" or "QCplot.pdf".
#' The command address can help to specify where to store the file as well as
#' how to modify the beginning of the file name.
#' If address=FALSE, plot will be not saved as pdf file but showed in window.
#' @return plot or pdf
#' @examples
#' 
#' # QCPlot
#' dataProcessPlotsPTM(summary.data,
#'                     type = 'QCPLOT',
#'                     which.PTM = "allonly",
#'                     address = FALSE)
#'                     
#' #ProfilePlot
#' dataProcessPlotsPTM(summary.data,
#'                     type = 'PROFILEPLOT',
#'                     which.PTM = "Q9UQ80_K376",
#'                     address = FALSE)
dataProcessPlotsPTM = function(data,
                                type = 'PROFILEPLOT',
                                ylimUp = FALSE,
                                ylimDown = FALSE,
                                x.axis.size = 10,
                                y.axis.size = 10,
                                text.size = 4,
                                text.angle = 90,
                                legend.size = 7,
                                dot.size.profile = 2,
                                ncol.guide = 5,
                                width = 10,
                                height = 12,
                                ptm.title = "All PTMs",
                                protein.title = "All Proteins",
                                which.PTM = "all", 
                                which.Protein = NULL,
                                originalPlot = TRUE,
                                summaryPlot = TRUE,
                                address = "",
                                isPlotly = FALSE) {
  
  type = toupper(type)
  label = .check.dataProcess.plotting.data(data, type, ylimUp, ylimDown, 
                                            x.axis.size,
                                y.axis.size, text.size, text.angle, legend.size,
                                dot.size.profile, ncol.guide, width, height,
                                ptm.title, protein.title, which.PTM, 
                                originalPlot, summaryPlot, address)
  
  data.table.list = .format.data.process.plots(data, label)
  
  if(isPlotly & address != FALSE) {
    print("Plots will be saved as .HTML file as plotly is selected, set isPlotly = FALSE, if 
            you want to generate PDF using ggplot2")
  }
  
  ## Filter for all PTMs in one proteins
  if (!is.null(which.Protein)){
    
    data.table.list[[1]] = data.table.list[[1]][PROTEINNAME %in% which.Protein]
    data.table.list[[2]] = data.table.list[[2]][
      GLOBALPROTEIN %in% which.Protein]
    data.table.list[[3]] = data.table.list[[3]][PROTEINNAME %in% which.Protein]
    data.table.list[[4]] = data.table.list[[4]][
      GLOBALPROTEIN %in% which.Protein]
    
    if (sum(nrow(data.table.list[[1]]), nrow(data.table.list[[2]]), 
           nrow(data.table.list[[3]]), nrow(data.table.list[[4]])) == 0){
      msg = paste0("The protein ", which.Protein, " specified in Which.Protein",
                   " is not in the global protein run. Please specify ",
                   "individual peptides only for this Protein.")
      stop(msg)
    }
  }
  
  ## Profile plot ##
  ## ---------------
  if (type == "PROFILEPLOT") {
    if (label == 'TMT'){
      plots <- .profile.tmt(data.table.list, type, ylimUp, ylimDown, 
                   x.axis.size, y.axis.size,text.size,text.angle, 
                   legend.size, dot.size.profile, ncol.guide, width, 
                   height, which.PTM, originalPlot, summaryPlot, 
                   address, isPlotly)
    } else if (label == 'LabelFree'){
      plots <- .profile.lf(data.table.list, type, ylimUp, ylimDown, 
                  x.axis.size, y.axis.size,text.size,text.angle, 
                  legend.size, dot.size.profile, ncol.guide, width, 
                  height, which.PTM, originalPlot, 
                  summaryPlot, address, isPlotly)
    }
    
    plotly_plots <- list()
    if(isPlotly) {
      og_plotly_plot = NULL
      summ_plotly_plot = NULL
      if("original_plot" %in% names(plots)) {
        for(i in seq_along(plots[["original_plot"]])) {
          plot_i <- plots[["original_plot"]][[paste("plot",i)]]
          plotly_plot_ptm <- .convertGgplot2Plotly(plot_i[["PTEMP.PTM"]])
          plotly_plot_ptm = .fixLegendPlotlyPlotsDataprocess(plotly_plot_ptm, "OriginalPlot")
          if (!is.null(plot_i[["PTEMP.PROTEIN"]])) {
            plotly_plot_protein <- .convertGgplot2Plotly(plot_i[["PTEMP.PROTEIN"]])
            plotly_plot_protein = .fixLegendPlotlyPlotsDataprocess(plotly_plot_protein, "OriginalPlot")
            plotly_plot_combined_original <- .combineSubPlotsPlotly(plotly_plot_ptm, plotly_plot_protein)
            plotly_plots = c(plotly_plots, list(plotly_plot_combined_original))
          } else {
            plotly_plots = c(plotly_plots, list(plotly_plot_ptm))
          }
        }
      }
      if("summary_plot" %in% names(plots)) {
        for(i in seq_along(plots[["summary_plot"]])) {
          plot_i <- plots[["summary_plot"]][[paste("plot",i)]]
          plotly_plot_ptm <- .convertGgplot2Plotly(plot_i[["PTEMP.PTM"]])
          plotly_plot_ptm = .fixLegendPlotlyPlotsDataprocess(plotly_plot_ptm, "SummaryPlot")
          if (!is.null(plot_i[["PTEMP.PROTEIN"]])) {
            plotly_plot_protein <- .convertGgplot2Plotly(plot_i[["PTEMP.PROTEIN"]])
            plotly_plot_combined_summary <- .combineSubPlotsPlotly(plotly_plot_ptm, plotly_plot_protein)
            plotly_plot_combined_summary = .fixLegendPlotlyPlotsDataprocess(plotly_plot_combined_summary, "SummaryPlot")
            plotly_plots = c(plotly_plots, list(plotly_plot_combined_summary))
          } else {
            plotly_plots = c(plotly_plots, list(plotly_plot_ptm))
          }
        }
      }
      if(address != FALSE) {
        .savePlotlyPlotHTML(plotly_plots,address,"ProfilePlot" ,width, height)
      }
      plotly_plots
    }
  }

  ## QC plot (Quality control plot) ##
  ## ---------------------------------
  else if (type == "QCPLOT") {
    if (label == 'TMT'){
      plots <- .qc.tmt(data.table.list, type, ylimUp, ylimDown, width, height, 
              x.axis.size, y.axis.size, text.size, text.angle,
              which.PTM, address, ptm.title, protein.title, isPlotly)
    } else if (label == 'LabelFree'){
      plots <- .qc.lf(data.table.list, type, ylimUp, ylimDown, width, height, 
             x.axis.size, y.axis.size, text.size,
             which.PTM, address, ptm.title, protein.title, isPlotly)
    } 
    
    plotly_plots <- list()
    if(isPlotly) {
      for(i in seq_along(plots)) {
        plot <- plots[[i]]
        plotly_plot_ptm <- .convertGgplot2Plotly(plot[["PTEMP.PTM"]])
        if(label == "TMT") {
          plotly_plot_ptm <- facet_strip_bigger(plotly_plot_ptm)
        }
        if (!is.null(plot[["PTEMP.PROTEIN"]])) {
          plotly_plot_protein <- .convertGgplot2Plotly(plot[["PTEMP.PROTEIN"]])
          if(label == "TMT") {
            plotly_plot_protein <- facet_strip_bigger(plotly_plot_protein)
          }
          plotly_plot_combined <- .combineSubPlotsPlotly(plotly_plot_ptm, plotly_plot_protein)
          plotly_plots[[i]] = list(plotly_plot_combined)
        } else {
          plotly_plots[[i]] = list(plotly_plot_ptm)
        }
        
      }
      if(address != FALSE) {
        .savePlotlyPlotHTML(plotly_plots,address,"QCPlot" ,width, height)
      }
      plotly_plots <- unlist(plotly_plots, recursive = FALSE)
      plotly_plots
    }
  }
}

.combineSubPlotsPlotly = function(plotly_plot_ptm, plotly_plot_protein) {
  title_ptm <- plotly_plot_ptm$x$layout$title$text
  title_protein <- plotly_plot_protein$x$layout$title$text
  
  plotly_plot_ptm <- plotly::layout(plotly_plot_ptm, title = "")
  plotly_plot_protein <- plotly::layout(plotly_plot_protein, title = "")
  
  plotly_plot_combined <- subplot(plotly_plot_ptm, plotly_plot_protein, nrows = 2, margin=0.1,titleX = TRUE, titleY = TRUE)
  plotly_plot_combined <- plotly::layout(plotly_plot_combined,
                                         annotations = list(
                                           list(
                                             x = 0.04,  # Centered horizontally
                                             y = 1.05,  # Above the first plot
                                             text = title_ptm,  # Use extracted title
                                             showarrow = FALSE,
                                             xref = 'paper',
                                             yref = 'paper',
                                             xanchor = 'center',
                                             yanchor = 'bottom',
                                             font = list(size = 16)
                                           ),
                                           list(
                                             x = 0.04,  # Centered horizontally
                                             y = 0.45,  # Above the second plot (adjust as necessary)
                                             text = title_protein,  # Use extracted title
                                             showarrow = FALSE,
                                             xref = 'paper',
                                             yref = 'paper',
                                             xanchor = 'center',
                                             yanchor = 'bottom',
                                             font = list(size = 16)
                                           )
                                         )
  )
  plotly_plot_combined
}

.fixLegendPlotlyPlotsDataprocess = function(plot, type) {
  df <- data.frame(id = seq_along(plot$x$data), legend_entries = unlist(lapply(plot$x$data, `[[`, "name")))
  df$legend_group <- gsub("^\\((.*?),.*", "\\1", df$legend_entries)
  df$is_first <- !duplicated(df$legend_group)
  df$is_bool <- ifelse(grepl("TRUE|FALSE", df$legend_group), TRUE, FALSE)
  df$is_valid_column <- ifelse(grepl("Processed feature-level data|Run summary", df$legend_entries), TRUE, FALSE)
  plot$x$data[[nrow(df)]]$showlegend <- FALSE # remove text legend
  
  for (i in df$id) {
    is_first <- df$is_first[[i]]
    is_bool <- df$is_bool[[i]]
    plot$x$data[[i]]$name <- df$legend_group[[i]]
    plot$x$data[[i]]$legendgroup <- plot$x$data[[i]]$name
    if (!is_first) {
      plot$x$data[[i]]$showlegend <- FALSE
    } 
    if(type == "SummaryPlot") {
      is_valid_column <- df$is_valid_column[[i]]
      if (!is_valid_column) plot$x$data[[i]]$showlegend <- FALSE
    }
    if(is_bool) plot$x$data[[i]]$showlegend <- FALSE
  }
  plot
}

facet_strip_bigger <- function(gp){
  if (!is.null(gp$x$layout$annotations)) {
    for (i in seq_along(gp$x$layout$annotations)) {
      if(gp$x$layout$annotations[[i]]$text != "Log2-intensities" && gp$x$layout$annotations[[i]]$text != "MS runs") {
        gp$x$layout$annotations[[i]]$font$size <- 7
        # gp$x$layout$annotations[[i]]$xanchor <- "center"
        gp$x$layout$annotations[[i]]$xshift <- 50
      }
    }
  }
  return(gp)
}

#' converter for plots from ggplot to plotly
#' @noRd
.convertGgplot2Plotly = function(plot, tips = "all", width = 1800, height = 600) {
  converted_plot <- ggplotly(plot,tooltip = tips)
  converted_plot <- plotly::layout(
    converted_plot,
    width = width,   # Set the width of the chart in pixels
    height = height,  # Set the height of the chart in pixels
    title = list(
      font = list(
        size = 18
      )
    ),
    legend = list(
      x = 0,     # Set the x position of the legend
      y = -0.25,    # Set the y position of the legend (negative value to move below the plot)
      orientation = "h",  # Horizontal orientation
      font = list(
        size = 12  # Set the font size for legend item labels
      ),
      title = list(
        font = list(
          size = 12  # Set the font size for the legend title
        )
      )
    )
  ) 
  converted_plot
}

.savePlotlyPlotHTML = function(plots, address, file_name, width, height) {
  print("Saving plots as HTML")
  pb <- txtProgressBar(min = 0, max = 4, style = 3)
  
  setTxtProgressBar(pb, 1)
  file_name = getFileName(address, file_name, width, height)
  file_name = paste0(file_name,".html")
  
  setTxtProgressBar(pb, 2)
  doc <- .getPlotlyPlotHTML(plots, width, height)
  
  setTxtProgressBar(pb, 3)
  htmltools::save_html(html = doc, file = file_name) # works but lib same folder
  
  setTxtProgressBar(pb, 4)
  zip(paste0(gsub("\\.html$", "", file_name),".zip"), c(file_name, "lib"))
  unlink(file_name)
  unlink("lib",recursive = T)
  
  close(pb)
}

getFileName = function(name_base, file_name, width, height) {
  all_files = list.files(".")
  if(file_name == 'ProfilePlot'){
    num_same_name = sum(grepl(paste0("^", name_base, file_name, "_[0-9]?"), all_files))
  } else {
    num_same_name = sum(grepl(paste0("^", name_base, file_name, "[0-9]?"), all_files))
  }
  if (num_same_name > 0) {
    file_name = paste(file_name, num_same_name + 1, sep = "_")
  }
  file_path = paste0(name_base, file_name)
  return(file_path)
}

.getPlotlyPlotHTML = function(plots, width, height) {
  doc <- htmltools::tagList(lapply(plots,function(x) htmltools::div(x, style = "float:left;width:100%;")))
  # Set a specific width for each plot
  plot_width <- 800
  plot_height <- 600
  
  # Create a div for each plot with style settings
  divs <- lapply(plots, function(x) {
    htmltools::div(x, style = paste0("width:", plot_width, "px; height:", plot_height, "px; margin: 10px;"))
  })
  
  # Combine the divs into a tagList
  doc <- htmltools::tagList(divs)
  doc
}
