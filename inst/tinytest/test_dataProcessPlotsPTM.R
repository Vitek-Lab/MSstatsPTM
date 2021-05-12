
data("summary.data", package = "MSstatsPTM")
data("summary.data.tmt", package = "MSstatsPTM")

## Plot dataProcessPlots expected output
expect_warning(dataProcessPlotsPTM(summary.data, "ProfilePlot", 
                                  which.PTM = 1:10, address = FALSE))
# expect_silent(dataProcessPlotsPTM(summary.data.tmt, "ProfilePlot", 
#                                   which.PTM = 1:10, address = FALSE))
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   dataProcessPlotsPTM(data,
#                                 type = 'PROFILEPLOT',
#                                 ylimUp = FALSE,
#                                 ylimDown = FALSE,
#                                 x.axis.size = 10,
#                                 y.axis.size = 10,
#                                 text.size = 4,
#                                 text.angle = 90,
#                                 legend.size = 7,
#                                 dot.size.profile = 2,
#                                 ncol.guide = 5,
#                                 width = 10,
#                                 height = 12,
#                                 ptm.title = "All PTMs",
#                                 protein.title = "All Proteins",
#                                 which.PTM = "all", 
#                                 which.Protein = NULL,
#                                 originalPlot = TRUE,
#                                 summaryPlot = TRUE,
#                                 address = "")