data("raw.input", package = "MSstatsPTM")

## Test missing columns are handled
expect_error(dataSummarizationPTM())
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM[, -1],
                                       PROTEIN = raw.input$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM,
                                       PROTEIN = raw.input$PROTEIN[, -1])))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM[, -2],
                                       PROTEIN = raw.input$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM,
                                       PROTEIN = raw.input$PROTEIN[, -2])))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM[, -3],
                                       PROTEIN = raw.input$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM,
                                       PROTEIN = raw.input$PROTEIN[, -3])))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM[, -4],
                                       PROTEIN = raw.input$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM,
                                       PROTEIN = raw.input$PROTEIN[, -4])))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM[, -5],
                                       PROTEIN = raw.input$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM,
                                       PROTEIN = raw.input$PROTEIN[, -5])))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM[, -6],
                                       PROTEIN = raw.input$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM,
                                       PROTEIN = raw.input$PROTEIN[, -6])))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM[, -7],
                                       PROTEIN = raw.input$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM,
                                       PROTEIN = raw.input$PROTEIN[, -7])))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM[, -8],
                                       PROTEIN = raw.input$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM,
                                       PROTEIN = raw.input$PROTEIN[, -8])))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM[, -9],
                                       PROTEIN = raw.input$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM,
                                       PROTEIN = raw.input$PROTEIN[, -9])))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM[, -10],
                                       PROTEIN = raw.input$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input$PTM,
                                       PROTEIN = raw.input$PROTEIN[, -10])))

## Test misspecified parameters
expect_error(dataSummarizationPTM(raw.input, logTrans = 20))
expect_error(dataSummarizationPTM(raw.input, logTrans = "test"))
expect_error(dataSummarizationPTM(raw.input, normalization = TRUE))
expect_error(dataSummarizationPTM(raw.input, normalization.PTM = TRUE))
expect_error(dataSummarizationPTM(raw.input, fillIncompleteRows = 20))
expect_error(dataSummarizationPTM(raw.input, summaryMethod = TRUE))
expect_error(dataSummarizationPTM(raw.input, censoredInt = ""))
expect_error(dataSummarizationPTM(raw.input, cutoffCensored = ""))
expect_error(dataSummarizationPTM(raw.input, MBimpute = ""))
expect_error(dataSummarizationPTM(raw.input, MBimpute.PTM = ""))