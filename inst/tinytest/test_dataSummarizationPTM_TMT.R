data("raw.input.tmt", package = "MSstatsPTM")

## Test missing columns are handled
expect_error(dataSummarizationPTM_TMT())
expect_error(dataSummarizationPTM_TMT(list(PTM = raw.input.tmt$PTM[, -1],
                                       PROTEIN = raw.input.tmt$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM,
                                       PROTEIN = raw.input.tmt$PROTEIN[, -1])))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM[, -2],
                                       PROTEIN = raw.input.tmt$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM,
                                       PROTEIN = raw.input.tmt$PROTEIN[, -2])))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM[, -3],
                                       PROTEIN = raw.input.tmt$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM,
                                       PROTEIN = raw.input.tmt$PROTEIN[, -3])))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM[, -4],
                                       PROTEIN = raw.input.tmt$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM,
                                       PROTEIN = raw.input.tmt$PROTEIN[, -4])))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM[, -5],
                                       PROTEIN = raw.input.tmt$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM,
                                       PROTEIN = raw.input$PROTEIN[, -5])))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM[, -6],
                                       PROTEIN = raw.input.tmt$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM,
                                       PROTEIN = raw.input.tmt$PROTEIN[, -6])))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM[, -7],
                                       PROTEIN = raw.input.tmt$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM,
                                       PROTEIN = raw.input.tmt$PROTEIN[, -7])))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM[, -8],
                                       PROTEIN = raw.input.tmt$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM,
                                       PROTEIN = raw.input.tmt$PROTEIN[, -8])))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM[, -9],
                                       PROTEIN = raw.input.tmt$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM,
                                       PROTEIN = raw.input.tmt$PROTEIN[, -9])))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM[, -10],
                                       PROTEIN = raw.input.tmt$PROTEIN)))
expect_error(dataSummarizationPTM(list(PTM = raw.input.tmt$PTM,
                                       PROTEIN = raw.input.tmt$PROTEIN[, -10])))

## Test misspecified parameters
expect_error(dataSummarizationPTM(raw.input.tmt, method  = ""))
expect_error(dataSummarizationPTM(raw.input.tmt, global_norm = 12))
expect_error(dataSummarizationPTM(raw.input.tmt, global_norm.PTM = 12))
expect_error(dataSummarizationPTM(raw.input.tmt, reference_norm = 12))
expect_error(dataSummarizationPTM(raw.input.tmt, reference_norm.PTM = 12))
expect_error(dataSummarizationPTM(raw.input.tmt, remove_norm_channel = 12))
expect_error(dataSummarizationPTM(raw.input.tmt, remove_empty_channel = ""))
expect_error(dataSummarizationPTM(raw.input.tmt, MBimpute = ""))
expect_error(dataSummarizationPTM(raw.input.tmt, MBimpute.PTM = ""))
expect_error(dataSummarizationPTM(raw.input.tmt, maxQuantileforCensored  = ""))
