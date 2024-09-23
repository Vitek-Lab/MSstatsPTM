test_that("Medians are equalized in ptm data to QC plot helper function", {
    mock_qc_all_plot_lf <- mock()
    stub(
        dataProcessPlotsPTM, ".qc.all.plot.lf",
        mock_qc_all_plot_lf, depth=2
    )
    
    dataProcessPlotsPTM(summary.data,
                        type = 'QCPLOT',
                        which.PTM = "allonly",
                        address = FALSE)
    
    calls <- mock_args(mock_qc_all_plot_lf)
    ptm_data <- calls[[1]][[1]]
    run_data = split(ptm_data, ptm_data$RUN)
    for (i in 1:length(run_data)){
        expect_equal(median(run_data[[i]]$ABUNDANCE, na.rm = TRUE), 22.028431)
    }
})
