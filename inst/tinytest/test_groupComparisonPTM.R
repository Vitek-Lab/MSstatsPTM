
data("summary.data", package = "MSstatsPTM")
data("summary.data.tmt", package = "MSstatsPTM")

## Test missing data is handled
expect_error(groupComparisonPTM())
expect_error(groupComparisonPTM(list(PTM = NULL, 
                                     PROTEIN = summary.data$PROTEIN), 
                                ptm_label_type="LF",
                                protein_label_type="LF"))

expect_silent(groupComparisonPTM(list(PTM = summary.data$PTM, 
                                     PROTEIN = NULL), 
                                 ptm_label_type="LF",
                                 protein_label_type="LF"))

## Test expected output - LF
group.comp <- MSstatsPTM::groupComparisonPTM(summary.data, 
                                             ptm_label_type="LF",
                                             protein_label_type="LF")
expect_inherits(group.comp, "list")
expect_inherits(group.comp$PTM.Model, "data.table")
expect_inherits(group.comp$PROTEIN.Model, "data.table")
expect_inherits(group.comp$ADJUSTED.Model, "data.table")

expect_equal(colnames(group.comp$PTM.Model), 
             c("Protein", "Label", "log2FC", "SE", "Tvalue", "DF", "pvalue",
               "adj.pvalue", "issue", "MissingPercentage", 
               "ImputationPercentage"))
expect_equal(colnames(group.comp$PROTEIN.Model), 
             c("Protein", "Label", "log2FC", "SE", "Tvalue", "DF", "pvalue",
               "adj.pvalue", "issue", "MissingPercentage", 
               "ImputationPercentage"))
expect_equal(colnames(group.comp$ADJUSTED.Model), 
             c("Protein", "Label", "log2FC", "SE", "Tvalue", "DF", "pvalue", 
               "adj.pvalue", "GlobalProtein", "Adjusted"))

## Test expected output - TMT
group.comp.tmt <- MSstatsPTM::groupComparisonPTM(summary.data.tmt, 
                                                 ptm_label_type="TMT",
                                                 protein_label_type="TMT")
expect_inherits(group.comp.tmt, "list")
expect_inherits(group.comp.tmt$PTM.Model, "data.table")
expect_inherits(group.comp.tmt$PROTEIN.Model, "data.table")
expect_inherits(group.comp.tmt$ADJUSTED.Model, "data.table")

expect_equal(colnames(group.comp.tmt$PTM.Model), 
             c("Protein", "Label", "log2FC", "SE", "DF", "pvalue", 
               "adj.pvalue", "issue"))
expect_equal(colnames(group.comp.tmt$PROTEIN.Model), 
             c("Protein", "Label", "log2FC", "SE", "DF", "pvalue", 
               "adj.pvalue", "issue"))
expect_equal(colnames(group.comp$ADJUSTED.Model), 
             c("Protein", "Label", "log2FC", "SE", "Tvalue", "DF", "pvalue", 
               "adj.pvalue", "GlobalProtein", "Adjusted"))
