## Utility functions

#' Validate positive number of rows from converter output
#' 
#' @author Anthony Wu
#' 
#' @param msstats_ptm_input A list containing PTM and PROTEIN data tables
.validatePositiveNumberOfRows = function(msstats_ptm_input) {
    expect_true(nrow(msstats_ptm_input$PTM) > 0)
    expect_true(nrow(msstats_ptm_input$PROTEIN) > 0)
}

#' Validate protein ID count in converter output
#' 
#' @author Anthony Wu
#' 
#' @param msstats_input Either PTM or PROTEIN data tables
#' @param protein_id Protein ID to validate
#' @param count Expected count of protein ID in ProteinName column
.validateProteinId = function(msstats_input, protein_id, count) {
    expect_equal(length(which(msstats_input$ProteinName==protein_id)), count)
}

#' Validate ptm ID count in converter output
#' 
#' @author Anthony Wu
#' 
#' @param msstats_input Either PTM or PROTEIN data tables
#' @param ptm_substring PTM to validate
#' @param count Expected count of PTM in PeptideSequence column
.validatePtmSubstring = function(msstats_input, ptm_substring, count) {
    expect_equal(
        length(which(grepl(ptm_substring, msstats_input$PeptideSequence))),
        count
    )
}

## MaxQ TMT
data("maxq_tmt_evidence", package = "MSstatsPTM")
data("maxq_tmt_annotation", package = "MSstatsPTM")

expect_silent(MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                                     annotation=maxq_tmt_annotation,
                                     fasta=system.file("extdata", "maxq_tmt_fasta.fasta", package="MSstatsPTM"),
                                     fasta_protein_name="uniprot_ac",
                                     use_unmod_peptides=TRUE,
                                     labeling_type = "TMT"))

expect_silent(MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                                     annotation=maxq_tmt_annotation,
                                     fasta=system.file("extdata", "maxq_tmt_fasta.fasta", package="MSstatsPTM"),
                                     fasta_protein_name="uniprot_ac",
                                     use_unmod_peptides=FALSE,
                                     labeling_type = "TMT"))

expect_error(MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                                    annotation=maxq_tmt_annotation,
                                     fasta_protein_name="uniprot_ac",
                                     use_unmod_peptides=FALSE,
                                     labeling_type = "TMT"))

expect_error(MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                                    fasta_protein_name="uniprot_ac",
                                    fasta=system.file("extdata", "maxq_tmt_fasta.fasta", package="MSstatsPTM"),
                                    use_unmod_peptides=FALSE,
                                    labeling_type = "TMT"))

expect_silent(MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                                     annotation=maxq_tmt_annotation,
                                     fasta=system.file("extdata", "maxq_tmt_fasta.fasta", package="MSstatsPTM"),
                                     fasta_protein_name="uniprot_ac",
                                     use_unmod_peptides=FALSE,
                                     labeling_type = "TMT",
                                     removeOxidationMpeptides=TRUE))

## MaxQ LF
data("maxq_lf_evidence", package = "MSstatsPTM")
data("maxq_lf_annotation", package = "MSstatsPTM")

expect_silent(MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                                     annotation=maxq_lf_annotation,
                                    fasta=system.file("extdata", "maxq_lf_fasta.fasta", package="MSstatsPTM"),
                                    fasta_protein_name="uniprot_ac",
                                    mod_id="\\(Phospho \\(STY\\)\\)",
                                    use_unmod_peptides=TRUE,
                                    labeling_type = "LF",
                                    which_proteinid_ptm = "Proteins"))

expect_error(MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                                    annotation=maxq_lf_annotation,
                                     fasta=system.file("extdata", "maxq_lf_fasta.fasta", package="MSstatsPTM"),
                                     fasta_protein_name="uniprot_ac",
                                     mod_id="\\(Phospho \\(STY\\)\\)",
                                     use_unmod_peptides=TRUE,
                                     labeling_type = "TMT",
                                     which_proteinid_ptm = "Proteins"))

expect_error(MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                                    annotation=maxq_lf_annotation,
                                    fasta_protein_name="uniprot_ac",
                                    mod_id="\\(Phospho \\(STY\\)\\)",
                                    use_unmod_peptides=TRUE,
                                    labeling_type = "TMT",
                                    which_proteinid_ptm = "Proteins"))

expect_error(MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                                    annotation=maxq_lf_annotation,
                                    fasta=system.file("extdata", "maxq_lf_fasta.fasta", package="MSstatsPTM"),
                                    fasta_protein_name="uniprot_ac",
                                    mod_id="\\(Phospho \\(STY\\)\\)",
                                    use_unmod_peptides=TRUE,
                                    labeling_type = "TMT",
                                    which_proteinid_ptm = "Leading.proteins"))

expect_silent(MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                                     annotation=maxq_lf_annotation,
                                     fasta=system.file("extdata", "maxq_lf_fasta.fasta", package="MSstatsPTM"),
                                     fasta_protein_name="uniprot_ac",
                                     mod_id="\\(Phospho \\(STY\\)\\)",
                                     use_unmod_peptides=FALSE,
                                     labeling_type = "LF",
                                     which_proteinid_ptm = "Proteins"))

## Spectronaut
data("spectronaut_input", package = "MSstatsPTM")
data("spectronaut_annotation", package = "MSstatsPTM")
expect_silent(SpectronauttoMSstatsPTMFormat(spectronaut_input,
                  annotation=spectronaut_annotation,
                  fasta_path=system.file("extdata", "spectronaut_fasta.fasta", package="MSstatsPTM"),
                  use_unmod_peptides=TRUE,
                  mod_id = "\\[Phospho \\(STY\\)\\]",
                  fasta_protein_name = "uniprot_iso"
                  ))

expect_error(SpectronauttoMSstatsPTMFormat(spectronaut_input,
                                            annotation=spectronaut_annotation,
                                            use_unmod_peptides=TRUE,
                                            mod_id = "[[[Phospho \\(STY\\)\\]",
                                            fasta_protein_name = "uniprot_iso"
))

expect_silent(SpectronauttoMSstatsPTMFormat(spectronaut_input,
                                           fasta_path=system.file("extdata", "spectronaut_fasta.fasta", package="MSstatsPTM"),
                                           use_unmod_peptides=TRUE,
                                           mod_id = "\\[Phospho \\(STY\\)\\]",
                                           fasta_protein_name = "uniprot_iso"
))

expect_silent(SpectronauttoMSstatsPTMFormat(spectronaut_input,
                                            annotation=spectronaut_annotation,
                                            fasta_path=system.file("extdata", "spectronaut_fasta.fasta", package="MSstatsPTM"),
                                            mod_id = "\\[Phospho \\(STY\\)\\]",
                                            fasta_protein_name = "uniprot_iso"
))

## PD
data("pd_psm_input", package = "MSstatsPTM")
data("pd_annotation", package = "MSstatsPTM")
data("pd_testing_output", package = "MSstatsPTM")

expect_message(PDtoMSstatsPTMFormat(pd_psm_input, 
    pd_annotation,
    system.file("extdata", "pd_fasta.fasta", package="MSstatsPTM"),
    use_unmod_peptides=TRUE,
    which_proteinid = "Master.Protein.Accessions"))

expect_equal(pd_testing_output, PDtoMSstatsPTMFormat(pd_psm_input, 
                                   pd_annotation,
                                   system.file("extdata", "pd_fasta.fasta", 
                                               package="MSstatsPTM"),
                                   use_unmod_peptides=TRUE,
                                   which_proteinid = "Master.Protein.Accessions"))
 
expect_equal(pd_testing_output, PDtoMSstatsPTMFormat(pd_psm_input, 
                                    pd_annotation,
                                    system.file("extdata", "pd_fasta.fasta", 
                                                package="MSstatsPTM"),
                                    use_unmod_peptides=TRUE,
                                    which_proteinid = "Master.Protein.Accessions"))

expect_message(PDtoMSstatsPTMFormat(pd_psm_input, 
                                   pd_annotation,
                                   system.file("extdata", "pd_fasta.fasta", 
                                               package="MSstatsPTM"),
                                   use_unmod_peptides=TRUE,
                                   which_proteinid = "Master.Protein.Accessions",
                                   use_localization_cutoff = FALSE))

expect_error(PDtoMSstatsPTMFormat(pd_psm_input, 
                                   pd_annotation,
                                   use_unmod_peptides=TRUE,
                                   which_proteinid = "Master.Protein.Accessions"))

expect_message(PDtoMSstatsPTMFormat(pd_psm_input, 
                                  pd_annotation,
                                  system.file("extdata", "pd_fasta.fasta", 
                                              package="MSstatsPTM"),
                                  which_proteinid = "Master.Protein.Accessions"))

expect_error(PDtoMSstatsPTMFormat(pd_psm_input, 
                                  NULL,
                                  system.file("extdata", "pd_fasta.fasta", 
                                              package="MSstatsPTM"),
                                  which_proteinid = "Master.Protein.Accessions"))

expect_error(PDtoMSstatsPTMFormat(pd_psm_input, 
                                  data.frame(columns=c("Condition", "BioReplicate")),
                                  system.file("extdata", "pd_fasta.fasta", 
                                              package="MSstatsPTM"),
                                  which_proteinid = "Master.Protein.Accessions"))


## Metamorpheus
input = system.file("tinytest/raw_data/Metamorpheus/AllQuantifiedPeaks.tsv", 
                                package = "MSstatsPTM")
input = data.table::fread(input)
annot = system.file("tinytest/raw_data/Metamorpheus/ExperimentalDesign.tsv",
                                package = "MSstatsPTM")
annot = data.table::fread(annot)
input_protein = system.file("tinytest/raw_data/Metamorpheus/AllQuantifiedPeaksGlobalProteome.tsv",
                                package = "MSstatsPTM")
input_protein = data.table::fread(input_protein)
annot_protein = system.file("tinytest/raw_data/Metamorpheus/ExperimentalDesignGlobalProteome.tsv",
                                package = "MSstatsPTM")
annot_protein = data.table::fread(annot_protein)
fasta_path=system.file("extdata", "metamorpheus_fasta.fasta", 
                       package="MSstatsPTM")
metamorpheus_imported = MetamorpheusToMSstatsPTMFormat(
        input,
        annot,
        fasta_path=fasta_path,
        input_protein=input_protein,
        annotation_protein=annot_protein,
        mod_ids = c("\\[Common Fixed:Carbamidomethyl on C\\]")
    )
.validatePositiveNumberOfRows(metamorpheus_imported)
.validateProteinId(metamorpheus_imported$PROTEIN, "O95817", 12)
.validateProteinId(metamorpheus_imported$PTM, "O95817_C207", 6)
.validateProteinId(metamorpheus_imported$PTM, "O95817", 0)
.validatePtmSubstring(
    metamorpheus_imported$PTM, "\\[Common Fixed:Carbamidomethyl on C\\]", 
    length(metamorpheus_imported$PTM$PeptideSequence))
.validatePtmSubstring(
    metamorpheus_imported$PROTEIN, "\\[Common Fixed:Carbamidomethyl on C\\]", 0)

metamorpheus_imported = MetamorpheusToMSstatsPTMFormat(
    input,
    annot,
    fasta_path=fasta_path,
    use_unmod_peptides = TRUE,
    mod_ids = c("\\[Common Fixed:Carbamidomethyl on C\\]")
)
.validatePositiveNumberOfRows(metamorpheus_imported)
.validateProteinId(metamorpheus_imported$PROTEIN, "O95817", 6)
.validateProteinId(metamorpheus_imported$PTM, "O95817_C207", 6)
.validateProteinId(metamorpheus_imported$PTM, "O95817", 0)
.validatePtmSubstring(
    metamorpheus_imported$PTM, "\\[Common Fixed:Carbamidomethyl on C\\]", 
    length(metamorpheus_imported$PTM$PeptideSequence))
.validatePtmSubstring(
    metamorpheus_imported$PROTEIN, "\\[Common Fixed:Carbamidomethyl on C\\]", 0)


# Progenesis
input = system.file("tinytest/raw_data/Progenesis/progenesis_peptide.csv", 
                    package = "MSstatsPTM")
input = data.table::fread(input)
colnames(input) = unlist(input[1,])
input = input[-1,]
annot = system.file("tinytest/raw_data/Progenesis/phospho_annotation.csv",
                    package = "MSstatsPTM")
annot = data.table::fread(annot)
prog_imported = ProgenesistoMSstatsPTMFormat(
    input,
    annot,
)

.validateProteinId(prog_imported$PTM, "sp|A2ASS6|TITIN_MOUSE_AVTSPPRVKSPEPR_[4] Phospho (ST)|[10] Phospho (ST)", 18)