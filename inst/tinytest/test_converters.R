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

mq_imported = MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                       annotation=maxq_tmt_annotation,
                       fasta=system.file("extdata", "maxq_tmt_fasta.fasta", package="MSstatsPTM"),
                       fasta_protein_name="uniprot_ac",
                       use_unmod_peptides=TRUE,
                       labeling_type = "TMT")
.validatePositiveNumberOfRows(mq_imported)
.validateProteinId(mq_imported$PROTEIN, "P29966", 10)
.validateProteinId(mq_imported$PTM, "P29966_T150", 70)
.validateProteinId(mq_imported$PTM, "P29966_T143_S145", 10)
.validatePtmSubstring(
    mq_imported$PTM, "Phospho \\(STY\\)", 
    length(mq_imported$PTM$PeptideSequence))
.validatePtmSubstring(
    mq_imported$PROTEIN, "Phospho \\(STY\\)", 0)

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

mq_imported = MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                                     annotation=maxq_lf_annotation,
                                     fasta=system.file("extdata", "maxq_lf_fasta.fasta", package="MSstatsPTM"),
                                     fasta_protein_name="uniprot_ac",
                                     mod_id="\\(Phospho \\(STY\\)\\)",
                                     use_unmod_peptides=TRUE,
                                     labeling_type = "LF",
                                     which_proteinid_ptm = "Proteins")
.validatePositiveNumberOfRows(mq_imported)
.validateProteinId(mq_imported$PROTEIN, "P36578", 144)
.validateProteinId(mq_imported$PTM, "P36578_S295", 72)
.validateProteinId(mq_imported$PTM, "Q13523_S431_S437", 36)
.validatePtmSubstring(
    mq_imported$PTM, "Phospho \\(STY\\)", 
    length(mq_imported$PTM$PeptideSequence))
.validatePtmSubstring(
    mq_imported$PROTEIN, "Phospho \\(STY\\)", 0)

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

spectronaut_imported = SpectronauttoMSstatsPTMFormat(spectronaut_input,
                                            annotation=spectronaut_annotation,
                                            fasta_path=system.file("extdata", "spectronaut_fasta.fasta", package="MSstatsPTM"),
                                            use_unmod_peptides=TRUE,
                                            mod_id = "\\[Phospho \\(STY\\)\\]",
                                            fasta_protein_name = "uniprot_iso"
)
.validatePositiveNumberOfRows(spectronaut_imported)
.validateProteinId(spectronaut_imported$PROTEIN, "P36578", 36)
.validateProteinId(spectronaut_imported$PTM, "P36578_S295", 36)
.validateProteinId(spectronaut_imported$PTM, "P09938_S15_S22", 36)
.validatePtmSubstring(
    spectronaut_imported$PTM, "Phospho \\(STY\\)", 
    length(spectronaut_imported$PTM$PeptideSequence))
.validatePtmSubstring(
    spectronaut_imported$PROTEIN, "Phospho \\(STY\\)", 0)

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

input = system.file("tinytest/raw_data/PD/pd-ptm-input.csv", 
                    package = "MSstatsPTM")
input = data.table::fread(input)
annot = system.file("tinytest/raw_data/PD/pd-ptm-annot.csv",
                    package = "MSstatsPTM")
annot = data.table::fread(annot)
input_protein = system.file("tinytest/raw_data/PD/protein-input.csv",
                            package = "MSstatsPTM")
input_protein = data.table::fread(input_protein)
annot_protein = system.file("tinytest/raw_data/PD/protein-annot.csv",
                            package = "MSstatsPTM")
annot_protein = data.table::fread(annot_protein)
fasta_path=system.file("extdata", "pd_with_proteome.fasta", 
                       package="MSstatsPTM")
pd_imported = PDtoMSstatsPTMFormat(
    input,
    annotation = annot,
    protein_input = input_protein,
    annotation_protein = annot_protein,
    fasta_path = fasta_path,
    mod_id = "\\(GG\\)",
    labeling_type = "TMT",
    use_localization_cutoff = FALSE,
    which_proteinid = "Master.Protein.Accessions")

.validatePositiveNumberOfRows(pd_imported)
.validateProteinId(pd_imported$PROTEIN, "P52480", 30)
.validateProteinId(pd_imported$PTM, "P52480_K305", 12)
.validateProteinId(pd_imported$PTM, "P52480", 0)
.validatePtmSubstring(
    pd_imported$PTM, "\\*",
    length(pd_imported$PTM$PeptideSequence))
.validatePtmSubstring(
    pd_imported$PROTEIN, "\\*", 0)

## Metamorpheus

#### Setup
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

#### Happy Case
expect_silent(MetamorpheusToMSstatsPTMFormat(
    input,
    annot,
    fasta_path=fasta_path,
    input_protein=input_protein,
    annotation_protein=annot_protein,
    mod_ids = c("\\[Common Fixed:Carbamidomethyl on C\\]")
))

#### Global profiling run removed
expect_silent(MetamorpheusToMSstatsPTMFormat(
    input,
    annot,
    fasta_path=fasta_path,
    use_unmod_peptides = FALSE,
    mod_ids = c("\\[Common Fixed:Carbamidomethyl on C\\]")
))

#### Global profiling run removed with use_unmod_peptides = TRUE
expect_silent(MetamorpheusToMSstatsPTMFormat(
    input,
    annot,
    fasta_path=fasta_path,
    use_unmod_peptides = TRUE,
    mod_ids = c("\\[Common Fixed:Carbamidomethyl on C\\]")
))

#### FASTA file missing
expect_error(MetamorpheusToMSstatsPTMFormat(
    input,
    annot,
    input_protein=input_protein,
    annotation_protein=annot_protein,
    mod_ids = c("\\[Common Fixed:Carbamidomethyl on C\\]")
))

#### Annotation missing
expect_error(MetamorpheusToMSstatsPTMFormat(
    input,
    fasta_path=fasta_path,
    input_protein=input_protein,
    annotation_protein=annot_protein,
    mod_ids = c("\\[Common Fixed:Carbamidomethyl on C\\]")
))

#### Output validation tests
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

expect_silent(ProgenesistoMSstatsPTMFormat(
    input,
    annot,
))

expect_error(ProgenesistoMSstatsPTMFormat(
    input
))

prog_imported = ProgenesistoMSstatsPTMFormat(
    input,
    annot,
)

.validateProteinId(prog_imported$PTM, 
    "sp|A2ASS6|TITIN_MOUSE_AVTSPPRVKSPEPR_[4] Phospho (ST)|[10] Phospho (ST)", 
    18
)


# Fragpipe
#### Happy Case
expect_silent(FragPipetoMSstatsPTMFormat(fragpipe_input,
                                         fragpipe_annotation,
                                         fragpipe_input_protein,
                                         fragpipe_annotation_protein,
                                         label_type="TMT",
                                         mod_id_col = "STY",
                                         localization_cutoff=.75,
                                         remove_unlocalized_peptides=TRUE)
              )

#### Global profiling run removed
expect_silent(FragPipetoMSstatsPTMFormat(fragpipe_input,
                                         fragpipe_annotation,
                                         label_type="TMT",
                                         mod_id_col = "STY",
                                         localization_cutoff=.75,
                                         remove_unlocalized_peptides=TRUE)
              )

#### use_unmod_peptides = TRUE, no unmodified proteins in dataset
expect_error(FragPipetoMSstatsPTMFormat(fragpipe_input,
                                         fragpipe_annotation,
                                         label_type="TMT",
                                         mod_id_col = "STY",
                                         localization_cutoff=.75,
                                         remove_unlocalized_peptides=TRUE,
                                         use_unmod_peptides=TRUE)
              )

#### Annotation missing
expect_error(FragPipetoMSstatsPTMFormat(fragpipe_input,
                                        fragpipe_input_protein,
                                        fragpipe_annotation_protein,
                                        label_type="TMT",
                                        mod_id_col = "STY",
                                        localization_cutoff=.75,
                                        remove_unlocalized_peptides=TRUE)
             )
fragpipe_input = FragPipetoMSstatsPTMFormat(fragpipe_input,
                                          fragpipe_annotation,
                                          fragpipe_input_protein,
                                          fragpipe_annotation_protein,
                                          label_type="TMT",
                                          mod_id_col = "STY",
                                          localization_cutoff=.75,
                                          remove_unlocalized_peptides=TRUE)
.validatePositiveNumberOfRows(fragpipe_input)
.validateProteinId(fragpipe_input$PROTEIN, "sp|Q9Y2W1|TR150_HUMAN", 110)
.validateProteinId(fragpipe_input$PTM, "sp|Q9Y2W1|TR150_HUMAN_S805", 20)
.validateProteinId(fragpipe_input$PTM, "sp|Q9Y2W1|TR150_HUMAN_S248_S253", 20)
.validatePtmSubstring(
    fragpipe_input$PTM, "\\*", 
    length(fragpipe_input$PTM$PeptideSequence))
.validatePtmSubstring(
    fragpipe_input$PROTEIN, "\\*", 0)

# LFQ Example (w/out global profiling run)
input = system.file("tinytest/raw_data/Fragpipe/MSstats.csv",
                                        package = "MSstatsPTM")
input = data.table::fread(input)
annot = system.file("tinytest/raw_data/Fragpipe/experiment_annotation.tsv",
                                        package = "MSstatsPTM")
annot = data.table::fread(annot)

#### Happy Case
expect_silent(FragPipetoMSstatsPTMFormat(input,
                                         annot,
                                         label_type="LF",
                                         mod_id_col = "STY",
                                         localization_cutoff=.75,
                                         protein_id_col = "ProteinName",
                                         peptide_id_col = "PeptideSequence")
)

#### Annotation file not provided but present in input
expect_silent(FragPipetoMSstatsPTMFormat(input,
                                        label_type="LF",
                                        mod_id_col = "STY",
                                        localization_cutoff=.75,
                                        protein_id_col = "ProteinName",
                                        peptide_id_col = "PeptideSequence")
)

fragpipe_input = FragPipetoMSstatsPTMFormat(input,
                                          annot,
                                          label_type="LF",
                                          mod_id_col = "STY",
                                          localization_cutoff=.75,
                                          protein_id_col = "ProteinName",
                                          peptide_id_col = "PeptideSequence")

.validateProteinId(fragpipe_input$PTM, "sp|P02400|RLA4_YEAST_S100", 8)
.validateProteinId(fragpipe_input$PTM, "sp|O13563|RPN13_YEAST_S132", 4)
.validatePtmSubstring(
    fragpipe_input$PTM, "\\*", 
    length(fragpipe_input$PTM$PeptideSequence))