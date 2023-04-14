
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

expect_silent(PDtoMSstatsPTMFormat(pd_psm_input, 
    pd_annotation,
    system.file("extdata", "pd_fasta.fasta", package="MSstatsPTM"),
    use_unmod_peptides=TRUE,
    which_proteinid = "Master.Protein.Accessions"))

expect_silent(PDtoMSstatsPTMFormat(pd_psm_input, 
                                   pd_annotation,
                                   system.file("extdata", "pd_fasta.fasta", package="MSstatsPTM"),
                                   use_unmod_peptides=TRUE,
                                   which_proteinid = "Master.Protein.Accessions",
                                   use_localization_cutoff = FALSE))

expect_error(PDtoMSstatsPTMFormat(pd_psm_input, 
                                   pd_annotation,
                                   use_unmod_peptides=TRUE,
                                   which_proteinid = "Master.Protein.Accessions"))

expect_silent(PDtoMSstatsPTMFormat(pd_psm_input, 
                                  pd_annotation,
                                  system.file("extdata", "pd_fasta.fasta", package="MSstatsPTM"),
                                  which_proteinid = "Master.Protein.Accessions"))

expect_error(PDtoMSstatsPTMFormat(pd_psm_input, 
                                  NULL,
                                  system.file("extdata", "pd_fasta.fasta", package="MSstatsPTM"),
                                  which_proteinid = "Master.Protein.Accessions"))

expect_error(PDtoMSstatsPTMFormat(pd_psm_input, 
                                  data.frame(columns=c("Condition", "BioReplicate")),
                                  system.file("extdata", "pd_fasta.fasta", package="MSstatsPTM"),
                                  which_proteinid = "Master.Protein.Accessions"))
