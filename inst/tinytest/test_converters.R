
## MaxQ TMT
data("maxq_tmt_evidence", package = "MSstatsPTM")
data("maxq_tmt_annotation", package = "MSstatsPTM")

expect_silent(MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                                     annotation_ptm=maxq_tmt_annotation,
                                     fasta=system.file("extdata", "maxq_tmt_fasta.fasta", package="MSstatsPTM"),
                                     fasta_protein_name="uniprot_ac",
                                     use_unmod_peptides=TRUE,
                                     labeling_type = "TMT"))

expect_silent(MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                                     annotation_ptm=maxq_tmt_annotation,
                                     fasta=system.file("extdata", "maxq_tmt_fasta.fasta", package="MSstatsPTM"),
                                     fasta_protein_name="uniprot_ac",
                                     use_unmod_peptides=FALSE,
                                     labeling_type = "TMT"))

expect_error(MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                                     annotation_ptm=maxq_tmt_annotation,
                                     fasta_protein_name="uniprot_ac",
                                     use_unmod_peptides=FALSE,
                                     labeling_type = "TMT"))

expect_error(MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                                    fasta_protein_name="uniprot_ac",
                                    fasta=system.file("extdata", "maxq_tmt_fasta.fasta", package="MSstatsPTM"),
                                    use_unmod_peptides=FALSE,
                                    labeling_type = "TMT"))

expect_silent(MaxQtoMSstatsPTMFormat(evidence=maxq_tmt_evidence,
                                     annotation_ptm=maxq_tmt_annotation,
                                     fasta=system.file("extdata", "maxq_tmt_fasta.fasta", package="MSstatsPTM"),
                                     fasta_protein_name="uniprot_ac",
                                     use_unmod_peptides=FALSE,
                                     labeling_type = "TMT",
                                     removeOxidationMpeptides=TRUE))

## MaxQ LF
data("maxq_lf_evidence", package = "MSstatsPTM")
data("maxq_lf_annotation", package = "MSstatsPTM")

expect_silent(MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                                    annotation_ptm=maxq_lf_annotation,
                                    fasta=system.file("extdata", "maxq_lf_fasta.fasta", package="MSstatsPTM"),
                                    fasta_protein_name="uniprot_ac",
                                    mod_id="\\(Phospho \\(STY\\)\\)",
                                    use_unmod_peptides=TRUE,
                                    labeling_type = "LF",
                                    which_proteinid_ptm = "Proteins"))

expect_error(MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                                     annotation_ptm=maxq_lf_annotation,
                                     fasta=system.file("extdata", "maxq_lf_fasta.fasta", package="MSstatsPTM"),
                                     fasta_protein_name="uniprot_ac",
                                     mod_id="\\(Phospho \\(STY\\)\\)",
                                     use_unmod_peptides=TRUE,
                                     labeling_type = "TMT",
                                     which_proteinid_ptm = "Proteins"))

expect_error(MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                                    annotation_ptm=maxq_lf_annotation,
                                    fasta_protein_name="uniprot_ac",
                                    mod_id="\\(Phospho \\(STY\\)\\)",
                                    use_unmod_peptides=TRUE,
                                    labeling_type = "TMT",
                                    which_proteinid_ptm = "Proteins"))

expect_error(MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                                    annotation_ptm=maxq_lf_annotation,
                                    fasta=system.file("extdata", "maxq_lf_fasta.fasta", package="MSstatsPTM"),
                                    fasta_protein_name="uniprot_ac",
                                    mod_id="\\(Phospho \\(STY\\)\\)",
                                    use_unmod_peptides=TRUE,
                                    labeling_type = "TMT",
                                    which_proteinid_ptm = "Leading.proteins"))

expect_silent(MaxQtoMSstatsPTMFormat(evidence=maxq_lf_evidence,
                                     annotation_ptm=maxq_lf_annotation,
                                     fasta=system.file("extdata", "maxq_lf_fasta.fasta", package="MSstatsPTM"),
                                     fasta_protein_name="uniprot_ac",
                                     mod_id="\\(Phospho \\(STY\\)\\)",
                                     use_unmod_peptides=FALSE,
                                     labeling_type = "LF",
                                     which_proteinid_ptm = "Proteins"))