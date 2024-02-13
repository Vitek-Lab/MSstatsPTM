#' Example annotation file for a label-free MaxQuant experiment.
#' 
#' Must be manually created by the user and input into the 
#' MaxQtoMSstatsPTMFormat converter. Requires the correct columns and maps the 
#' experimental desing into the MSstats format. Specify unique bioreplicates for
#' group comparison designs, and the same bioreplicate for repeated measure 
#' designs. The columns and descriptions are below.
#' 
#' \itemize{
#'   \item Run : Run name that matches exactly with MaxQuant run. Used to join 
#'   evidence and metadata in annotation file.
#'   \item Condition : Name of condition that was used for each run.
#'   \item BioReplicate : Name of biological replicate. Repeating the same name 
#'   here will tell MSstatsPTM that the experiment is a repeated measure design.
#'   \item Raw.file : Run name that matches exactly with MaxQuant run. Used to join 
#'   evidence and metadata in annotation file.
#'   \item IsotopeLabelType: Name of isotope label. May be all `L` or unique 
#'   depending on experimental design.
#' }
#' 
#' @format A data.table with 5 columns.
#' @examples
#' head(maxq_lf_annotation)
"maxq_lf_annotation"

#' Example MaxQuant evidence file from the output of a label free experiment
#' 
#' Experiment was performed by the Olsen lab and published on Nat. Commun.
#' (citation below).
#' 
#' Bekker-Jensen, D.B., Bernhardt, O.M., Hogrebe, A. et al. Rapid and 
#' site-specific deep phosphoproteome profiling by data-independent acquisition 
#' without the need for spectral libraries. Nat Commun 11, 787 (2020). 
#' https://doi.org/10.1038/s41467-020-14609-1
#' 
#' The experiment was processed using MaxQuant by the computational proteomics 
#' team at Pfizer (Liang Xue and Pierre Jean). 
#' 
#' The experiment did not contain a global profiling run, but we show an example 
#' of extracting the unmodified peptides and using them in place of the 
#' profiling run.
#' 
#' @format a data.table with 63 columns and 511 rows, the output of MaxQuant
#' @examples
#' head(maxq_lf_evidence)
"maxq_lf_evidence"

#' Example annotation file for a TMT MaxQuant experiment.
#' 
#' Must be manually created by the user and input into the 
#' MaxQtoMSstatsPTMFormat converter. Requires the correct columns and maps the 
#' experimental desing into the MSstats format. Specify unique bioreplicates for
#' group comparison designs, and the same bioreplicate for repeated measure 
#' designs. The columns and descriptions are below.
#' 
#' \itemize{
#'   \item Run : Run name that matches exactly with MaxQuant run. Used to join 
#'   evidence and metadata in annotation file.
#'   \item Fraction : If multiple fractions were used (i.e. the same mixture 
#'   split into multiple fractions) enter that here.
#'   TechRepMixture : Multiple runs using the same bioreplicate
#'   \item Channel : Mixture channel used
#'   \item Condition : Name of condition that was used for each run.
#'   \item Mixture : The unique mixture (plex) name
#'   \item BioReplicate : Name of biological replicate. Repeating the same name 
#'   here will tell MSstatsPTM that the experiment is a repeated measure design.
#' }
#' 
#' @format A data.table with 7 columns.
#' @examples
#' head(maxq_tmt_annotation)
"maxq_tmt_annotation"

#' Example MaxQuant evidence file from the output of a TMT experiment
#' 
#' Experiment was performed by the Olsen lab and published on Nat. Commun.
#' (citation below).
#' 
#' Hogrebe, A., von Stechow, L., Bekker-Jensen, D.B. et al. Benchmarking common 
#' quantification strategies for large-scale phosphoproteomics. Nat Commun 9, 
#' 1045 (2018). https://doi.org/10.1038/s41467-018-03309-6
#' 
#' The experiment was processed using MaxQuant by the computational proteomics 
#' team at Pfizer (Liang Xue and Pierre Jean). 
#' 
#' The experiment did not contain a global profiling run, but we show an example 
#' of extracting the unmodified peptides and using them in place of the 
#' profiling run.
#' 
#' @format a data.table with 96 columns and 199 rows, the output of MaxQuant
#' @examples
#' head(maxq_tmt_evidence)
"maxq_tmt_evidence"

#' Example annotation file for a label-free Proteome Discoverer experiment.
#' 
#' Must be manually created by the user and input into the 
#' PDtoMSstatsPTMFormat converter. Requires the correct columns and maps the 
#' experimental desing into the MSstats format. Specify unique bioreplicates for
#' group comparison designs, and the same bioreplicate for repeated measure 
#' designs. The columns and descriptions are below.
#' 
#' \itemize{
#'   \item Run : Run name that matches exactly with PD run. Used to join 
#'   evidence and metadata in annotation file.
#'   \item Condition : Name of condition that was used for each run.
#'   \item BioReplicate : Name of biological replicate. Repeating the same name 
#'   here will tell MSstatsPTM that the experiment is a repeated measure design.
#' }
#' 
#' @format A data.table with 3 columns.
#' @examples
#' head(pd_annotation)
"pd_annotation"

#' Example Proteome Discoverer evidence file from the output of a label free experiment
#' 
#' Experiment was performed by the Olsen lab and published on Nat. Commun.
#' (citation below).
#' 
#' Bekker-Jensen, D.B., Bernhardt, O.M., Hogrebe, A. et al. Rapid and 
#' site-specific deep phosphoproteome profiling by data-independent acquisition 
#' without the need for spectral libraries. Nat Commun 11, 787 (2020). 
#' https://doi.org/10.1038/s41467-020-14609-1
#' 
#' The experiment was processed using Proteome Discoverer by the computational 
#' proteomics team at Pfizer (Liang Xue and Pierre Jean). 
#' 
#' The experiment did not contain a global profiling run, but we show an example 
#' of extracting the unmodified peptides and using them in place of the 
#' profiling run.
#' 
#' @format a data.table with 60 columns and 1657 rows, the output of PD
#' @examples
#' head(pd_psm_input)
"pd_psm_input"

#' Example output of Proteome Discoverer converter
#' 
#' output using example data provided in package
#' 
#' The experiment did not contain a global profiling run, but we show an example 
#' of extracting the unmodified peptides and using them in place of the 
#' profiling run.
#' 
#' @format a list with 2 data.frames
#' @examples
#' head(pd_testing_output)
"pd_testing_output"

#' Example annotation file for a label-free Spectronaut experiment.
#' 
#' Must be manually created by the user and input into the 
#' SpectronauttoMSstatsPTMFormat converter. Requires the correct columns and 
#' maps the experimental desing into the MSstats format. Specify unique 
#' bioreplicates for group comparison designs, and the same bioreplicate for 
#' repeated measure designs. The columns and descriptions are below.
#' 
#' \itemize{
#'   \item Run : Run name that matches exactly with Spectronaut run. Used to join 
#'   evidence and metadata in annotation file.
#'   \item Condition : Name of condition that was used for each run.
#'   \item BioReplicate : Name of biological replicate. Repeating the same name 
#'   here will tell MSstatsPTM that the experiment is a repeated measure design.
#'   \item Raw.file : Run name that matches exactly with Spectronaut run. Used to join 
#'   evidence and metadata in annotation file.
#' }
#' 
#' @format A data.table with 5 columns.
#' @examples
#' head(spectronaut_annotation)
"spectronaut_annotation"

#' Example Spectronaut evidence file from the output of a label free experiment
#' 
#' Experiment was performed by the Olsen lab and published on Nat. Commun.
#' (citation below).
#' 
#' Bekker-Jensen, D.B., Bernhardt, O.M., Hogrebe, A. et al. Rapid and 
#' site-specific deep phosphoproteome profiling by data-independent acquisition 
#' without the need for spectral libraries. Nat Commun 11, 787 (2020). 
#' https://doi.org/10.1038/s41467-020-14609-1
#' 
#' The experiment was processed using Spectronaut by the computational proteomics 
#' team at Pfizer (Liang Xue and Pierre Jean). 
#' 
#' The experiment did not contain a global profiling run, but we show an example 
#' of extracting the unmodified peptides and using them in place of the 
#' profiling run.
#' 
#' @format a data.table with 23 columns and 2683 rows, the output of Spectronaut
#' @examples
#' head(spectronaut_input)
"spectronaut_input"

#' Example annotation file for a TMT FragPipe experiment.
#' 
#' Automatically created by FragPipe, manually checked by the user and input 
#' into the FragPipetoMSstatsPTMFormat converter. Requires the correct columns 
#' and maps the experimental desing into the MSstats format. Specify unique 
#' bioreplicates for group comparison designs, and the same bioreplicate for 
#' repeated measure designs. The columns and descriptions are below.
#' 
#' \itemize{
#'   \item Run : Run name that matches exactly with FragPipe run. Used to join 
#'   evidence and metadata in annotation file.
#'   \item Fraction : If multiple fractions were used (i.e. the same mixture 
#'   split into multiple fractions) enter that here.
#'   TechRepMixture : Multiple runs using the same bioreplicate
#'   \item Channel : Mixture channel used
#'   \item Condition : Name of condition that was used for each run.
#'   \item Mixture : The unique mixture (plex) name
#'   \item BioReplicate : Name of biological replicate. Repeating the same name 
#'   here will tell MSstatsPTM that the experiment is a repeated measure design.
#' }
#' 
#' @format A data.table with 7 columns.
#' @examples
#' head(fragpipe_annotation)
"fragpipe_annotation"

#' Example annotation file for a global profiling run TMT FragPipe experiment.
#' 
#' Automatically created by FragPipe, manually checked by the user and input 
#' into the FragPipetoMSstatsPTMFormat converter. Requires the correct columns 
#' and maps the experimental desing into the MSstats format. Specify unique 
#' bioreplicates for group comparison designs, and the same bioreplicate for 
#' repeated measure designs. The columns and descriptions are below.
#' 
#' \itemize{
#'   \item Run : Run name that matches exactly with FragPipe run. Used to join 
#'   evidence and metadata in annotation file.
#'   \item Fraction : If multiple fractions were used (i.e. the same mixture 
#'   split into multiple fractions) enter that here.
#'   TechRepMixture : Multiple runs using the same bioreplicate
#'   \item Channel : Mixture channel used
#'   \item Condition : Name of condition that was used for each run.
#'   \item Mixture : The unique mixture (plex) name
#'   \item BioReplicate : Name of biological replicate. Repeating the same name 
#'   here will tell MSstatsPTM that the experiment is a repeated measure design.
#' }
#' 
#' @format A data.table with 7 columns.
#' @examples
#' head(fragpipe_annotation_protein)
"fragpipe_annotation_protein"

#' Output of FragPipe TMT PTM experiment
#' 
#' This dataset was provided by the FragPipe team at the Nesvilab. It was 
#' processed using Philosopher and targeted Phosphorylation.
#' 
#' @format A data.table with 29 columns and 246 rows.
#' @examples
#' head(fragpipe_input)
"fragpipe_input"

#' Output of FragPipe TMT global profiling experiment
#' 
#' This dataset was provided by the FragPipe team at the Nesvilab. It was 
#' processed using Philosopher and targeted Phosphorylation.
#' 
#' @format A data.table with 27 columns and 47 rows.
#' @examples
#' head(fragpipe_input_protein)
"fragpipe_input_protein"

#' Example of input PTM dataset for TMT experiments.
#' 
#' It can be the output of MSstatsPTM converter MaxQtoMSstatsPTMFormat or other 
#' MSstatsTMT converter functions (Please see MSstatsPTM_TMT_Workflow vignette).
#' The dataset is formatted as a list with two data.tables named PTM and 
#' PROTEIN. In each data.table the variables are as follows:
#' 
#' \itemize{
#'   \item ProteinName : Name of protein with modification site mapped in with
#'    an underscore. ie "Protein_4_Y474"
#'   \item PeptideSequence
#'   \item Charge
#'   \item PSM
#'   \item Mixture : Mixture of samples labeled with different TMT reagents,
#'    which can be analyzed in
#'   a single mass spectrometry experiment. If the channal doesn't have sample,
#'    please add `Empty' under Condition.
#'   \item TechRepMixture : Technical replicate of one mixture. One mixture may
#'    have multiple technical replicates.
#'   For example, if `TechRepMixture' = 1, 2 are the two technical replicates of
#'    one mixture, then they should match
#'   with same `Mixture' value.
#'   \item Run : MS run ID.
#'   \item Channel : Labeling information (126, ... 131).
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item BioReplicate : Unique ID for biological subject. If the channal
#'   doesn't have sample, please add `Empty' under BioReplicate.
#'   \item Intensity
#' }
#'
#' @format A list of two data.tables named PTM and PROTEIN with 1716 and 29221 
#' rows respectively.
#' @examples
#' head(raw.input.tmt$PTM)
#' head(raw.input.tmt$PROTEIN)
#'
"raw.input.tmt"

#' Example of input PTM dataset for LabelFree/DDA/DIA experiments.
#' 
#' It can be the output of MSstatsPTM converter ProgenesistoMSstatsPTMFormat or 
#' other MSstats converter functions (Please see MSstatsPTM_LabelFree_Workflow 
#' vignette). The dataset is formatted as a list with two data.tables named PTM 
#' and PROTEIN. In each data.table the variables are as follows:
#' 
#' \itemize{
#'   #'   \item ProteinName : Name of protein with modification site mapped in with
#'    an underscore. ie "Protein_4_Y474"
#'   \item PeptideSequence
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item BioReplicate : Unique ID for biological subject.
#'   \item Run : MS run ID.
#'   \item Intensity
#'   \item PrecursorCharge 
#'   \item FragmentIon
#'   \item ProductCharge
#'   \item IsotopeLabelType
#' }
#'
#' @format A list of two data.tables named PTM and PROTEIN with 1745 and 478 
#' rows respectively.
#' @examples
#' head(raw.input$PTM)
#' head(raw.input$PROTEIN)
#'
"raw.input"

#' Example of output from dataSummarizationPTM function for non-TMT data
#'
#' It is made from \code{\link{raw.input}}.
#' It is the output of dataSummarizationPTM function from MSstatsPTM.
#' It should include a list with two names \code{PTM} and \code{PROTEIN}. Each 
#' of these list values is also a list with two names \code{ProteinLevelData} 
#' and \code{FeatureLevelData}, which correspond to two data.tables.The columns 
#' in these two data.tables are listed below. The variables are as follows:
#' \itemize{
#'   \item FeatureLevelData : \itemize{
#'     \item PROTEIN : Protein ID with modification site mapped in. Ex.
#'     Protein_1002_S836
#'     \item PEPTIDE : Full peptide with charge
#'     \item TRANSITION: Charge
#'     \item FEATURE : Combination of Protien, Peptide, and Transition Columns
#'     \item LABEL : 
#'     \item GROUP : Condition (ex. Healthy, Cancer, Time0)
#'     \item RUN : Unique ID for technical replicate of one TMT
#'     mixture.
#'     \item SUBJECT : Unique ID for biological subject.
#'     \item FRACTION : Unique Fraction ID
#'     \item originalRUN : Run name
#'     \item censored : 
#'     \item INTENSITY : Unique ID for TMT mixture.
#'     \item ABUNDANCE : Unique ID for TMT mixture.
#'     \item newABUNDANCE : Unique ID for TMT mixture.
#'     \item predicted : Unique ID for TMT mixture.
#'   }
#'   \item ProteinLevelData : \itemize{
#'     \item RUN : MS run ID
#'     \item Protein : Protein ID with modification site mapped in. Ex.
#'     Protein_1002_S836
#'     \item LogIntensities: Protein-level summarized abundance
#'     \item originalRUN : Labeling information (126, ... 131)
#'     \item GROUP : Condition (ex. Healthy, Cancer, Time0)
#'     \item SUBJECT : Unique ID for biological subject.
#'     \item TotalGroupMeasurements : Unique ID for technical replicate of one TMT
#'     mixture.
#'     \item NumMeasuredFeature : Unique ID for TMT mixture.
#'     \item MissingPercentage : Unique ID for TMT mixture.
#'     \item more50missing : Unique ID for TMT mixture.
#'     \item NumImputedFeature : Unique ID for TMT mixture.
#'   }
#'    
#'   }
#'
#' @format A list of two lists with four data.tables.
#' @examples
#' head(summary.data)
#'
"summary.data"

#' Example of output from dataSummarizationPTM_TMT function for TMT data
#'
#' It is made from \code{\link{raw.input.tmt}}.
#' It is the output of dataSummarizationPTM_TMT function from MSstatsPTM.
#' It should include a list with two names \code{PTM} and \code{PROTEIN}. Each 
#' of these list values is also a list with two names \code{ProteinLevelData} 
#' and \code{FeatureLevelData}, which correspond to two data.tables.The columns 
#' in these two data.tables are listed below. The variables are as follows:
#' \itemize{
#'   \item FeatureLevelData : \itemize{
#'     \item ProteinName : MS run ID
#'     \item PSM : Protein ID with modification site mapped in. Ex.
#'     Protein_1002_S836
#'     \item censored: Protein-level summarized abundance
#'     \item predicted : Labeling information (126, ... 131)
#'     \item log2Intensity : Condition (ex. Healthy, Cancer, Time0)
#'     \item Run : Unique ID for biological subject.
#'     \item Channel : Unique ID for technical replicate of one TMT
#'     mixture.
#'     \item BioReplicate : Unique ID for TMT mixture.
#'     \item Condition : Unique ID for TMT mixture.
#'     \item Mixture : Unique ID for TMT mixture.
#'     \item TechRepMixture : Unique ID for TMT mixture.
#'     \item PeptideSequence : Unique ID for TMT mixture.
#'     \item Charge : Unique ID for TMT mixture.
#'   }
#'   \item ProteinLevelData : \itemize{
#'     \item Mixture : MS run ID
#'     \item TechRepMixture : Protein ID with modification site mapped in. Ex.
#'     Protein_1002_S836
#'     \item Run: Protein-level summarized abundance
#'     \item Channel : Labeling information (126, ... 131)
#'     \item Protein : Condition (ex. Healthy, Cancer, Time0)
#'     \item Abundance : Unique ID for biological subject.
#'     \item BioReplicate : Unique ID for technical replicate of one TMT
#'     mixture.
#'     \item Condition : Unique ID for TMT mixture.
#'   }
#'    
#'   }
#'
#' @format A list of two lists with four data.tables.
#' @examples
#' head(summary.data.tmt)
#'
"summary.data.tmt"