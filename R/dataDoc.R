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
#'     \item PROTEIN : MS run ID
#'     \item PEPTIDE : Protein ID with modification site mapped in. Ex.
#'     Protein_1002_S836
#'     \item TRANSITION: Protein-level summarized abundance
#'     \item FEATURE : Labeling information (126, ... 131)
#'     \item LABEL : Condition (ex. Healthy, Cancer, Time0)
#'     \item GROUP : Unique ID for biological subject.
#'     \item RUN : Unique ID for technical replicate of one TMT
#'     mixture.
#'     \item SUBJECT : Unique ID for TMT mixture.
#'     \item FRACTION : Unique ID for TMT mixture.
#'     \item originalRUN : Unique ID for TMT mixture.
#'     \item censored : Unique ID for TMT mixture.
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