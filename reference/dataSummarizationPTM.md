# Data summarization function for label-free MS experiments targeting PTMs.

Utilizes functionality from MSstats to clean, summarize, and normalize
PTM and protein level data. Imputes missing values, performs
normalization, and summarizes data. PTM data is summarized up to the
modification and protein data is summarized up to the protein level.
Takes as input the output of the included converters (see included
`raw.input` data object for required input format).

## Usage

``` r
dataSummarizationPTM(
  data,
  logTrans = 2,
  normalization = "equalizeMedians",
  normalization.PTM = "equalizeMedians",
  nameStandards = NULL,
  nameStandards.PTM = NULL,
  featureSubset = "all",
  featureSubset.PTM = "all",
  remove_uninformative_feature_outlier = FALSE,
  remove_uninformative_feature_outlier.PTM = FALSE,
  min_feature_count = 2,
  min_feature_count.PTM = 1,
  n_top_feature = 3,
  n_top_feature.PTM = 3,
  summaryMethod = "TMP",
  equalFeatureVar = TRUE,
  censoredInt = "NA",
  MBimpute = TRUE,
  MBimpute.PTM = TRUE,
  remove50missing = FALSE,
  fix_missing = NULL,
  maxQuantileforCensored = 0.999,
  use_log_file = TRUE,
  append = TRUE,
  verbose = TRUE,
  log_file_path = NULL,
  base = "MSstatsPTM_log_"
)
```

## Arguments

- data:

  name of the list with PTM and (optionally) unmodified protein
  data.tables, which can be the output of the MSstatsPTM converter
  functions

- logTrans:

  logarithm transformation with base 2(default) or 10

- normalization:

  normalization for the protein level dataset, to remove systematic bias
  between MS runs. There are three different normalizations supported.
  'equalizeMedians'(default) represents constant normalization
  (equalizing the medians) based on reference signals is performed.
  'quantile' represents quantile normalization based on reference
  signals is performed. 'globalStandards' represents normalization with
  global standards proteins. FALSE represents no normalization is
  performed

- normalization.PTM:

  normalization for PTM level dataset. Default is "equalizeMedians" Can
  be adjusted to any of the options described above.

- nameStandards:

  vector of global standard peptide names for protein dataset. only for
  normalization with global standard peptides.

- nameStandards.PTM:

  Same as above for PTM dataset.

- featureSubset:

  "all" (default) uses all features that the data set has. "top3" uses
  top 3 features which have highest average of log-intensity across
  runs. "topN" uses top N features which has highest average of
  log-intensity across runs. It needs the input for n_top_feature
  option. "highQuality" flags uninformative feature and outliers.

- featureSubset.PTM:

  For PTM dataset only. Options same as above.

- remove_uninformative_feature_outlier:

  For protein dataset only. It only works after users used
  featureSubset="highQuality" in dataProcess. TRUE allows to remove 1)
  the features are flagged in the column,
  feature_quality="Uninformative" which are features with bad
  quality, 2) outliers that are flagged in the column, is_outlier=TRUE,
  for run-level summarization. FALSE (default) uses all features and
  intensities for run-level summarization.

- remove_uninformative_feature_outlier.PTM:

  For PTM dataset only. Options same as above.

- min_feature_count:

  optional. Only required if featureSubset = "highQuality". Defines a
  minimum number of informative features a protein needs to be
  considered in the feature selection algorithm.

- min_feature_count.PTM:

  For PTM dataset only. Options the same as above. Default is 1 due to
  low average feature count for PTMs.

- n_top_feature:

  For protein dataset only. The number of top features for
  featureSubset='topN'. Default is 3, which means to use top 3 features.

- n_top_feature.PTM:

  For PTM dataset only. Options same as above.

- summaryMethod:

  "TMP"(default) means Tukey's median polish, which is robust estimation
  method. "linear" uses linear mixed model.

- equalFeatureVar:

  only for summaryMethod="linear". default is TRUE. Logical variable for
  whether the model should account for heterogeneous variation among
  intensities from different features. Default is TRUE, which assume
  equal variance among intensities from features. FALSE means that we
  cannot assume equal variance among intensities from features, then we
  will account for heterogeneous variation from different features.

- censoredInt:

  Missing values are censored or at random. 'NA' (default) assumes that
  all 'NA's in 'Intensity' column are censored. '0' uses zero
  intensities as censored intensity. In this case, NA intensities are
  missing at random. The output from Skyline should use '0'. Null
  assumes that all NA intensites are randomly missing.

- MBimpute:

  For protein dataset only. only for summaryMethod="TMP" and
  censoredInt='NA' or '0'. TRUE (default) imputes 'NA' or '0' (depending
  on censoredInt option) by Accelated failure model. FALSE uses the
  values assigned by cutoffCensored.

- MBimpute.PTM:

  For PTM dataset only. Options same as above.

- remove50missing:

  only for summaryMethod="TMP". TRUE removes the runs which have more
  than 50% missing values. FALSE is default.

- fix_missing:

  Default is Null. Optional, same as the 'fix_missing' parameter in
  MSstatsConvert::MSstatsBalancedDesign function

- maxQuantileforCensored:

  Maximum quantile for deciding censored missing values. default is
  0.999

- use_log_file:

  logical. If TRUE, information about data processing will be saved to a
  file.

- append:

  logical. If TRUE, information about data processing will be added to
  an existing log file.

- verbose:

  logical. If TRUE, information about data processing will be printed to
  the console.

- log_file_path:

  character. Path to a file to which information about data processing
  will be saved. If not provided, such a file will be created
  automatically. If `append = TRUE`, has to be a valid path to a file.

- base:

  start of the file name.

## Value

list of summarized PTM and Protein results. These results contain the
reformatted input to the summarization function, as well as run-level
summarization results.

## Examples

``` r
head(raw.input$PTM)
#> # A tibble: 6 × 10
#>   ProteinName PeptideSequence Condition BioReplicate Run        Intensity
#>   <chr>       <chr>           <chr>     <chr>        <chr>          <dbl>
#> 1 Q9UHD8_K262 DAGLK*QAPASR    CCCP      BCH1         CCCP-B1T1   1423906.
#> 2 Q9UHD8_K262 DAGLK*QAPASR    CCCP      BCH1         CCCP-B1T2    877045.
#> 3 Q9UHD8_K262 DAGLK*QAPASR    CCCP      BCH2         CCCP-B2T1    384418.
#> 4 Q9UHD8_K262 DAGLK*QAPASR    CCCP      BCH2         CCCP-B2T2    454858.
#> 5 Q9UHD8_K262 DAGLK*QAPASR    Combo     BCH1         Combo-B1T1  1603377.
#> 6 Q9UHD8_K262 DAGLK*QAPASR    Combo     BCH1         Combo-B1T2   676555.
#> # ℹ 4 more variables: PrecursorCharge <chr>, FragmentIon <lgl>,
#> #   ProductCharge <lgl>, IsotopeLabelType <chr>
head(raw.input$PROTEIN)
#> # A tibble: 6 × 10
#>   ProteinName PeptideSequence Condition BioReplicate Run           Intensity
#>   <chr>       <chr>           <chr>     <chr>        <chr>             <dbl>
#> 1 Q9UHD8      STLINTLFK       CCCP      BCH2         CCCP-B2T1       367944.
#> 2 Q9UHD8      STLINTLFK       CCCP      BCH2         CCCP-B2T2       341207.
#> 3 Q9UHD8      STLINTLFK       Combo     BCH2         Combo-B2T1      185843.
#> 4 Q9UHD8      STLINTLFK       Ctrl      BCH2         Ctrl-B2T1       529224.
#> 5 Q9UHD8      STLINTLFK       Ctrl      BCH2         Ctrl-B2T2       483355.
#> 6 Q9UHD8      STLINTLFK       USP30_OE  BCH2         USP30_OE-B2T1   447795.
#> # ℹ 4 more variables: PrecursorCharge <chr>, FragmentIon <lgl>,
#> #   ProductCharge <lgl>, IsotopeLabelType <chr>

quant.lf.msstatsptm = dataSummarizationPTM(raw.input, verbose = FALSE)
#> Starting PTM summarization...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%
#> Convergence warning caught: Ran out of iterations and did not converge
#> Warning: Ran out of iterations and did not converge
#>   |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
#> Starting Protein summarization...
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |=====                                                                 |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  15%  |                                                                              |=============                                                         |  19%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  35%  |                                                                              |===========================                                           |  38%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  58%  |                                                                              |===========================================                           |  62%  |                                                                              |==============================================                        |  65%  |                                                                              |================================================                      |  69%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |=========================================================             |  81%  |                                                                              |===========================================================           |  85%  |                                                                              |==============================================================        |  88%  |                                                                              |=================================================================     |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
head(quant.lf.msstatsptm$PTM$ProteinLevelData)
#>   RUN     Protein LogIntensities originalRUN GROUP SUBJECT
#> 1   3 Q9UHD8_K028       20.40124   CCCP-B2T1  CCCP    BCH2
#> 2   4 Q9UHD8_K028       20.48536   CCCP-B2T2  CCCP    BCH2
#> 3   7 Q9UHD8_K028       20.64447  Combo-B2T1 Combo    BCH2
#> 4   8 Q9UHD8_K028       20.73128  Combo-B2T2 Combo    BCH2
#> 5  11 Q9UHD8_K028       20.41225   Ctrl-B2T1  Ctrl    BCH2
#> 6  12 Q9UHD8_K028       20.65940   Ctrl-B2T2  Ctrl    BCH2
#>   TotalGroupMeasurements NumMeasuredFeature MissingPercentage more50missing
#> 1                      4                  1                 0         FALSE
#> 2                      4                  1                 0         FALSE
#> 3                      4                  1                 0         FALSE
#> 4                      4                  1                 0         FALSE
#> 5                      4                  1                 0         FALSE
#> 6                      4                  1                 0         FALSE
#>   NumImputedFeature
#> 1                 0
#> 2                 0
#> 3                 0
#> 4                 0
#> 5                 0
#> 6                 0
```
