# Generate MSstatsPTM required input format from Protein Prospector output

Generate MSstatsPTM required input format from Protein Prospector output

## Usage

``` r
ProteinProspectortoMSstatsPTMFormat(
  input,
  annotation,
  input_protein = NULL,
  annotation_protein = NULL,
  use_unmod_peptides = FALSE,
  mod_ids = c("Phospho"),
  useUniquePeptide = TRUE,
  removeFewMeasurements = TRUE,
  removeProtein_with1Feature = FALSE,
  summaryforMultipleRows = sum,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = NULL
)
```

## Arguments

- input:

  Input txt peptide report file from Protein Prospector with "Keep
  Replicates", "Mods in Peptide", and "Protein Mods" options selected.

- annotation:

  data frame which contains column Run, Fraction, TechRepMixture,
  Mixture, Channel, BioReplicate, Condition.

- input_protein:

  same as `input` for global profiling run. Default is NULL.

- annotation_protein:

  same as `annotation` for global profiling run. Default is NULL.

- use_unmod_peptides:

  If `protein_input` is not provided, unmodified peptides can be
  extracted from `input` to be used in place of a global profiling run.
  Default is `FALSE`.

- mod_ids:

  List of modifications of interest. Default is a list with only
  `Phospho`. Please note that the 'mod_ids' parameter currently supports
  lists of size 1 only. Future updates aim to extend its functionality
  to accommodate lists of greater sizes.

- useUniquePeptide:

  TRUE (default) removes peptides that are assigned for more than one
  proteins. We assume to use unique peptide for each protein.

- removeFewMeasurements:

  TRUE (default) will remove the features that have 1 or 2 measurements
  across runs.

- removeProtein_with1Feature:

  TRUE will remove the proteins which have only 1 feature, which is the
  combination of peptide, precursor charge, fragment and charge. FALSE
  is default.

- summaryforMultipleRows:

  sum(default) or max - when there are multiple measurements for certain
  feature and certain run, use highest or sum of multiple intensities.

- use_log_file:

  logical. If TRUE, information about data processing will be saved to a
  file.

- append:

  logical. If TRUE, information about data processing will be added to
  an existing log file.

- verbose:

  logical. If TRUE, information about data processing wil be printed to
  the console.

- log_file_path:

  character. Path to a file to which information about data processing
  will be saved. If not provided, such a file will be created
  automatically. If 'append = TRUE', has to be a valid path to a file.

## Value

a list of two data.tables named 'PTM' and 'PROTEIN' in the format
required by MSstatsPTM.

## Author

Anthony Wu

## Examples

``` r
input = system.file("tinytest/raw_data/ProteinProspector/Prospector_PhosphoTMT.txt",
    package = "MSstatsPTM")
input = data.table::fread(input)
annot = system.file("tinytest/raw_data/ProteinProspector/Annotation.csv",
                                package = "MSstatsPTM")
annot = data.table::fread(annot)
input_protein = system.file("tinytest/raw_data/ProteinProspector/Prospector_TotalTMT.txt",
    package = "MSstatsConvert")
input_protein = data.table::fread(input_protein)
annot_protein = system.file("tinytest/raw_data/ProteinProspector/Annotation.csv",
                                package = "MSstatsConvert")
annot_protein = data.table::fread(annot_protein)
output <- ProteinProspectortoMSstatsPTMFormat(
    input, 
    annot, 
    input_protein,
    annot_protein
)
#> INFO  [2026-04-25 15:36:38] ** Raw data from ProteinProspector imported successfully.
#> INFO  [2026-04-25 15:36:38] ** Raw data from ProteinProspector cleaned successfully.
#> INFO  [2026-04-25 15:36:38] ** Using provided annotation.
#> INFO  [2026-04-25 15:36:38] ** Run and Channel labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-25 15:36:38] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements within each run will be removed.
#> INFO  [2026-04-25 15:36:38] ** Features with all missing measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:38] ** Shared peptides are removed.
#> INFO  [2026-04-25 15:36:38] ** Features with one or two measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:38] ** PSMs have been aggregated to peptide ions.
#> INFO  [2026-04-25 15:36:38] ** Run annotation merged with quantification data.
#> INFO  [2026-04-25 15:36:38] ** For peptides overlapped between fractions of Mixture1_1 use the fraction with maximal average abundance.
#> INFO  [2026-04-25 15:36:38] ** Fractions belonging to same mixture have been combined.
#> INFO  [2026-04-25 15:36:38] ** Features with one or two measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:38] ** Fractionation handled.
#> INFO  [2026-04-25 15:36:38] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-25 15:36:38] ** Finished preprocessing. The dataset is ready to be processed by the proteinSummarization function.
#> INFO  [2026-04-25 15:36:38] ** Raw data from ProteinProspector imported successfully.
#> INFO  [2026-04-25 15:36:38] ** Raw data from ProteinProspector cleaned successfully.
#> INFO  [2026-04-25 15:36:38] ** Using provided annotation.
#> INFO  [2026-04-25 15:36:38] ** Run and Channel labels were standardized to remove symbols such as '.' or '%'.
#> INFO  [2026-04-25 15:36:38] ** The following options are used:
#>   - Features will be defined by the columns: PeptideSequence, PrecursorCharge
#>   - Shared peptides will be removed.
#>   - Proteins with single feature will not be removed.
#>   - Features with less than 3 measurements within each run will be removed.
#> INFO  [2026-04-25 15:36:38] ** Features with all missing measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:38] ** Shared peptides are removed.
#> INFO  [2026-04-25 15:36:38] ** Features with one or two measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:38] ** PSMs have been aggregated to peptide ions.
#> INFO  [2026-04-25 15:36:38] ** Run annotation merged with quantification data.
#> INFO  [2026-04-25 15:36:38] ** For peptides overlapped between fractions of Mixture1_1 use the fraction with maximal average abundance.
#> INFO  [2026-04-25 15:36:38] ** Fractions belonging to same mixture have been combined.
#> INFO  [2026-04-25 15:36:38] ** Features with one or two measurements across channels within each run are removed.
#> INFO  [2026-04-25 15:36:38] ** Fractionation handled.
#> INFO  [2026-04-25 15:36:38] ** Updated quantification data to make balanced design. Missing values are marked by NA
#> INFO  [2026-04-25 15:36:38] ** Finished preprocessing. The dataset is ready to be processed by the proteinSummarization function.
head(output)
#> $PTM
#>          ProteinName                          PeptideSequence Charge
#> 1        Q9QYX7_S608                     AAS(Phospho)VQPATASK      2
#> 2        Q9QYX7_S608                     AAS(Phospho)VQPATASK      2
#> 3        Q9QYX7_S608                     AAS(Phospho)VQPATASK      2
#> 4        Q9QYX7_S608                     AAS(Phospho)VQPATASK      2
#> 5        Q9QYX7_S608                     AAS(Phospho)VQPATASK      2
#> 6        Q9QYX7_S608                     AAS(Phospho)VQPATASK      2
#> 7   Q9QYX7_S617_S625 AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK      4
#> 8   Q9QYX7_S617_S625 AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK      4
#> 9   Q9QYX7_S617_S625 AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK      4
#> 10  Q9QYX7_S617_S625 AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK      4
#> 11  Q9QYX7_S617_S625 AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK      4
#> 12  Q9QYX7_S617_S625 AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK      4
#> 13       Q9QYX7_S625          AASVQPATASKSPVPSQQAS(Phospho)PK      4
#> 14       Q9QYX7_S625          AASVQPATASKSPVPSQQAS(Phospho)PK      4
#> 15       Q9QYX7_S625          AASVQPATASKSPVPSQQAS(Phospho)PK      4
#> 16       Q9QYX7_S625          AASVQPATASKSPVPSQQAS(Phospho)PK      4
#> 17       Q9QYX7_S625          AASVQPATASKSPVPSQQAS(Phospho)PK      4
#> 18       Q9QYX7_S625          AASVQPATASKSPVPSQQAS(Phospho)PK      4
#> 19       P16546_S866                       AKLS(Phospho)ELNQK      3
#> 20       P16546_S866                       AKLS(Phospho)ELNQK      3
#> 21       P16546_S866                       AKLS(Phospho)ELNQK      3
#> 22       P16546_S866                       AKLS(Phospho)ELNQK      3
#> 23       P16546_S866                       AKLS(Phospho)ELNQK      3
#> 24       P16546_S866                       AKLS(Phospho)ELNQK      3
#> 25       Q9QYX7_S599              ALGGELAAIPSS(Phospho)PQPTPK      3
#> 26       Q9QYX7_S599              ALGGELAAIPSS(Phospho)PQPTPK      3
#> 27       Q9QYX7_S599              ALGGELAAIPSS(Phospho)PQPTPK      3
#> 28       Q9QYX7_S599              ALGGELAAIPSS(Phospho)PQPTPK      3
#> 29       Q9QYX7_S599              ALGGELAAIPSS(Phospho)PQPTPK      3
#> 30       Q9QYX7_S599              ALGGELAAIPSS(Phospho)PQPTPK      3
#> 31       P16546_S572                  AQLADS(Phospho)FHLQQFFR      3
#> 32       P16546_S572                  AQLADS(Phospho)FHLQQFFR      3
#> 33       P16546_S572                  AQLADS(Phospho)FHLQQFFR      3
#> 34       P16546_S572                  AQLADS(Phospho)FHLQQFFR      3
#> 35       P16546_S572                  AQLADS(Phospho)FHLQQFFR      3
#> 36       P16546_S572                  AQLADS(Phospho)FHLQQFFR      3
#> 37       Q9QYX7_S289            ATVQQPGPAKS(Phospho)PAQPAGTGK      4
#> 38       Q9QYX7_S289            ATVQQPGPAKS(Phospho)PAQPAGTGK      4
#> 39       Q9QYX7_S289            ATVQQPGPAKS(Phospho)PAQPAGTGK      4
#> 40       Q9QYX7_S289            ATVQQPGPAKS(Phospho)PAQPAGTGK      4
#> 41       Q9QYX7_S289            ATVQQPGPAKS(Phospho)PAQPAGTGK      4
#> 42       Q9QYX7_S289            ATVQQPGPAKS(Phospho)PAQPAGTGK      4
#> 43       P16546_S553                          DALLS(Phospho)R      2
#> 44       P16546_S553                          DALLS(Phospho)R      2
#> 45       P16546_S553                          DALLS(Phospho)R      2
#> 46       P16546_S553                          DALLS(Phospho)R      2
#> 47       P16546_S553                          DALLS(Phospho)R      2
#> 48       P16546_S553                          DALLS(Phospho)R      2
#> 49       P16546_S289                      DLAS(Phospho)VQALLR      2
#> 50       P16546_S289                      DLAS(Phospho)VQALLR      2
#> 51       P16546_S289                      DLAS(Phospho)VQALLR      2
#> 52       P16546_S289                      DLAS(Phospho)VQALLR      2
#> 53       P16546_S289                      DLAS(Phospho)VQALLR      2
#> 54       P16546_S289                      DLAS(Phospho)VQALLR      2
#> 55       Q9QYX7_S212             DQGKSEGITKPSLQQPS(Phospho)PK      5
#> 56       Q9QYX7_S212             DQGKSEGITKPSLQQPS(Phospho)PK      5
#> 57       Q9QYX7_S212             DQGKSEGITKPSLQQPS(Phospho)PK      5
#> 58       Q9QYX7_S212             DQGKSEGITKPSLQQPS(Phospho)PK      5
#> 59       Q9QYX7_S212             DQGKSEGITKPSLQQPS(Phospho)PK      5
#> 60       Q9QYX7_S212             DQGKSEGITKPSLQQPS(Phospho)PK      5
#> 61       P16546_S587                    DSDELKS(Phospho)WVNEK      3
#> 62       P16546_S587                    DSDELKS(Phospho)WVNEK      3
#> 63       P16546_S587                    DSDELKS(Phospho)WVNEK      3
#> 64       P16546_S587                    DSDELKS(Phospho)WVNEK      3
#> 65       P16546_S587                    DSDELKS(Phospho)WVNEK      3
#> 66       P16546_S587                    DSDELKS(Phospho)WVNEK      3
#> 67       P16546_S484                     DTEQVDNWMS(Phospho)K      3
#> 68       P16546_S484                     DTEQVDNWMS(Phospho)K      3
#> 69       P16546_S484                     DTEQVDNWMS(Phospho)K      3
#> 70       P16546_S484                     DTEQVDNWMS(Phospho)K      3
#> 71       P16546_S484                     DTEQVDNWMS(Phospho)K      3
#> 72       P16546_S484                     DTEQVDNWMS(Phospho)K      3
#> 73       P16546_S809                     EKEPIAAS(Phospho)TNR      3
#> 74       P16546_S809                     EKEPIAAS(Phospho)TNR      3
#> 75       P16546_S809                     EKEPIAAS(Phospho)TNR      3
#> 76       P16546_S809                     EKEPIAAS(Phospho)TNR      3
#> 77       P16546_S809                     EKEPIAAS(Phospho)TNR      3
#> 78       P16546_S809                     EKEPIAAS(Phospho)TNR      3
#> 79       Q9QYX7_S235             EVIPQDIPSKS(Phospho)VSSQQAEK      3
#> 80       Q9QYX7_S235             EVIPQDIPSKS(Phospho)VSSQQAEK      3
#> 81       Q9QYX7_S235             EVIPQDIPSKS(Phospho)VSSQQAEK      3
#> 82       Q9QYX7_S235             EVIPQDIPSKS(Phospho)VSSQQAEK      3
#> 83       Q9QYX7_S235             EVIPQDIPSKS(Phospho)VSSQQAEK      3
#> 84       Q9QYX7_S235             EVIPQDIPSKS(Phospho)VSSQQAEK      3
#> 85       Q9QYX7_S177             FNPFDLIS(Phospho)DSEAVQEETTK      3
#> 86       Q9QYX7_S177             FNPFDLIS(Phospho)DSEAVQEETTK      3
#> 87       Q9QYX7_S177             FNPFDLIS(Phospho)DSEAVQEETTK      3
#> 88       Q9QYX7_S177             FNPFDLIS(Phospho)DSEAVQEETTK      3
#> 89       Q9QYX7_S177             FNPFDLIS(Phospho)DSEAVQEETTK      3
#> 90       Q9QYX7_S177             FNPFDLIS(Phospho)DSEAVQEETTK      3
#> 91       Q9QYX7_S179             FNPFDLISDS(Phospho)EAVQEETTK      3
#> 92       Q9QYX7_S179             FNPFDLISDS(Phospho)EAVQEETTK      3
#> 93       Q9QYX7_S179             FNPFDLISDS(Phospho)EAVQEETTK      3
#> 94       Q9QYX7_S179             FNPFDLISDS(Phospho)EAVQEETTK      3
#> 95       Q9QYX7_S179             FNPFDLISDS(Phospho)EAVQEETTK      3
#> 96       Q9QYX7_S179             FNPFDLISDS(Phospho)EAVQEETTK      3
#> 97        Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      2
#> 98        Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      2
#> 99        Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      2
#> 100       Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      2
#> 101       Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      2
#> 102       Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      2
#> 103       Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      3
#> 104       Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      3
#> 105       Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      3
#> 106       Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      3
#> 107       Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      3
#> 108       Q9QYX7_S75                 GSVPAAAAES(Phospho)PSMHR      3
#> 109       Q9QYX7_S77                 GSVPAAAAESPS(Phospho)MHR      3
#> 110       Q9QYX7_S77                 GSVPAAAAESPS(Phospho)MHR      3
#> 111       Q9QYX7_S77                 GSVPAAAAESPS(Phospho)MHR      3
#> 112       Q9QYX7_S77                 GSVPAAAAESPS(Phospho)MHR      3
#> 113       Q9QYX7_S77                 GSVPAAAAESPS(Phospho)MHR      3
#> 114       Q9QYX7_S77                 GSVPAAAAESPS(Phospho)MHR      3
#> 115      Q6PFD5_S965          HSSATESADS(Phospho)IEIYIPEAQTRL      3
#> 116      Q6PFD5_S965          HSSATESADS(Phospho)IEIYIPEAQTRL      3
#> 117      Q6PFD5_S965          HSSATESADS(Phospho)IEIYIPEAQTRL      3
#> 118      Q6PFD5_S965          HSSATESADS(Phospho)IEIYIPEAQTRL      3
#> 119      Q6PFD5_S965          HSSATESADS(Phospho)IEIYIPEAQTRL      3
#> 120      Q6PFD5_S965          HSSATESADS(Phospho)IEIYIPEAQTRL      3
#> 121      Q6PFD5_S965           HSSATESADS(Phospho)IEIYIPEAQTR      3
#> 122      Q6PFD5_S965           HSSATESADS(Phospho)IEIYIPEAQTR      3
#> 123      Q6PFD5_S965           HSSATESADS(Phospho)IEIYIPEAQTR      3
#> 124      Q6PFD5_S965           HSSATESADS(Phospho)IEIYIPEAQTR      3
#> 125      Q6PFD5_S965           HSSATESADS(Phospho)IEIYIPEAQTR      3
#> 126      Q6PFD5_S965           HSSATESADS(Phospho)IEIYIPEAQTR      3
#> 127      Q9QYX7_S636                     KELPSKQDS(Phospho)PK      4
#> 128      Q9QYX7_S636                     KELPSKQDS(Phospho)PK      4
#> 129      Q9QYX7_S636                     KELPSKQDS(Phospho)PK      4
#> 130      Q9QYX7_S636                     KELPSKQDS(Phospho)PK      4
#> 131      Q9QYX7_S636                     KELPSKQDS(Phospho)PK      4
#> 132      Q9QYX7_S636                     KELPSKQDS(Phospho)PK      4
#> 133      Q9QYX7_S128                    LPGRSPS(Phospho)TISLK      3
#> 134      Q9QYX7_S128                    LPGRSPS(Phospho)TISLK      3
#> 135      Q9QYX7_S128                    LPGRSPS(Phospho)TISLK      3
#> 136      Q9QYX7_S128                    LPGRSPS(Phospho)TISLK      3
#> 137      Q9QYX7_S128                    LPGRSPS(Phospho)TISLK      3
#> 138      Q9QYX7_S128                    LPGRSPS(Phospho)TISLK      3
#> 139      P16546_S324                 LQQS(Phospho)HPLSASQIQVK      3
#> 140      P16546_S324                 LQQS(Phospho)HPLSASQIQVK      3
#> 141      P16546_S324                 LQQS(Phospho)HPLSASQIQVK      3
#> 142      P16546_S324                 LQQS(Phospho)HPLSASQIQVK      3
#> 143      P16546_S324                 LQQS(Phospho)HPLSASQIQVK      3
#> 144      P16546_S324                 LQQS(Phospho)HPLSASQIQVK      3
#> 145      P16546_S330                 LQQSHPLSAS(Phospho)QIQVK      3
#> 146      P16546_S330                 LQQSHPLSAS(Phospho)QIQVK      3
#> 147      P16546_S330                 LQQSHPLSAS(Phospho)QIQVK      3
#> 148      P16546_S330                 LQQSHPLSAS(Phospho)QIQVK      3
#> 149      P16546_S330                 LQQSHPLSAS(Phospho)QIQVK      3
#> 150      P16546_S330                 LQQSHPLSAS(Phospho)QIQVK      3
#> 151      Q9QYX7_S636                      QDS(Phospho)PKAPESK      3
#> 152      Q9QYX7_S636                      QDS(Phospho)PKAPESK      3
#> 153      Q9QYX7_S636                      QDS(Phospho)PKAPESK      3
#> 154      Q9QYX7_S636                      QDS(Phospho)PKAPESK      3
#> 155      Q9QYX7_S636                      QDS(Phospho)PKAPESK      3
#> 156      Q9QYX7_S636                      QDS(Phospho)PKAPESK      3
#> 157      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      3
#> 158      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      3
#> 159      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      3
#> 160      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      3
#> 161      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      3
#> 162      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      3
#> 163      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      4
#> 164      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      4
#> 165      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      4
#> 166      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      4
#> 167      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      4
#> 168      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      4
#> 169      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      5
#> 170      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      5
#> 171      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      5
#> 172      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      5
#> 173      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      5
#> 174      Q9QYX7_S212                 SEGITKPSLQQPS(Phospho)PK      5
#> 175 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      3
#> 176 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      3
#> 177 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      3
#> 178 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      3
#> 179 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      3
#> 180 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      3
#> 181 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      4
#> 182 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      4
#> 183 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      4
#> 184 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      4
#> 185 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      4
#> 186 Q9QYX7_S126_S131       SEQKLPGRS(Phospho)PSTIS(Phospho)LK      4
#> 187      Q9QYX7_S128                        SPS(Phospho)TISLK      2
#> 188      Q9QYX7_S128                        SPS(Phospho)TISLK      2
#> 189      Q9QYX7_S128                        SPS(Phospho)TISLK      2
#> 190      Q9QYX7_S128                        SPS(Phospho)TISLK      2
#> 191      Q9QYX7_S128                        SPS(Phospho)TISLK      2
#> 192      Q9QYX7_S128                        SPS(Phospho)TISLK      2
#> 193      Q9QYX7_S131                        SPSTIS(Phospho)LK      2
#> 194      Q9QYX7_S131                        SPSTIS(Phospho)LK      2
#> 195      Q9QYX7_S131                        SPSTIS(Phospho)LK      2
#> 196      Q9QYX7_S131                        SPSTIS(Phospho)LK      2
#> 197      Q9QYX7_S131                        SPSTIS(Phospho)LK      2
#> 198      Q9QYX7_S131                        SPSTIS(Phospho)LK      2
#> 199      Q9QYX7_S131                        SPSTIS(Phospho)LK      3
#> 200      Q9QYX7_S131                        SPSTIS(Phospho)LK      3
#> 201      Q9QYX7_S131                        SPSTIS(Phospho)LK      3
#> 202      Q9QYX7_S131                        SPSTIS(Phospho)LK      3
#> 203      Q9QYX7_S131                        SPSTIS(Phospho)LK      3
#> 204      Q9QYX7_S131                        SPSTIS(Phospho)LK      3
#> 205      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      3
#> 206      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      3
#> 207      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      3
#> 208      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      3
#> 209      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      3
#> 210      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      3
#> 211      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      4
#> 212      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      4
#> 213      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      4
#> 214      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      4
#> 215      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      4
#> 216      Q9QYX7_S625                    SPVPSQQAS(Phospho)PKK      4
#> 217      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      2
#> 218      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      2
#> 219      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      2
#> 220      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      2
#> 221      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      2
#> 222      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      2
#> 223      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      3
#> 224      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      3
#> 225      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      3
#> 226      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      3
#> 227      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      3
#> 228      Q9QYX7_S625                     SPVPSQQAS(Phospho)PK      3
#> 229      Q9QYX7_T112                        SRT(Phospho)TDTFR      3
#> 230      Q9QYX7_T112                        SRT(Phospho)TDTFR      3
#> 231      Q9QYX7_T112                        SRT(Phospho)TDTFR      3
#> 232      Q9QYX7_T112                        SRT(Phospho)TDTFR      3
#> 233      Q9QYX7_T112                        SRT(Phospho)TDTFR      3
#> 234      Q9QYX7_T112                        SRT(Phospho)TDTFR      3
#> 235      Q9QYX7_S165         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK      3
#> 236      Q9QYX7_S165         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK      3
#> 237      Q9QYX7_S165         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK      3
#> 238      Q9QYX7_S165         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK      3
#> 239      Q9QYX7_S165         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK      3
#> 240      Q9QYX7_S165         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK      3
#> 241      Q9QYX7_S237                       SVS(Phospho)SQQAEK      2
#> 242      Q9QYX7_S237                       SVS(Phospho)SQQAEK      2
#> 243      Q9QYX7_S237                       SVS(Phospho)SQQAEK      2
#> 244      Q9QYX7_S237                       SVS(Phospho)SQQAEK      2
#> 245      Q9QYX7_S237                       SVS(Phospho)SQQAEK      2
#> 246      Q9QYX7_S237                       SVS(Phospho)SQQAEK      2
#> 247      Q9QYX7_S238                       SVSS(Phospho)QQAEK      2
#> 248      Q9QYX7_S238                       SVSS(Phospho)QQAEK      2
#> 249      Q9QYX7_S238                       SVSS(Phospho)QQAEK      2
#> 250      Q9QYX7_S238                       SVSS(Phospho)QQAEK      2
#> 251      Q9QYX7_S238                       SVSS(Phospho)QQAEK      2
#> 252      Q9QYX7_S238                       SVSS(Phospho)QQAEK      2
#>                                            PSM  Mixture TechRepMixture
#> 1                       AAS(Phospho)VQPATASK_2 Mixture1              1
#> 2                       AAS(Phospho)VQPATASK_2 Mixture1              1
#> 3                       AAS(Phospho)VQPATASK_2 Mixture1              1
#> 4                       AAS(Phospho)VQPATASK_2 Mixture1              1
#> 5                       AAS(Phospho)VQPATASK_2 Mixture1              1
#> 6                       AAS(Phospho)VQPATASK_2 Mixture1              1
#> 7   AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK_4 Mixture1              1
#> 8   AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK_4 Mixture1              1
#> 9   AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK_4 Mixture1              1
#> 10  AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK_4 Mixture1              1
#> 11  AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK_4 Mixture1              1
#> 12  AASVQPATASKS(Phospho)PVPSQQAS(Phospho)PK_4 Mixture1              1
#> 13           AASVQPATASKSPVPSQQAS(Phospho)PK_4 Mixture1              1
#> 14           AASVQPATASKSPVPSQQAS(Phospho)PK_4 Mixture1              1
#> 15           AASVQPATASKSPVPSQQAS(Phospho)PK_4 Mixture1              1
#> 16           AASVQPATASKSPVPSQQAS(Phospho)PK_4 Mixture1              1
#> 17           AASVQPATASKSPVPSQQAS(Phospho)PK_4 Mixture1              1
#> 18           AASVQPATASKSPVPSQQAS(Phospho)PK_4 Mixture1              1
#> 19                        AKLS(Phospho)ELNQK_3 Mixture1              1
#> 20                        AKLS(Phospho)ELNQK_3 Mixture1              1
#> 21                        AKLS(Phospho)ELNQK_3 Mixture1              1
#> 22                        AKLS(Phospho)ELNQK_3 Mixture1              1
#> 23                        AKLS(Phospho)ELNQK_3 Mixture1              1
#> 24                        AKLS(Phospho)ELNQK_3 Mixture1              1
#> 25               ALGGELAAIPSS(Phospho)PQPTPK_3 Mixture1              1
#> 26               ALGGELAAIPSS(Phospho)PQPTPK_3 Mixture1              1
#> 27               ALGGELAAIPSS(Phospho)PQPTPK_3 Mixture1              1
#> 28               ALGGELAAIPSS(Phospho)PQPTPK_3 Mixture1              1
#> 29               ALGGELAAIPSS(Phospho)PQPTPK_3 Mixture1              1
#> 30               ALGGELAAIPSS(Phospho)PQPTPK_3 Mixture1              1
#> 31                   AQLADS(Phospho)FHLQQFFR_3 Mixture1              1
#> 32                   AQLADS(Phospho)FHLQQFFR_3 Mixture1              1
#> 33                   AQLADS(Phospho)FHLQQFFR_3 Mixture1              1
#> 34                   AQLADS(Phospho)FHLQQFFR_3 Mixture1              1
#> 35                   AQLADS(Phospho)FHLQQFFR_3 Mixture1              1
#> 36                   AQLADS(Phospho)FHLQQFFR_3 Mixture1              1
#> 37             ATVQQPGPAKS(Phospho)PAQPAGTGK_4 Mixture1              1
#> 38             ATVQQPGPAKS(Phospho)PAQPAGTGK_4 Mixture1              1
#> 39             ATVQQPGPAKS(Phospho)PAQPAGTGK_4 Mixture1              1
#> 40             ATVQQPGPAKS(Phospho)PAQPAGTGK_4 Mixture1              1
#> 41             ATVQQPGPAKS(Phospho)PAQPAGTGK_4 Mixture1              1
#> 42             ATVQQPGPAKS(Phospho)PAQPAGTGK_4 Mixture1              1
#> 43                           DALLS(Phospho)R_2 Mixture1              1
#> 44                           DALLS(Phospho)R_2 Mixture1              1
#> 45                           DALLS(Phospho)R_2 Mixture1              1
#> 46                           DALLS(Phospho)R_2 Mixture1              1
#> 47                           DALLS(Phospho)R_2 Mixture1              1
#> 48                           DALLS(Phospho)R_2 Mixture1              1
#> 49                       DLAS(Phospho)VQALLR_2 Mixture1              1
#> 50                       DLAS(Phospho)VQALLR_2 Mixture1              1
#> 51                       DLAS(Phospho)VQALLR_2 Mixture1              1
#> 52                       DLAS(Phospho)VQALLR_2 Mixture1              1
#> 53                       DLAS(Phospho)VQALLR_2 Mixture1              1
#> 54                       DLAS(Phospho)VQALLR_2 Mixture1              1
#> 55              DQGKSEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 56              DQGKSEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 57              DQGKSEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 58              DQGKSEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 59              DQGKSEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 60              DQGKSEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 61                     DSDELKS(Phospho)WVNEK_3 Mixture1              1
#> 62                     DSDELKS(Phospho)WVNEK_3 Mixture1              1
#> 63                     DSDELKS(Phospho)WVNEK_3 Mixture1              1
#> 64                     DSDELKS(Phospho)WVNEK_3 Mixture1              1
#> 65                     DSDELKS(Phospho)WVNEK_3 Mixture1              1
#> 66                     DSDELKS(Phospho)WVNEK_3 Mixture1              1
#> 67                      DTEQVDNWMS(Phospho)K_3 Mixture1              1
#> 68                      DTEQVDNWMS(Phospho)K_3 Mixture1              1
#> 69                      DTEQVDNWMS(Phospho)K_3 Mixture1              1
#> 70                      DTEQVDNWMS(Phospho)K_3 Mixture1              1
#> 71                      DTEQVDNWMS(Phospho)K_3 Mixture1              1
#> 72                      DTEQVDNWMS(Phospho)K_3 Mixture1              1
#> 73                      EKEPIAAS(Phospho)TNR_3 Mixture1              1
#> 74                      EKEPIAAS(Phospho)TNR_3 Mixture1              1
#> 75                      EKEPIAAS(Phospho)TNR_3 Mixture1              1
#> 76                      EKEPIAAS(Phospho)TNR_3 Mixture1              1
#> 77                      EKEPIAAS(Phospho)TNR_3 Mixture1              1
#> 78                      EKEPIAAS(Phospho)TNR_3 Mixture1              1
#> 79              EVIPQDIPSKS(Phospho)VSSQQAEK_3 Mixture1              1
#> 80              EVIPQDIPSKS(Phospho)VSSQQAEK_3 Mixture1              1
#> 81              EVIPQDIPSKS(Phospho)VSSQQAEK_3 Mixture1              1
#> 82              EVIPQDIPSKS(Phospho)VSSQQAEK_3 Mixture1              1
#> 83              EVIPQDIPSKS(Phospho)VSSQQAEK_3 Mixture1              1
#> 84              EVIPQDIPSKS(Phospho)VSSQQAEK_3 Mixture1              1
#> 85              FNPFDLIS(Phospho)DSEAVQEETTK_3 Mixture1              1
#> 86              FNPFDLIS(Phospho)DSEAVQEETTK_3 Mixture1              1
#> 87              FNPFDLIS(Phospho)DSEAVQEETTK_3 Mixture1              1
#> 88              FNPFDLIS(Phospho)DSEAVQEETTK_3 Mixture1              1
#> 89              FNPFDLIS(Phospho)DSEAVQEETTK_3 Mixture1              1
#> 90              FNPFDLIS(Phospho)DSEAVQEETTK_3 Mixture1              1
#> 91              FNPFDLISDS(Phospho)EAVQEETTK_3 Mixture1              1
#> 92              FNPFDLISDS(Phospho)EAVQEETTK_3 Mixture1              1
#> 93              FNPFDLISDS(Phospho)EAVQEETTK_3 Mixture1              1
#> 94              FNPFDLISDS(Phospho)EAVQEETTK_3 Mixture1              1
#> 95              FNPFDLISDS(Phospho)EAVQEETTK_3 Mixture1              1
#> 96              FNPFDLISDS(Phospho)EAVQEETTK_3 Mixture1              1
#> 97                  GSVPAAAAES(Phospho)PSMHR_2 Mixture1              1
#> 98                  GSVPAAAAES(Phospho)PSMHR_2 Mixture1              1
#> 99                  GSVPAAAAES(Phospho)PSMHR_2 Mixture1              1
#> 100                 GSVPAAAAES(Phospho)PSMHR_2 Mixture1              1
#> 101                 GSVPAAAAES(Phospho)PSMHR_2 Mixture1              1
#> 102                 GSVPAAAAES(Phospho)PSMHR_2 Mixture1              1
#> 103                 GSVPAAAAES(Phospho)PSMHR_3 Mixture1              1
#> 104                 GSVPAAAAES(Phospho)PSMHR_3 Mixture1              1
#> 105                 GSVPAAAAES(Phospho)PSMHR_3 Mixture1              1
#> 106                 GSVPAAAAES(Phospho)PSMHR_3 Mixture1              1
#> 107                 GSVPAAAAES(Phospho)PSMHR_3 Mixture1              1
#> 108                 GSVPAAAAES(Phospho)PSMHR_3 Mixture1              1
#> 109                 GSVPAAAAESPS(Phospho)MHR_3 Mixture1              1
#> 110                 GSVPAAAAESPS(Phospho)MHR_3 Mixture1              1
#> 111                 GSVPAAAAESPS(Phospho)MHR_3 Mixture1              1
#> 112                 GSVPAAAAESPS(Phospho)MHR_3 Mixture1              1
#> 113                 GSVPAAAAESPS(Phospho)MHR_3 Mixture1              1
#> 114                 GSVPAAAAESPS(Phospho)MHR_3 Mixture1              1
#> 115          HSSATESADS(Phospho)IEIYIPEAQTRL_3 Mixture1              1
#> 116          HSSATESADS(Phospho)IEIYIPEAQTRL_3 Mixture1              1
#> 117          HSSATESADS(Phospho)IEIYIPEAQTRL_3 Mixture1              1
#> 118          HSSATESADS(Phospho)IEIYIPEAQTRL_3 Mixture1              1
#> 119          HSSATESADS(Phospho)IEIYIPEAQTRL_3 Mixture1              1
#> 120          HSSATESADS(Phospho)IEIYIPEAQTRL_3 Mixture1              1
#> 121           HSSATESADS(Phospho)IEIYIPEAQTR_3 Mixture1              1
#> 122           HSSATESADS(Phospho)IEIYIPEAQTR_3 Mixture1              1
#> 123           HSSATESADS(Phospho)IEIYIPEAQTR_3 Mixture1              1
#> 124           HSSATESADS(Phospho)IEIYIPEAQTR_3 Mixture1              1
#> 125           HSSATESADS(Phospho)IEIYIPEAQTR_3 Mixture1              1
#> 126           HSSATESADS(Phospho)IEIYIPEAQTR_3 Mixture1              1
#> 127                     KELPSKQDS(Phospho)PK_4 Mixture1              1
#> 128                     KELPSKQDS(Phospho)PK_4 Mixture1              1
#> 129                     KELPSKQDS(Phospho)PK_4 Mixture1              1
#> 130                     KELPSKQDS(Phospho)PK_4 Mixture1              1
#> 131                     KELPSKQDS(Phospho)PK_4 Mixture1              1
#> 132                     KELPSKQDS(Phospho)PK_4 Mixture1              1
#> 133                    LPGRSPS(Phospho)TISLK_3 Mixture1              1
#> 134                    LPGRSPS(Phospho)TISLK_3 Mixture1              1
#> 135                    LPGRSPS(Phospho)TISLK_3 Mixture1              1
#> 136                    LPGRSPS(Phospho)TISLK_3 Mixture1              1
#> 137                    LPGRSPS(Phospho)TISLK_3 Mixture1              1
#> 138                    LPGRSPS(Phospho)TISLK_3 Mixture1              1
#> 139                 LQQS(Phospho)HPLSASQIQVK_3 Mixture1              1
#> 140                 LQQS(Phospho)HPLSASQIQVK_3 Mixture1              1
#> 141                 LQQS(Phospho)HPLSASQIQVK_3 Mixture1              1
#> 142                 LQQS(Phospho)HPLSASQIQVK_3 Mixture1              1
#> 143                 LQQS(Phospho)HPLSASQIQVK_3 Mixture1              1
#> 144                 LQQS(Phospho)HPLSASQIQVK_3 Mixture1              1
#> 145                 LQQSHPLSAS(Phospho)QIQVK_3 Mixture1              1
#> 146                 LQQSHPLSAS(Phospho)QIQVK_3 Mixture1              1
#> 147                 LQQSHPLSAS(Phospho)QIQVK_3 Mixture1              1
#> 148                 LQQSHPLSAS(Phospho)QIQVK_3 Mixture1              1
#> 149                 LQQSHPLSAS(Phospho)QIQVK_3 Mixture1              1
#> 150                 LQQSHPLSAS(Phospho)QIQVK_3 Mixture1              1
#> 151                      QDS(Phospho)PKAPESK_3 Mixture1              1
#> 152                      QDS(Phospho)PKAPESK_3 Mixture1              1
#> 153                      QDS(Phospho)PKAPESK_3 Mixture1              1
#> 154                      QDS(Phospho)PKAPESK_3 Mixture1              1
#> 155                      QDS(Phospho)PKAPESK_3 Mixture1              1
#> 156                      QDS(Phospho)PKAPESK_3 Mixture1              1
#> 157                 SEGITKPSLQQPS(Phospho)PK_3 Mixture1              1
#> 158                 SEGITKPSLQQPS(Phospho)PK_3 Mixture1              1
#> 159                 SEGITKPSLQQPS(Phospho)PK_3 Mixture1              1
#> 160                 SEGITKPSLQQPS(Phospho)PK_3 Mixture1              1
#> 161                 SEGITKPSLQQPS(Phospho)PK_3 Mixture1              1
#> 162                 SEGITKPSLQQPS(Phospho)PK_3 Mixture1              1
#> 163                 SEGITKPSLQQPS(Phospho)PK_4 Mixture1              1
#> 164                 SEGITKPSLQQPS(Phospho)PK_4 Mixture1              1
#> 165                 SEGITKPSLQQPS(Phospho)PK_4 Mixture1              1
#> 166                 SEGITKPSLQQPS(Phospho)PK_4 Mixture1              1
#> 167                 SEGITKPSLQQPS(Phospho)PK_4 Mixture1              1
#> 168                 SEGITKPSLQQPS(Phospho)PK_4 Mixture1              1
#> 169                 SEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 170                 SEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 171                 SEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 172                 SEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 173                 SEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 174                 SEGITKPSLQQPS(Phospho)PK_5 Mixture1              1
#> 175       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_3 Mixture1              1
#> 176       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_3 Mixture1              1
#> 177       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_3 Mixture1              1
#> 178       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_3 Mixture1              1
#> 179       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_3 Mixture1              1
#> 180       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_3 Mixture1              1
#> 181       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_4 Mixture1              1
#> 182       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_4 Mixture1              1
#> 183       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_4 Mixture1              1
#> 184       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_4 Mixture1              1
#> 185       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_4 Mixture1              1
#> 186       SEQKLPGRS(Phospho)PSTIS(Phospho)LK_4 Mixture1              1
#> 187                        SPS(Phospho)TISLK_2 Mixture1              1
#> 188                        SPS(Phospho)TISLK_2 Mixture1              1
#> 189                        SPS(Phospho)TISLK_2 Mixture1              1
#> 190                        SPS(Phospho)TISLK_2 Mixture1              1
#> 191                        SPS(Phospho)TISLK_2 Mixture1              1
#> 192                        SPS(Phospho)TISLK_2 Mixture1              1
#> 193                        SPSTIS(Phospho)LK_2 Mixture1              1
#> 194                        SPSTIS(Phospho)LK_2 Mixture1              1
#> 195                        SPSTIS(Phospho)LK_2 Mixture1              1
#> 196                        SPSTIS(Phospho)LK_2 Mixture1              1
#> 197                        SPSTIS(Phospho)LK_2 Mixture1              1
#> 198                        SPSTIS(Phospho)LK_2 Mixture1              1
#> 199                        SPSTIS(Phospho)LK_3 Mixture1              1
#> 200                        SPSTIS(Phospho)LK_3 Mixture1              1
#> 201                        SPSTIS(Phospho)LK_3 Mixture1              1
#> 202                        SPSTIS(Phospho)LK_3 Mixture1              1
#> 203                        SPSTIS(Phospho)LK_3 Mixture1              1
#> 204                        SPSTIS(Phospho)LK_3 Mixture1              1
#> 205                    SPVPSQQAS(Phospho)PKK_3 Mixture1              1
#> 206                    SPVPSQQAS(Phospho)PKK_3 Mixture1              1
#> 207                    SPVPSQQAS(Phospho)PKK_3 Mixture1              1
#> 208                    SPVPSQQAS(Phospho)PKK_3 Mixture1              1
#> 209                    SPVPSQQAS(Phospho)PKK_3 Mixture1              1
#> 210                    SPVPSQQAS(Phospho)PKK_3 Mixture1              1
#> 211                    SPVPSQQAS(Phospho)PKK_4 Mixture1              1
#> 212                    SPVPSQQAS(Phospho)PKK_4 Mixture1              1
#> 213                    SPVPSQQAS(Phospho)PKK_4 Mixture1              1
#> 214                    SPVPSQQAS(Phospho)PKK_4 Mixture1              1
#> 215                    SPVPSQQAS(Phospho)PKK_4 Mixture1              1
#> 216                    SPVPSQQAS(Phospho)PKK_4 Mixture1              1
#> 217                     SPVPSQQAS(Phospho)PK_2 Mixture1              1
#> 218                     SPVPSQQAS(Phospho)PK_2 Mixture1              1
#> 219                     SPVPSQQAS(Phospho)PK_2 Mixture1              1
#> 220                     SPVPSQQAS(Phospho)PK_2 Mixture1              1
#> 221                     SPVPSQQAS(Phospho)PK_2 Mixture1              1
#> 222                     SPVPSQQAS(Phospho)PK_2 Mixture1              1
#> 223                     SPVPSQQAS(Phospho)PK_3 Mixture1              1
#> 224                     SPVPSQQAS(Phospho)PK_3 Mixture1              1
#> 225                     SPVPSQQAS(Phospho)PK_3 Mixture1              1
#> 226                     SPVPSQQAS(Phospho)PK_3 Mixture1              1
#> 227                     SPVPSQQAS(Phospho)PK_3 Mixture1              1
#> 228                     SPVPSQQAS(Phospho)PK_3 Mixture1              1
#> 229                        SRT(Phospho)TDTFR_3 Mixture1              1
#> 230                        SRT(Phospho)TDTFR_3 Mixture1              1
#> 231                        SRT(Phospho)TDTFR_3 Mixture1              1
#> 232                        SRT(Phospho)TDTFR_3 Mixture1              1
#> 233                        SRT(Phospho)TDTFR_3 Mixture1              1
#> 234                        SRT(Phospho)TDTFR_3 Mixture1              1
#> 235         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK_3 Mixture1              1
#> 236         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK_3 Mixture1              1
#> 237         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK_3 Mixture1              1
#> 238         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK_3 Mixture1              1
#> 239         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK_3 Mixture1              1
#> 240         SSMMPGFFSDVNPLSAVSS(Phospho)VVNK_3 Mixture1              1
#> 241                       SVS(Phospho)SQQAEK_2 Mixture1              1
#> 242                       SVS(Phospho)SQQAEK_2 Mixture1              1
#> 243                       SVS(Phospho)SQQAEK_2 Mixture1              1
#> 244                       SVS(Phospho)SQQAEK_2 Mixture1              1
#> 245                       SVS(Phospho)SQQAEK_2 Mixture1              1
#> 246                       SVS(Phospho)SQQAEK_2 Mixture1              1
#> 247                       SVSS(Phospho)QQAEK_2 Mixture1              1
#> 248                       SVSS(Phospho)QQAEK_2 Mixture1              1
#> 249                       SVSS(Phospho)QQAEK_2 Mixture1              1
#> 250                       SVSS(Phospho)QQAEK_2 Mixture1              1
#> 251                       SVSS(Phospho)QQAEK_2 Mixture1              1
#> 252                       SVSS(Phospho)QQAEK_2 Mixture1              1
#>            Run Channel BioReplicate Condition Intensity
#> 1   Mixture1_1  Int126           S1     Young     35677
#> 2   Mixture1_1  Int127           S2     Young     39390
#> 3   Mixture1_1  Int128           S3     Young     31921
#> 4   Mixture1_1  Int129           S4       Old     42021
#> 5   Mixture1_1  Int130           S5       Old     23273
#> 6   Mixture1_1  Int131           S6       Old     34693
#> 7   Mixture1_1  Int126           S1     Young     32103
#> 8   Mixture1_1  Int127           S2     Young     45075
#> 9   Mixture1_1  Int128           S3     Young     65793
#> 10  Mixture1_1  Int129           S4       Old     53048
#> 11  Mixture1_1  Int130           S5       Old     43588
#> 12  Mixture1_1  Int131           S6       Old     49907
#> 13  Mixture1_1  Int126           S1     Young     11003
#> 14  Mixture1_1  Int127           S2     Young     24600
#> 15  Mixture1_1  Int128           S3     Young     23359
#> 16  Mixture1_1  Int129           S4       Old     29631
#> 17  Mixture1_1  Int130           S5       Old     32580
#> 18  Mixture1_1  Int131           S6       Old     26484
#> 19  Mixture1_1  Int126           S1     Young     21352
#> 20  Mixture1_1  Int127           S2     Young     20708
#> 21  Mixture1_1  Int128           S3     Young     24370
#> 22  Mixture1_1  Int129           S4       Old     27013
#> 23  Mixture1_1  Int130           S5       Old     21308
#> 24  Mixture1_1  Int131           S6       Old     26018
#> 25  Mixture1_1  Int126           S1     Young     39471
#> 26  Mixture1_1  Int127           S2     Young     48551
#> 27  Mixture1_1  Int128           S3     Young     62397
#> 28  Mixture1_1  Int129           S4       Old     54429
#> 29  Mixture1_1  Int130           S5       Old     35338
#> 30  Mixture1_1  Int131           S6       Old     49178
#> 31  Mixture1_1  Int126           S1     Young     10241
#> 32  Mixture1_1  Int127           S2     Young     16523
#> 33  Mixture1_1  Int128           S3     Young     17101
#> 34  Mixture1_1  Int129           S4       Old     16818
#> 35  Mixture1_1  Int130           S5       Old     12515
#> 36  Mixture1_1  Int131           S6       Old     14038
#> 37  Mixture1_1  Int126           S1     Young      5112
#> 38  Mixture1_1  Int127           S2     Young      6876
#> 39  Mixture1_1  Int128           S3     Young      6456
#> 40  Mixture1_1  Int129           S4       Old      9866
#> 41  Mixture1_1  Int130           S5       Old      5502
#> 42  Mixture1_1  Int131           S6       Old      7183
#> 43  Mixture1_1  Int126           S1     Young    140924
#> 44  Mixture1_1  Int127           S2     Young    147982
#> 45  Mixture1_1  Int128           S3     Young    219516
#> 46  Mixture1_1  Int129           S4       Old    187403
#> 47  Mixture1_1  Int130           S5       Old    175154
#> 48  Mixture1_1  Int131           S6       Old    175012
#> 49  Mixture1_1  Int126           S1     Young     14510
#> 50  Mixture1_1  Int127           S2     Young     20882
#> 51  Mixture1_1  Int128           S3     Young     22037
#> 52  Mixture1_1  Int129           S4       Old     20808
#> 53  Mixture1_1  Int130           S5       Old     15400
#> 54  Mixture1_1  Int131           S6       Old     19311
#> 55  Mixture1_1  Int126           S1     Young     43809
#> 56  Mixture1_1  Int127           S2     Young     58956
#> 57  Mixture1_1  Int128           S3     Young     76681
#> 58  Mixture1_1  Int129           S4       Old     66692
#> 59  Mixture1_1  Int130           S5       Old     82919
#> 60  Mixture1_1  Int131           S6       Old     77978
#> 61  Mixture1_1  Int126           S1     Young    107793
#> 62  Mixture1_1  Int127           S2     Young    107894
#> 63  Mixture1_1  Int128           S3     Young    136424
#> 64  Mixture1_1  Int129           S4       Old    137355
#> 65  Mixture1_1  Int130           S5       Old    133689
#> 66  Mixture1_1  Int131           S6       Old    131555
#> 67  Mixture1_1  Int126           S1     Young     72649
#> 68  Mixture1_1  Int127           S2     Young     69923
#> 69  Mixture1_1  Int128           S3     Young     82825
#> 70  Mixture1_1  Int129           S4       Old     96140
#> 71  Mixture1_1  Int130           S5       Old     62012
#> 72  Mixture1_1  Int131           S6       Old     91102
#> 73  Mixture1_1  Int126           S1     Young     31963
#> 74  Mixture1_1  Int127           S2     Young     30481
#> 75  Mixture1_1  Int128           S3     Young     39043
#> 76  Mixture1_1  Int129           S4       Old     36968
#> 77  Mixture1_1  Int130           S5       Old     26465
#> 78  Mixture1_1  Int131           S6       Old     33886
#> 79  Mixture1_1  Int126           S1     Young     30719
#> 80  Mixture1_1  Int127           S2     Young     32962
#> 81  Mixture1_1  Int128           S3     Young     26548
#> 82  Mixture1_1  Int129           S4       Old     35883
#> 83  Mixture1_1  Int130           S5       Old     25994
#> 84  Mixture1_1  Int131           S6       Old     31855
#> 85  Mixture1_1  Int126           S1     Young     11785
#> 86  Mixture1_1  Int127           S2     Young     15112
#> 87  Mixture1_1  Int128           S3     Young     15933
#> 88  Mixture1_1  Int129           S4       Old     13880
#> 89  Mixture1_1  Int130           S5       Old     11916
#> 90  Mixture1_1  Int131           S6       Old     13206
#> 91  Mixture1_1  Int126           S1     Young      5357
#> 92  Mixture1_1  Int127           S2     Young      9176
#> 93  Mixture1_1  Int128           S3     Young     10395
#> 94  Mixture1_1  Int129           S4       Old      8729
#> 95  Mixture1_1  Int130           S5       Old      5227
#> 96  Mixture1_1  Int131           S6       Old     11059
#> 97  Mixture1_1  Int126           S1     Young      5308
#> 98  Mixture1_1  Int127           S2     Young      8466
#> 99  Mixture1_1  Int128           S3     Young      9011
#> 100 Mixture1_1  Int129           S4       Old      9498
#> 101 Mixture1_1  Int130           S5       Old      5131
#> 102 Mixture1_1  Int131           S6       Old     10197
#> 103 Mixture1_1  Int126           S1     Young     72828
#> 104 Mixture1_1  Int127           S2     Young     80934
#> 105 Mixture1_1  Int128           S3     Young     90459
#> 106 Mixture1_1  Int129           S4       Old     78753
#> 107 Mixture1_1  Int130           S5       Old     76561
#> 108 Mixture1_1  Int131           S6       Old     94021
#> 109 Mixture1_1  Int126           S1     Young     65974
#> 110 Mixture1_1  Int127           S2     Young     76727
#> 111 Mixture1_1  Int128           S3     Young     82841
#> 112 Mixture1_1  Int129           S4       Old     80603
#> 113 Mixture1_1  Int130           S5       Old     59127
#> 114 Mixture1_1  Int131           S6       Old     74477
#> 115 Mixture1_1  Int126           S1     Young     37368
#> 116 Mixture1_1  Int127           S2     Young     44992
#> 117 Mixture1_1  Int128           S3     Young     50119
#> 118 Mixture1_1  Int129           S4       Old     49058
#> 119 Mixture1_1  Int130           S5       Old     47124
#> 120 Mixture1_1  Int131           S6       Old     47436
#> 121 Mixture1_1  Int126           S1     Young     25756
#> 122 Mixture1_1  Int127           S2     Young     33050
#> 123 Mixture1_1  Int128           S3     Young     30903
#> 124 Mixture1_1  Int129           S4       Old     32438
#> 125 Mixture1_1  Int130           S5       Old     25555
#> 126 Mixture1_1  Int131           S6       Old     31425
#> 127 Mixture1_1  Int126           S1     Young     44399
#> 128 Mixture1_1  Int127           S2     Young     46672
#> 129 Mixture1_1  Int128           S3     Young     55586
#> 130 Mixture1_1  Int129           S4       Old     52200
#> 131 Mixture1_1  Int130           S5       Old     38990
#> 132 Mixture1_1  Int131           S6       Old     46493
#> 133 Mixture1_1  Int126           S1     Young     36736
#> 134 Mixture1_1  Int127           S2     Young     41448
#> 135 Mixture1_1  Int128           S3     Young     40465
#> 136 Mixture1_1  Int129           S4       Old     45374
#> 137 Mixture1_1  Int130           S5       Old     27285
#> 138 Mixture1_1  Int131           S6       Old     34630
#> 139 Mixture1_1  Int126           S1     Young     17696
#> 140 Mixture1_1  Int127           S2     Young     20819
#> 141 Mixture1_1  Int128           S3     Young     26085
#> 142 Mixture1_1  Int129           S4       Old     27945
#> 143 Mixture1_1  Int130           S5       Old     21373
#> 144 Mixture1_1  Int131           S6       Old     25470
#> 145 Mixture1_1  Int126           S1     Young     69061
#> 146 Mixture1_1  Int127           S2     Young     67346
#> 147 Mixture1_1  Int128           S3     Young    103392
#> 148 Mixture1_1  Int129           S4       Old     86397
#> 149 Mixture1_1  Int130           S5       Old     87316
#> 150 Mixture1_1  Int131           S6       Old     81399
#> 151 Mixture1_1  Int126           S1     Young     67959
#> 152 Mixture1_1  Int127           S2     Young     56744
#> 153 Mixture1_1  Int128           S3     Young     80015
#> 154 Mixture1_1  Int129           S4       Old     64970
#> 155 Mixture1_1  Int130           S5       Old     60874
#> 156 Mixture1_1  Int131           S6       Old     65219
#> 157 Mixture1_1  Int126           S1     Young     59892
#> 158 Mixture1_1  Int127           S2     Young     55796
#> 159 Mixture1_1  Int128           S3     Young     89693
#> 160 Mixture1_1  Int129           S4       Old     62793
#> 161 Mixture1_1  Int130           S5       Old     62563
#> 162 Mixture1_1  Int131           S6       Old     67739
#> 163 Mixture1_1  Int126           S1     Young     37715
#> 164 Mixture1_1  Int127           S2     Young     37284
#> 165 Mixture1_1  Int128           S3     Young     47498
#> 166 Mixture1_1  Int129           S4       Old     38339
#> 167 Mixture1_1  Int130           S5       Old     28866
#> 168 Mixture1_1  Int131           S6       Old     39553
#> 169 Mixture1_1  Int126           S1     Young     13559
#> 170 Mixture1_1  Int127           S2     Young     16787
#> 171 Mixture1_1  Int128           S3     Young     25183
#> 172 Mixture1_1  Int129           S4       Old     16857
#> 173 Mixture1_1  Int130           S5       Old     16661
#> 174 Mixture1_1  Int131           S6       Old     19321
#> 175 Mixture1_1  Int126           S1     Young     79078
#> 176 Mixture1_1  Int127           S2     Young    109026
#> 177 Mixture1_1  Int128           S3     Young    126574
#> 178 Mixture1_1  Int129           S4       Old    112539
#> 179 Mixture1_1  Int130           S5       Old    108233
#> 180 Mixture1_1  Int131           S6       Old    118963
#> 181 Mixture1_1  Int126           S1     Young     27374
#> 182 Mixture1_1  Int127           S2     Young     36066
#> 183 Mixture1_1  Int128           S3     Young     45088
#> 184 Mixture1_1  Int129           S4       Old     43403
#> 185 Mixture1_1  Int130           S5       Old     35079
#> 186 Mixture1_1  Int131           S6       Old     41550
#> 187 Mixture1_1  Int126           S1     Young     74020
#> 188 Mixture1_1  Int127           S2     Young     84756
#> 189 Mixture1_1  Int128           S3     Young     84232
#> 190 Mixture1_1  Int129           S4       Old     81831
#> 191 Mixture1_1  Int130           S5       Old     57146
#> 192 Mixture1_1  Int131           S6       Old     67194
#> 193 Mixture1_1  Int126           S1     Young    233957
#> 194 Mixture1_1  Int127           S2     Young    238535
#> 195 Mixture1_1  Int128           S3     Young    357265
#> 196 Mixture1_1  Int129           S4       Old    268053
#> 197 Mixture1_1  Int130           S5       Old    262127
#> 198 Mixture1_1  Int131           S6       Old    280330
#> 199 Mixture1_1  Int126           S1     Young     18445
#> 200 Mixture1_1  Int127           S2     Young     21808
#> 201 Mixture1_1  Int128           S3     Young     27893
#> 202 Mixture1_1  Int129           S4       Old     24551
#> 203 Mixture1_1  Int130           S5       Old     23956
#> 204 Mixture1_1  Int131           S6       Old     23432
#> 205 Mixture1_1  Int126           S1     Young     53404
#> 206 Mixture1_1  Int127           S2     Young     52476
#> 207 Mixture1_1  Int128           S3     Young     62672
#> 208 Mixture1_1  Int129           S4       Old     57006
#> 209 Mixture1_1  Int130           S5       Old     48006
#> 210 Mixture1_1  Int131           S6       Old     54185
#> 211 Mixture1_1  Int126           S1     Young     24960
#> 212 Mixture1_1  Int127           S2     Young     20458
#> 213 Mixture1_1  Int128           S3     Young     29394
#> 214 Mixture1_1  Int129           S4       Old     27968
#> 215 Mixture1_1  Int130           S5       Old     23712
#> 216 Mixture1_1  Int131           S6       Old     24336
#> 217 Mixture1_1  Int126           S1     Young     73451
#> 218 Mixture1_1  Int127           S2     Young     69715
#> 219 Mixture1_1  Int128           S3     Young     98877
#> 220 Mixture1_1  Int129           S4       Old     79049
#> 221 Mixture1_1  Int130           S5       Old     77695
#> 222 Mixture1_1  Int131           S6       Old     81938
#> 223 Mixture1_1  Int126           S1     Young    110067
#> 224 Mixture1_1  Int127           S2     Young    108569
#> 225 Mixture1_1  Int128           S3     Young    129873
#> 226 Mixture1_1  Int129           S4       Old    108131
#> 227 Mixture1_1  Int130           S5       Old    107446
#> 228 Mixture1_1  Int131           S6       Old    101541
#> 229 Mixture1_1  Int126           S1     Young     28628
#> 230 Mixture1_1  Int127           S2     Young     32874
#> 231 Mixture1_1  Int128           S3     Young     38196
#> 232 Mixture1_1  Int129           S4       Old     35274
#> 233 Mixture1_1  Int130           S5       Old     26482
#> 234 Mixture1_1  Int131           S6       Old     28653
#> 235 Mixture1_1  Int126           S1     Young      5135
#> 236 Mixture1_1  Int127           S2     Young      5607
#> 237 Mixture1_1  Int128           S3     Young      3937
#> 238 Mixture1_1  Int129           S4       Old      4364
#> 239 Mixture1_1  Int130           S5       Old      2594
#> 240 Mixture1_1  Int131           S6       Old      3114
#> 241 Mixture1_1  Int126           S1     Young     24503
#> 242 Mixture1_1  Int127           S2     Young     23789
#> 243 Mixture1_1  Int128           S3     Young     29818
#> 244 Mixture1_1  Int129           S4       Old     24642
#> 245 Mixture1_1  Int130           S5       Old     23918
#> 246 Mixture1_1  Int131           S6       Old     22867
#> 247 Mixture1_1  Int126           S1     Young     32482
#> 248 Mixture1_1  Int127           S2     Young     32183
#> 249 Mixture1_1  Int128           S3     Young     41125
#> 250 Mixture1_1  Int129           S4       Old     35346
#> 251 Mixture1_1  Int130           S5       Old     33958
#> 252 Mixture1_1  Int131           S6       Old     36180
#> 
#> $PROTEIN
#> Object of class "MSstatsValidated"
#>     ProteinName                       PeptideSequence Charge
#> 1        Q9DCN2                    APDAWDYSQGFVNEEMIR      2
#> 2        Q9DCN2                    APDAWDYSQGFVNEEMIR      2
#> 3        Q9DCN2                    APDAWDYSQGFVNEEMIR      2
#> 4        Q9DCN2                    APDAWDYSQGFVNEEMIR      2
#> 5        Q9DCN2                    APDAWDYSQGFVNEEMIR      2
#> 6        Q9DCN2                    APDAWDYSQGFVNEEMIR      2
#> 7        Q9DCN2                    APDAWDYSQGFVNEEMIR      3
#> 8        Q9DCN2                    APDAWDYSQGFVNEEMIR      3
#> 9        Q9DCN2                    APDAWDYSQGFVNEEMIR      3
#> 10       Q9DCN2                    APDAWDYSQGFVNEEMIR      3
#> 11       Q9DCN2                    APDAWDYSQGFVNEEMIR      3
#> 12       Q9DCN2                    APDAWDYSQGFVNEEMIR      3
#> 13       Q60778 DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR      4
#> 14       Q60778 DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR      4
#> 15       Q60778 DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR      4
#> 16       Q60778 DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR      4
#> 17       Q60778 DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR      4
#> 18       Q60778 DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR      4
#> 19       Q6NS52                              DDGQHVWR      3
#> 20       Q6NS52                              DDGQHVWR      3
#> 21       Q6NS52                              DDGQHVWR      3
#> 22       Q6NS52                              DDGQHVWR      3
#> 23       Q6NS52                              DDGQHVWR      3
#> 24       Q6NS52                              DDGQHVWR      3
#> 25       Q6NS52           DHILPPTTICPVVLTMPSAGASVPEER      3
#> 26       Q6NS52           DHILPPTTICPVVLTMPSAGASVPEER      3
#> 27       Q6NS52           DHILPPTTICPVVLTMPSAGASVPEER      3
#> 28       Q6NS52           DHILPPTTICPVVLTMPSAGASVPEER      3
#> 29       Q6NS52           DHILPPTTICPVVLTMPSAGASVPEER      3
#> 30       Q6NS52           DHILPPTTICPVVLTMPSAGASVPEER      3
#> 31       Q9DCN2       DHLPTPGEEPLILMCGPPPMIQFACLPNLER      4
#> 32       Q9DCN2       DHLPTPGEEPLILMCGPPPMIQFACLPNLER      4
#> 33       Q9DCN2       DHLPTPGEEPLILMCGPPPMIQFACLPNLER      4
#> 34       Q9DCN2       DHLPTPGEEPLILMCGPPPMIQFACLPNLER      4
#> 35       Q9DCN2       DHLPTPGEEPLILMCGPPPMIQFACLPNLER      4
#> 36       Q9DCN2       DHLPTPGEEPLILMCGPPPMIQFACLPNLER      4
#> 37       Q6NS52                          DIESSTEIMLDR      2
#> 38       Q6NS52                          DIESSTEIMLDR      2
#> 39       Q6NS52                          DIESSTEIMLDR      2
#> 40       Q6NS52                          DIESSTEIMLDR      2
#> 41       Q6NS52                          DIESSTEIMLDR      2
#> 42       Q6NS52                          DIESSTEIMLDR      2
#> 43       Q6NS52                          DIESSTEIMLDR      3
#> 44       Q6NS52                          DIESSTEIMLDR      3
#> 45       Q6NS52                          DIESSTEIMLDR      3
#> 46       Q6NS52                          DIESSTEIMLDR      3
#> 47       Q6NS52                          DIESSTEIMLDR      3
#> 48       Q6NS52                          DIESSTEIMLDR      3
#> 49       Q9DCN2                          DILLRPELEELR      2
#> 50       Q9DCN2                          DILLRPELEELR      2
#> 51       Q9DCN2                          DILLRPELEELR      2
#> 52       Q9DCN2                          DILLRPELEELR      2
#> 53       Q9DCN2                          DILLRPELEELR      2
#> 54       Q9DCN2                          DILLRPELEELR      2
#> 55       Q9DCN2                          DILLRPELEELR      3
#> 56       Q9DCN2                          DILLRPELEELR      3
#> 57       Q9DCN2                          DILLRPELEELR      3
#> 58       Q9DCN2                          DILLRPELEELR      3
#> 59       Q9DCN2                          DILLRPELEELR      3
#> 60       Q9DCN2                          DILLRPELEELR      3
#> 61       P16546     DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 62       P16546     DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 63       P16546     DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 64       P16546     DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 65       P16546     DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 66       P16546     DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 67       Q6NS52                           DIVCYLSLLER      2
#> 68       Q6NS52                           DIVCYLSLLER      2
#> 69       Q6NS52                           DIVCYLSLLER      2
#> 70       Q6NS52                           DIVCYLSLLER      2
#> 71       Q6NS52                           DIVCYLSLLER      2
#> 72       Q6NS52                           DIVCYLSLLER      2
#> 73       Q6NS52                           DIVCYLSLLER      3
#> 74       Q6NS52                           DIVCYLSLLER      3
#> 75       Q6NS52                           DIVCYLSLLER      3
#> 76       Q6NS52                           DIVCYLSLLER      3
#> 77       Q6NS52                           DIVCYLSLLER      3
#> 78       Q6NS52                           DIVCYLSLLER      3
#> 79       P16546                 DLNSQADSLMTSSAFDTSQVK      3
#> 80       P16546                 DLNSQADSLMTSSAFDTSQVK      3
#> 81       P16546                 DLNSQADSLMTSSAFDTSQVK      3
#> 82       P16546                 DLNSQADSLMTSSAFDTSQVK      3
#> 83       P16546                 DLNSQADSLMTSSAFDTSQVK      3
#> 84       P16546                 DLNSQADSLMTSSAFDTSQVK      3
#> 85       Q6NS52                                DVPDFR      2
#> 86       Q6NS52                                DVPDFR      2
#> 87       Q6NS52                                DVPDFR      2
#> 88       Q6NS52                                DVPDFR      2
#> 89       Q6NS52                                DVPDFR      2
#> 90       Q6NS52                                DVPDFR      2
#> 91       P16546                 EAFLNTEDKGDSLDSVEALIK      3
#> 92       P16546                 EAFLNTEDKGDSLDSVEALIK      3
#> 93       P16546                 EAFLNTEDKGDSLDSVEALIK      3
#> 94       P16546                 EAFLNTEDKGDSLDSVEALIK      3
#> 95       P16546                 EAFLNTEDKGDSLDSVEALIK      3
#> 96       P16546                 EAFLNTEDKGDSLDSVEALIK      3
#> 97       P16546                 EAIVTSEELGQDLEHVEVLQK      3
#> 98       P16546                 EAIVTSEELGQDLEHVEVLQK      3
#> 99       P16546                 EAIVTSEELGQDLEHVEVLQK      3
#> 100      P16546                 EAIVTSEELGQDLEHVEVLQK      3
#> 101      P16546                 EAIVTSEELGQDLEHVEVLQK      3
#> 102      P16546                 EAIVTSEELGQDLEHVEVLQK      3
#> 103      P16546                 EAIVTSEELGQDLEHVEVLQK      4
#> 104      P16546                 EAIVTSEELGQDLEHVEVLQK      4
#> 105      P16546                 EAIVTSEELGQDLEHVEVLQK      4
#> 106      P16546                 EAIVTSEELGQDLEHVEVLQK      4
#> 107      P16546                 EAIVTSEELGQDLEHVEVLQK      4
#> 108      P16546                 EAIVTSEELGQDLEHVEVLQK      4
#> 109      P16546               EKEPIVGSTDYGKDEDSAEALLK      4
#> 110      P16546               EKEPIVGSTDYGKDEDSAEALLK      4
#> 111      P16546               EKEPIVGSTDYGKDEDSAEALLK      4
#> 112      P16546               EKEPIVGSTDYGKDEDSAEALLK      4
#> 113      P16546               EKEPIVGSTDYGKDEDSAEALLK      4
#> 114      P16546               EKEPIVGSTDYGKDEDSAEALLK      4
#> 115      P16546                     EQADYCVSHMKPYVDGK      4
#> 116      P16546                     EQADYCVSHMKPYVDGK      4
#> 117      P16546                     EQADYCVSHMKPYVDGK      4
#> 118      P16546                     EQADYCVSHMKPYVDGK      4
#> 119      P16546                     EQADYCVSHMKPYVDGK      4
#> 120      P16546                     EQADYCVSHMKPYVDGK      4
#> 121      Q9DCN2                              EVISPDTR      2
#> 122      Q9DCN2                              EVISPDTR      2
#> 123      Q9DCN2                              EVISPDTR      2
#> 124      Q9DCN2                              EVISPDTR      2
#> 125      Q9DCN2                              EVISPDTR      2
#> 126      Q9DCN2                              EVISPDTR      2
#> 127      P16546                        FEEFQTDLAAHEER      2
#> 128      P16546                        FEEFQTDLAAHEER      2
#> 129      P16546                        FEEFQTDLAAHEER      2
#> 130      P16546                        FEEFQTDLAAHEER      2
#> 131      P16546                        FEEFQTDLAAHEER      2
#> 132      P16546                        FEEFQTDLAAHEER      2
#> 133      P16546                        FEEFQTDLAAHEER      3
#> 134      P16546                        FEEFQTDLAAHEER      3
#> 135      P16546                        FEEFQTDLAAHEER      3
#> 136      P16546                        FEEFQTDLAAHEER      3
#> 137      P16546                        FEEFQTDLAAHEER      3
#> 138      P16546                        FEEFQTDLAAHEER      3
#> 139      P16546      GEIDAHEDSFKSADESGQALLAASHYASDEVR      5
#> 140      P16546      GEIDAHEDSFKSADESGQALLAASHYASDEVR      5
#> 141      P16546      GEIDAHEDSFKSADESGQALLAASHYASDEVR      5
#> 142      P16546      GEIDAHEDSFKSADESGQALLAASHYASDEVR      5
#> 143      P16546      GEIDAHEDSFKSADESGQALLAASHYASDEVR      5
#> 144      P16546      GEIDAHEDSFKSADESGQALLAASHYASDEVR      5
#> 145      Q9DCN2                              GFVDLVVK      2
#> 146      Q9DCN2                              GFVDLVVK      2
#> 147      Q9DCN2                              GFVDLVVK      2
#> 148      Q9DCN2                              GFVDLVVK      2
#> 149      Q9DCN2                              GFVDLVVK      2
#> 150      Q9DCN2                              GFVDLVVK      2
#> 151      P16546                      GNAMVEEGHFAAEDVK      2
#> 152      P16546                      GNAMVEEGHFAAEDVK      2
#> 153      P16546                      GNAMVEEGHFAAEDVK      2
#> 154      P16546                      GNAMVEEGHFAAEDVK      2
#> 155      P16546                      GNAMVEEGHFAAEDVK      2
#> 156      P16546                      GNAMVEEGHFAAEDVK      2
#> 157      P16546                      GNAMVEEGHFAAEDVK      3
#> 158      P16546                      GNAMVEEGHFAAEDVK      3
#> 159      P16546                      GNAMVEEGHFAAEDVK      3
#> 160      P16546                      GNAMVEEGHFAAEDVK      3
#> 161      P16546                      GNAMVEEGHFAAEDVK      3
#> 162      P16546                      GNAMVEEGHFAAEDVK      3
#> 163      Q6NS52                          GRPEDKLEFMFR      4
#> 164      Q6NS52                          GRPEDKLEFMFR      4
#> 165      Q6NS52                          GRPEDKLEFMFR      4
#> 166      Q6NS52                          GRPEDKLEFMFR      4
#> 167      Q6NS52                          GRPEDKLEFMFR      4
#> 168      Q6NS52                          GRPEDKLEFMFR      4
#> 169      Q6NS52                          GSDKRPTLTDAK      4
#> 170      Q6NS52                          GSDKRPTLTDAK      4
#> 171      Q6NS52                          GSDKRPTLTDAK      4
#> 172      Q6NS52                          GSDKRPTLTDAK      4
#> 173      Q6NS52                          GSDKRPTLTDAK      4
#> 174      Q6NS52                          GSDKRPTLTDAK      4
#> 175      Q6NS52                     HPPVAILPLGTGNDLAR      3
#> 176      Q6NS52                     HPPVAILPLGTGNDLAR      3
#> 177      Q6NS52                     HPPVAILPLGTGNDLAR      3
#> 178      Q6NS52                     HPPVAILPLGTGNDLAR      3
#> 179      Q6NS52                     HPPVAILPLGTGNDLAR      3
#> 180      Q6NS52                     HPPVAILPLGTGNDLAR      3
#> 181      P16546                        HQAFEAELHANADR      3
#> 182      P16546                        HQAFEAELHANADR      3
#> 183      P16546                        HQAFEAELHANADR      3
#> 184      P16546                        HQAFEAELHANADR      3
#> 185      P16546                        HQAFEAELHANADR      3
#> 186      P16546                        HQAFEAELHANADR      3
#> 187      P16546                        HQAFEAELHANADR      4
#> 188      P16546                        HQAFEAELHANADR      4
#> 189      P16546                        HQAFEAELHANADR      4
#> 190      P16546                        HQAFEAELHANADR      4
#> 191      P16546                        HQAFEAELHANADR      4
#> 192      P16546                        HQAFEAELHANADR      4
#> 193      P16546                   IAALQAFADQLIAVDHYAK      3
#> 194      P16546                   IAALQAFADQLIAVDHYAK      3
#> 195      P16546                   IAALQAFADQLIAVDHYAK      3
#> 196      P16546                   IAALQAFADQLIAVDHYAK      3
#> 197      P16546                   IAALQAFADQLIAVDHYAK      3
#> 198      P16546                   IAALQAFADQLIAVDHYAK      3
#> 199      P16546                   IAALQAFADQLIAVDHYAK      4
#> 200      P16546                   IAALQAFADQLIAVDHYAK      4
#> 201      P16546                   IAALQAFADQLIAVDHYAK      4
#> 202      P16546                   IAALQAFADQLIAVDHYAK      4
#> 203      P16546                   IAALQAFADQLIAVDHYAK      4
#> 204      P16546                   IAALQAFADQLIAVDHYAK      4
#> 205      Q9DCN2           IDGNLVIRPYTPVSSDDDKGFVDLVVK      4
#> 206      Q9DCN2           IDGNLVIRPYTPVSSDDDKGFVDLVVK      4
#> 207      Q9DCN2           IDGNLVIRPYTPVSSDDDKGFVDLVVK      4
#> 208      Q9DCN2           IDGNLVIRPYTPVSSDDDKGFVDLVVK      4
#> 209      Q9DCN2           IDGNLVIRPYTPVSSDDDKGFVDLVVK      4
#> 210      Q9DCN2           IDGNLVIRPYTPVSSDDDKGFVDLVVK      4
#> 211      Q9DCN2                   IDGNLVIRPYTPVSSDDDK      3
#> 212      Q9DCN2                   IDGNLVIRPYTPVSSDDDK      3
#> 213      Q9DCN2                   IDGNLVIRPYTPVSSDDDK      3
#> 214      Q9DCN2                   IDGNLVIRPYTPVSSDDDK      3
#> 215      Q9DCN2                   IDGNLVIRPYTPVSSDDDK      3
#> 216      Q9DCN2                   IDGNLVIRPYTPVSSDDDK      3
#> 217      Q9DCN2                              IGDTIEFR      2
#> 218      Q9DCN2                              IGDTIEFR      2
#> 219      Q9DCN2                              IGDTIEFR      2
#> 220      Q9DCN2                              IGDTIEFR      2
#> 221      Q9DCN2                              IGDTIEFR      2
#> 222      Q9DCN2                              IGDTIEFR      2
#> 223      Q6NS52                       ILKDIESSTEIMLDR      3
#> 224      Q6NS52                       ILKDIESSTEIMLDR      3
#> 225      Q6NS52                       ILKDIESSTEIMLDR      3
#> 226      Q6NS52                       ILKDIESSTEIMLDR      3
#> 227      Q6NS52                       ILKDIESSTEIMLDR      3
#> 228      Q6NS52                       ILKDIESSTEIMLDR      3
#> 229      Q9WV30                            IVVQPETQHR      3
#> 230      Q9WV30                            IVVQPETQHR      3
#> 231      Q9WV30                            IVVQPETQHR      3
#> 232      Q9WV30                            IVVQPETQHR      3
#> 233      Q9WV30                            IVVQPETQHR      3
#> 234      Q9WV30                            IVVQPETQHR      3
#> 235      Q6NS52                            LAQCSSVVIR      2
#> 236      Q6NS52                            LAQCSSVVIR      2
#> 237      Q6NS52                            LAQCSSVVIR      2
#> 238      Q6NS52                            LAQCSSVVIR      2
#> 239      Q6NS52                            LAQCSSVVIR      2
#> 240      Q6NS52                            LAQCSSVVIR      2
#> 241      P16546                  LDETGNLMISEGHFASETIR      3
#> 242      P16546                  LDETGNLMISEGHFASETIR      3
#> 243      P16546                  LDETGNLMISEGHFASETIR      3
#> 244      P16546                  LDETGNLMISEGHFASETIR      3
#> 245      P16546                  LDETGNLMISEGHFASETIR      3
#> 246      P16546                  LDETGNLMISEGHFASETIR      3
#> 247      Q6NS52                                LEFMFR      2
#> 248      Q6NS52                                LEFMFR      2
#> 249      Q6NS52                                LEFMFR      2
#> 250      Q6NS52                                LEFMFR      2
#> 251      Q6NS52                                LEFMFR      2
#> 252      Q6NS52                                LEFMFR      2
#> 253      Q9DCN2                         LIDKEVISPDTRR      4
#> 254      Q9DCN2                         LIDKEVISPDTRR      4
#> 255      Q9DCN2                         LIDKEVISPDTRR      4
#> 256      Q9DCN2                         LIDKEVISPDTRR      4
#> 257      Q9DCN2                         LIDKEVISPDTRR      4
#> 258      Q9DCN2                         LIDKEVISPDTRR      4
#> 259      Q9DCN2                          LIDKEVISPDTR      3
#> 260      Q9DCN2                          LIDKEVISPDTR      3
#> 261      Q9DCN2                          LIDKEVISPDTR      3
#> 262      Q9DCN2                          LIDKEVISPDTR      3
#> 263      Q9DCN2                          LIDKEVISPDTR      3
#> 264      Q9DCN2                          LIDKEVISPDTR      3
#> 265      Q9DCN2             LWYTVDKAPDAWDYSQGFVNEEMIR      4
#> 266      Q9DCN2             LWYTVDKAPDAWDYSQGFVNEEMIR      4
#> 267      Q9DCN2             LWYTVDKAPDAWDYSQGFVNEEMIR      4
#> 268      Q9DCN2             LWYTVDKAPDAWDYSQGFVNEEMIR      4
#> 269      Q9DCN2             LWYTVDKAPDAWDYSQGFVNEEMIR      4
#> 270      Q9DCN2             LWYTVDKAPDAWDYSQGFVNEEMIR      4
#> 271      Q9DCN2                               LWYTVDK      2
#> 272      Q9DCN2                               LWYTVDK      2
#> 273      Q9DCN2                               LWYTVDK      2
#> 274      Q9DCN2                               LWYTVDK      2
#> 275      Q9DCN2                               LWYTVDK      2
#> 276      Q9DCN2                               LWYTVDK      2
#> 277      Q9DCN2                             MSQYLENMK      2
#> 278      Q9DCN2                             MSQYLENMK      2
#> 279      Q9DCN2                             MSQYLENMK      2
#> 280      Q9DCN2                             MSQYLENMK      2
#> 281      Q9DCN2                             MSQYLENMK      2
#> 282      Q9DCN2                             MSQYLENMK      2
#> 283      Q9DCN2                             MSQYLENMK      3
#> 284      Q9DCN2                             MSQYLENMK      3
#> 285      Q9DCN2                             MSQYLENMK      3
#> 286      Q9DCN2                             MSQYLENMK      3
#> 287      Q9DCN2                             MSQYLENMK      3
#> 288      Q9DCN2                             MSQYLENMK      3
#> 289      P16546                 NQALNTDNYGHDLASVQALQR      3
#> 290      P16546                 NQALNTDNYGHDLASVQALQR      3
#> 291      P16546                 NQALNTDNYGHDLASVQALQR      3
#> 292      P16546                 NQALNTDNYGHDLASVQALQR      3
#> 293      P16546                 NQALNTDNYGHDLASVQALQR      3
#> 294      P16546                 NQALNTDNYGHDLASVQALQR      3
#> 295      P16546                 NQALNTDNYGHDLASVQALQR      4
#> 296      P16546                 NQALNTDNYGHDLASVQALQR      4
#> 297      P16546                 NQALNTDNYGHDLASVQALQR      4
#> 298      P16546                 NQALNTDNYGHDLASVQALQR      4
#> 299      P16546                 NQALNTDNYGHDLASVQALQR      4
#> 300      P16546                 NQALNTDNYGHDLASVQALQR      4
#> 301      Q6NS52                     NTDVMHHYWVEGNCPTK      4
#> 302      Q6NS52                     NTDVMHHYWVEGNCPTK      4
#> 303      Q6NS52                     NTDVMHHYWVEGNCPTK      4
#> 304      Q6NS52                     NTDVMHHYWVEGNCPTK      4
#> 305      Q6NS52                     NTDVMHHYWVEGNCPTK      4
#> 306      Q6NS52                     NTDVMHHYWVEGNCPTK      4
#> 307      Q6NS52                        QDILNQTIDFEGFK      2
#> 308      Q6NS52                        QDILNQTIDFEGFK      2
#> 309      Q6NS52                        QDILNQTIDFEGFK      2
#> 310      Q6NS52                        QDILNQTIDFEGFK      2
#> 311      Q6NS52                        QDILNQTIDFEGFK      2
#> 312      Q6NS52                        QDILNQTIDFEGFK      2
#> 313      Q6NS52                        QDILNQTIDFEGFK      3
#> 314      Q6NS52                        QDILNQTIDFEGFK      3
#> 315      Q6NS52                        QDILNQTIDFEGFK      3
#> 316      Q6NS52                        QDILNQTIDFEGFK      3
#> 317      Q6NS52                        QDILNQTIDFEGFK      3
#> 318      Q6NS52                        QDILNQTIDFEGFK      3
#> 319      P16546              QDLEDSLQAQQYFADANEAESWMR      3
#> 320      P16546              QDLEDSLQAQQYFADANEAESWMR      3
#> 321      P16546              QDLEDSLQAQQYFADANEAESWMR      3
#> 322      P16546              QDLEDSLQAQQYFADANEAESWMR      3
#> 323      P16546              QDLEDSLQAQQYFADANEAESWMR      3
#> 324      P16546              QDLEDSLQAQQYFADANEAESWMR      3
#> 325      P16546                    QEFAQHANAFHQWIQETR      4
#> 326      P16546                    QEFAQHANAFHQWIQETR      4
#> 327      P16546                    QEFAQHANAFHQWIQETR      4
#> 328      P16546                    QEFAQHANAFHQWIQETR      4
#> 329      P16546                    QEFAQHANAFHQWIQETR      4
#> 330      P16546                    QEFAQHANAFHQWIQETR      4
#> 331      P16546        QETFDAGLQAFQQEGIANITALKDQLLAAK      4
#> 332      P16546        QETFDAGLQAFQQEGIANITALKDQLLAAK      4
#> 333      P16546        QETFDAGLQAFQQEGIANITALKDQLLAAK      4
#> 334      P16546        QETFDAGLQAFQQEGIANITALKDQLLAAK      4
#> 335      P16546        QETFDAGLQAFQQEGIANITALKDQLLAAK      4
#> 336      P16546        QETFDAGLQAFQQEGIANITALKDQLLAAK      4
#> 337      P16546               QETFDAGLQAFQQEGIANITALK      3
#> 338      P16546               QETFDAGLQAFQQEGIANITALK      3
#> 339      P16546               QETFDAGLQAFQQEGIANITALK      3
#> 340      P16546               QETFDAGLQAFQQEGIANITALK      3
#> 341      P16546               QETFDAGLQAFQQEGIANITALK      3
#> 342      P16546               QETFDAGLQAFQQEGIANITALK      3
#> 343      P16546                        QFQDAGHFDAENIK      3
#> 344      P16546                        QFQDAGHFDAENIK      3
#> 345      P16546                        QFQDAGHFDAENIK      3
#> 346      P16546                        QFQDAGHFDAENIK      3
#> 347      P16546                        QFQDAGHFDAENIK      3
#> 348      P16546                        QFQDAGHFDAENIK      3
#> 349      Q6NS52                             QGLCCSFCK      2
#> 350      Q6NS52                             QGLCCSFCK      2
#> 351      Q6NS52                             QGLCCSFCK      2
#> 352      Q6NS52                             QGLCCSFCK      2
#> 353      Q6NS52                             QGLCCSFCK      2
#> 354      Q6NS52                             QGLCCSFCK      2
#> 355      Q6NS52                             QGLCCSFCK      3
#> 356      Q6NS52                             QGLCCSFCK      3
#> 357      Q6NS52                             QGLCCSFCK      3
#> 358      Q6NS52                             QGLCCSFCK      3
#> 359      Q6NS52                             QGLCCSFCK      3
#> 360      Q6NS52                             QGLCCSFCK      3
#> 361      P45591                              QIIVEEAK      2
#> 362      P45591                              QIIVEEAK      2
#> 363      P45591                              QIIVEEAK      2
#> 364      P45591                              QIIVEEAK      2
#> 365      P45591                              QIIVEEAK      2
#> 366      P45591                              QIIVEEAK      2
#> 367      P16546              QQVAPMDDETGKELVLALYDYQEK      3
#> 368      P16546              QQVAPMDDETGKELVLALYDYQEK      3
#> 369      P16546              QQVAPMDDETGKELVLALYDYQEK      3
#> 370      P16546              QQVAPMDDETGKELVLALYDYQEK      3
#> 371      P16546              QQVAPMDDETGKELVLALYDYQEK      3
#> 372      P16546              QQVAPMDDETGKELVLALYDYQEK      3
#> 373      P16546                        QQYEQCMDLQLFYR      2
#> 374      P16546                        QQYEQCMDLQLFYR      2
#> 375      P16546                        QQYEQCMDLQLFYR      2
#> 376      P16546                        QQYEQCMDLQLFYR      2
#> 377      P16546                        QQYEQCMDLQLFYR      2
#> 378      P16546                        QQYEQCMDLQLFYR      2
#> 379      P16546             RQDLEDSLQAQQYFADANEAESWMR      3
#> 380      P16546             RQDLEDSLQAQQYFADANEAESWMR      3
#> 381      P16546             RQDLEDSLQAQQYFADANEAESWMR      3
#> 382      P16546             RQDLEDSLQAQQYFADANEAESWMR      3
#> 383      P16546             RQDLEDSLQAQQYFADANEAESWMR      3
#> 384      P16546             RQDLEDSLQAQQYFADANEAESWMR      3
#> 385      P16546               SADESGQALLAASHYASDEVREK      4
#> 386      P16546               SADESGQALLAASHYASDEVREK      4
#> 387      P16546               SADESGQALLAASHYASDEVREK      4
#> 388      P16546               SADESGQALLAASHYASDEVREK      4
#> 389      P16546               SADESGQALLAASHYASDEVREK      4
#> 390      P16546               SADESGQALLAASHYASDEVREK      4
#> 391      P16546                 SADESGQALLAASHYASDEVR      2
#> 392      P16546                 SADESGQALLAASHYASDEVR      2
#> 393      P16546                 SADESGQALLAASHYASDEVR      2
#> 394      P16546                 SADESGQALLAASHYASDEVR      2
#> 395      P16546                 SADESGQALLAASHYASDEVR      2
#> 396      P16546                 SADESGQALLAASHYASDEVR      2
#> 397      P16546                 SADESGQALLAASHYASDEVR      3
#> 398      P16546                 SADESGQALLAASHYASDEVR      3
#> 399      P16546                 SADESGQALLAASHYASDEVR      3
#> 400      P16546                 SADESGQALLAASHYASDEVR      3
#> 401      P16546                 SADESGQALLAASHYASDEVR      3
#> 402      P16546                 SADESGQALLAASHYASDEVR      3
#> 403      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      3
#> 404      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      3
#> 405      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      3
#> 406      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      3
#> 407      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      3
#> 408      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      3
#> 409      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      4
#> 410      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      4
#> 411      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      4
#> 412      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      4
#> 413      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      4
#> 414      P16546         SLGYDLPMVEEGEPDPEFEAILDTVDPNR      4
#> 415      Q6NS52                   SLPMQIDGEPWMQTPCTIK      3
#> 416      Q6NS52                   SLPMQIDGEPWMQTPCTIK      3
#> 417      Q6NS52                   SLPMQIDGEPWMQTPCTIK      3
#> 418      Q6NS52                   SLPMQIDGEPWMQTPCTIK      3
#> 419      Q6NS52                   SLPMQIDGEPWMQTPCTIK      3
#> 420      Q6NS52                   SLPMQIDGEPWMQTPCTIK      3
#> 421      Q6NS52                   SLPMQIDGEPWMQTPCTIK      4
#> 422      Q6NS52                   SLPMQIDGEPWMQTPCTIK      4
#> 423      Q6NS52                   SLPMQIDGEPWMQTPCTIK      4
#> 424      Q6NS52                   SLPMQIDGEPWMQTPCTIK      4
#> 425      Q6NS52                   SLPMQIDGEPWMQTPCTIK      4
#> 426      Q6NS52                   SLPMQIDGEPWMQTPCTIK      4
#> 427      Q9DCN2                                SNPVVR      2
#> 428      Q9DCN2                                SNPVVR      2
#> 429      Q9DCN2                                SNPVVR      2
#> 430      Q9DCN2                                SNPVVR      2
#> 431      Q9DCN2                                SNPVVR      2
#> 432      Q9DCN2                                SNPVVR      2
#> 433      P16546                    SSLSSAQADFNQLAELDR      2
#> 434      P16546                    SSLSSAQADFNQLAELDR      2
#> 435      P16546                    SSLSSAQADFNQLAELDR      2
#> 436      P16546                    SSLSSAQADFNQLAELDR      2
#> 437      P16546                    SSLSSAQADFNQLAELDR      2
#> 438      P16546                    SSLSSAQADFNQLAELDR      2
#> 439      Q6NS52                       SSPANTCSPEVIHLK      3
#> 440      Q6NS52                       SSPANTCSPEVIHLK      3
#> 441      Q6NS52                       SSPANTCSPEVIHLK      3
#> 442      Q6NS52                       SSPANTCSPEVIHLK      3
#> 443      Q6NS52                       SSPANTCSPEVIHLK      3
#> 444      Q6NS52                       SSPANTCSPEVIHLK      3
#> 445      Q9DCN2                     STPAITLENPDIKYPLR      3
#> 446      Q9DCN2                     STPAITLENPDIKYPLR      3
#> 447      Q9DCN2                     STPAITLENPDIKYPLR      3
#> 448      Q9DCN2                     STPAITLENPDIKYPLR      3
#> 449      Q9DCN2                     STPAITLENPDIKYPLR      3
#> 450      Q9DCN2                     STPAITLENPDIKYPLR      3
#> 451      Q9DCN2                     STPAITLENPDIKYPLR      4
#> 452      Q9DCN2                     STPAITLENPDIKYPLR      4
#> 453      Q9DCN2                     STPAITLENPDIKYPLR      4
#> 454      Q9DCN2                     STPAITLENPDIKYPLR      4
#> 455      Q9DCN2                     STPAITLENPDIKYPLR      4
#> 456      Q9DCN2                     STPAITLENPDIKYPLR      4
#> 457      Q9DCN2                         STPAITLENPDIK      2
#> 458      Q9DCN2                         STPAITLENPDIK      2
#> 459      Q9DCN2                         STPAITLENPDIK      2
#> 460      Q9DCN2                         STPAITLENPDIK      2
#> 461      Q9DCN2                         STPAITLENPDIK      2
#> 462      Q9DCN2                         STPAITLENPDIK      2
#> 463      Q9DCN2                         STPAITLENPDIK      3
#> 464      Q9DCN2                         STPAITLENPDIK      3
#> 465      Q9DCN2                         STPAITLENPDIK      3
#> 466      Q9DCN2                         STPAITLENPDIK      3
#> 467      Q9DCN2                         STPAITLENPDIK      3
#> 468      Q9DCN2                         STPAITLENPDIK      3
#> 469      Q9DCN2                   SVGMIAGGTGITPMLQVIR      3
#> 470      Q9DCN2                   SVGMIAGGTGITPMLQVIR      3
#> 471      Q9DCN2                   SVGMIAGGTGITPMLQVIR      3
#> 472      Q9DCN2                   SVGMIAGGTGITPMLQVIR      3
#> 473      Q9DCN2                   SVGMIAGGTGITPMLQVIR      3
#> 474      Q9DCN2                   SVGMIAGGTGITPMLQVIR      3
#> 475      P16546             TYLLDGSCMVEESGTLESQLEATKR      3
#> 476      P16546             TYLLDGSCMVEESGTLESQLEATKR      3
#> 477      P16546             TYLLDGSCMVEESGTLESQLEATKR      3
#> 478      P16546             TYLLDGSCMVEESGTLESQLEATKR      3
#> 479      P16546             TYLLDGSCMVEESGTLESQLEATKR      3
#> 480      P16546             TYLLDGSCMVEESGTLESQLEATKR      3
#> 481      P16546              TYLLDGSCMVEESGTLESQLEATK      3
#> 482      P16546              TYLLDGSCMVEESGTLESQLEATK      3
#> 483      P16546              TYLLDGSCMVEESGTLESQLEATK      3
#> 484      P16546              TYLLDGSCMVEESGTLESQLEATK      3
#> 485      P16546              TYLLDGSCMVEESGTLESQLEATK      3
#> 486      P16546              TYLLDGSCMVEESGTLESQLEATK      3
#> 487      P16546   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK      4
#> 488      P16546   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK      4
#> 489      P16546   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK      4
#> 490      P16546   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK      4
#> 491      P16546   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK      4
#> 492      P16546   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK      4
#> 493      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      3
#> 494      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      3
#> 495      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      3
#> 496      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      3
#> 497      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      3
#> 498      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      3
#> 499      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 500      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 501      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 502      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 503      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 504      P16546         VAEDLESEGLMAEEVQAVQQQEVYGAMPR      4
#> 505      Q6NS52                          WGGGYEGENLMK      2
#> 506      Q6NS52                          WGGGYEGENLMK      2
#> 507      Q6NS52                          WGGGYEGENLMK      2
#> 508      Q6NS52                          WGGGYEGENLMK      2
#> 509      Q6NS52                          WGGGYEGENLMK      2
#> 510      Q6NS52                          WGGGYEGENLMK      2
#> 511      Q6NS52                          WGGGYEGENLMK      3
#> 512      Q6NS52                          WGGGYEGENLMK      3
#> 513      Q6NS52                          WGGGYEGENLMK      3
#> 514      Q6NS52                          WGGGYEGENLMK      3
#> 515      Q6NS52                          WGGGYEGENLMK      3
#> 516      Q6NS52                          WGGGYEGENLMK      3
#> 517      Q6NS52                               YAEYSTK      2
#> 518      Q6NS52                               YAEYSTK      2
#> 519      Q6NS52                               YAEYSTK      2
#> 520      Q6NS52                               YAEYSTK      2
#> 521      Q6NS52                               YAEYSTK      2
#> 522      Q6NS52                               YAEYSTK      2
#> 523      P16546                YTEHSTVGLAQQWDQLDQLGMR      3
#> 524      P16546                YTEHSTVGLAQQWDQLDQLGMR      3
#> 525      P16546                YTEHSTVGLAQQWDQLDQLGMR      3
#> 526      P16546                YTEHSTVGLAQQWDQLDQLGMR      3
#> 527      P16546                YTEHSTVGLAQQWDQLDQLGMR      3
#> 528      P16546                YTEHSTVGLAQQWDQLDQLGMR      3
#>                                         PSM  Mixture TechRepMixture        Run
#> 1                      APDAWDYSQGFVNEEMIR_2 Mixture1              1 Mixture1_1
#> 2                      APDAWDYSQGFVNEEMIR_2 Mixture1              1 Mixture1_1
#> 3                      APDAWDYSQGFVNEEMIR_2 Mixture1              1 Mixture1_1
#> 4                      APDAWDYSQGFVNEEMIR_2 Mixture1              1 Mixture1_1
#> 5                      APDAWDYSQGFVNEEMIR_2 Mixture1              1 Mixture1_1
#> 6                      APDAWDYSQGFVNEEMIR_2 Mixture1              1 Mixture1_1
#> 7                      APDAWDYSQGFVNEEMIR_3 Mixture1              1 Mixture1_1
#> 8                      APDAWDYSQGFVNEEMIR_3 Mixture1              1 Mixture1_1
#> 9                      APDAWDYSQGFVNEEMIR_3 Mixture1              1 Mixture1_1
#> 10                     APDAWDYSQGFVNEEMIR_3 Mixture1              1 Mixture1_1
#> 11                     APDAWDYSQGFVNEEMIR_3 Mixture1              1 Mixture1_1
#> 12                     APDAWDYSQGFVNEEMIR_3 Mixture1              1 Mixture1_1
#> 13  DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR_4 Mixture1              1 Mixture1_1
#> 14  DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR_4 Mixture1              1 Mixture1_1
#> 15  DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR_4 Mixture1              1 Mixture1_1
#> 16  DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR_4 Mixture1              1 Mixture1_1
#> 17  DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR_4 Mixture1              1 Mixture1_1
#> 18  DASDTYLTQSQDCTPDTSHAPAAVDSQPNPENEEEPR_4 Mixture1              1 Mixture1_1
#> 19                               DDGQHVWR_3 Mixture1              1 Mixture1_1
#> 20                               DDGQHVWR_3 Mixture1              1 Mixture1_1
#> 21                               DDGQHVWR_3 Mixture1              1 Mixture1_1
#> 22                               DDGQHVWR_3 Mixture1              1 Mixture1_1
#> 23                               DDGQHVWR_3 Mixture1              1 Mixture1_1
#> 24                               DDGQHVWR_3 Mixture1              1 Mixture1_1
#> 25            DHILPPTTICPVVLTMPSAGASVPEER_3 Mixture1              1 Mixture1_1
#> 26            DHILPPTTICPVVLTMPSAGASVPEER_3 Mixture1              1 Mixture1_1
#> 27            DHILPPTTICPVVLTMPSAGASVPEER_3 Mixture1              1 Mixture1_1
#> 28            DHILPPTTICPVVLTMPSAGASVPEER_3 Mixture1              1 Mixture1_1
#> 29            DHILPPTTICPVVLTMPSAGASVPEER_3 Mixture1              1 Mixture1_1
#> 30            DHILPPTTICPVVLTMPSAGASVPEER_3 Mixture1              1 Mixture1_1
#> 31        DHLPTPGEEPLILMCGPPPMIQFACLPNLER_4 Mixture1              1 Mixture1_1
#> 32        DHLPTPGEEPLILMCGPPPMIQFACLPNLER_4 Mixture1              1 Mixture1_1
#> 33        DHLPTPGEEPLILMCGPPPMIQFACLPNLER_4 Mixture1              1 Mixture1_1
#> 34        DHLPTPGEEPLILMCGPPPMIQFACLPNLER_4 Mixture1              1 Mixture1_1
#> 35        DHLPTPGEEPLILMCGPPPMIQFACLPNLER_4 Mixture1              1 Mixture1_1
#> 36        DHLPTPGEEPLILMCGPPPMIQFACLPNLER_4 Mixture1              1 Mixture1_1
#> 37                           DIESSTEIMLDR_2 Mixture1              1 Mixture1_1
#> 38                           DIESSTEIMLDR_2 Mixture1              1 Mixture1_1
#> 39                           DIESSTEIMLDR_2 Mixture1              1 Mixture1_1
#> 40                           DIESSTEIMLDR_2 Mixture1              1 Mixture1_1
#> 41                           DIESSTEIMLDR_2 Mixture1              1 Mixture1_1
#> 42                           DIESSTEIMLDR_2 Mixture1              1 Mixture1_1
#> 43                           DIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 44                           DIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 45                           DIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 46                           DIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 47                           DIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 48                           DIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 49                           DILLRPELEELR_2 Mixture1              1 Mixture1_1
#> 50                           DILLRPELEELR_2 Mixture1              1 Mixture1_1
#> 51                           DILLRPELEELR_2 Mixture1              1 Mixture1_1
#> 52                           DILLRPELEELR_2 Mixture1              1 Mixture1_1
#> 53                           DILLRPELEELR_2 Mixture1              1 Mixture1_1
#> 54                           DILLRPELEELR_2 Mixture1              1 Mixture1_1
#> 55                           DILLRPELEELR_3 Mixture1              1 Mixture1_1
#> 56                           DILLRPELEELR_3 Mixture1              1 Mixture1_1
#> 57                           DILLRPELEELR_3 Mixture1              1 Mixture1_1
#> 58                           DILLRPELEELR_3 Mixture1              1 Mixture1_1
#> 59                           DILLRPELEELR_3 Mixture1              1 Mixture1_1
#> 60                           DILLRPELEELR_3 Mixture1              1 Mixture1_1
#> 61      DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 62      DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 63      DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 64      DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 65      DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 66      DINKVAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 67                            DIVCYLSLLER_2 Mixture1              1 Mixture1_1
#> 68                            DIVCYLSLLER_2 Mixture1              1 Mixture1_1
#> 69                            DIVCYLSLLER_2 Mixture1              1 Mixture1_1
#> 70                            DIVCYLSLLER_2 Mixture1              1 Mixture1_1
#> 71                            DIVCYLSLLER_2 Mixture1              1 Mixture1_1
#> 72                            DIVCYLSLLER_2 Mixture1              1 Mixture1_1
#> 73                            DIVCYLSLLER_3 Mixture1              1 Mixture1_1
#> 74                            DIVCYLSLLER_3 Mixture1              1 Mixture1_1
#> 75                            DIVCYLSLLER_3 Mixture1              1 Mixture1_1
#> 76                            DIVCYLSLLER_3 Mixture1              1 Mixture1_1
#> 77                            DIVCYLSLLER_3 Mixture1              1 Mixture1_1
#> 78                            DIVCYLSLLER_3 Mixture1              1 Mixture1_1
#> 79                  DLNSQADSLMTSSAFDTSQVK_3 Mixture1              1 Mixture1_1
#> 80                  DLNSQADSLMTSSAFDTSQVK_3 Mixture1              1 Mixture1_1
#> 81                  DLNSQADSLMTSSAFDTSQVK_3 Mixture1              1 Mixture1_1
#> 82                  DLNSQADSLMTSSAFDTSQVK_3 Mixture1              1 Mixture1_1
#> 83                  DLNSQADSLMTSSAFDTSQVK_3 Mixture1              1 Mixture1_1
#> 84                  DLNSQADSLMTSSAFDTSQVK_3 Mixture1              1 Mixture1_1
#> 85                                 DVPDFR_2 Mixture1              1 Mixture1_1
#> 86                                 DVPDFR_2 Mixture1              1 Mixture1_1
#> 87                                 DVPDFR_2 Mixture1              1 Mixture1_1
#> 88                                 DVPDFR_2 Mixture1              1 Mixture1_1
#> 89                                 DVPDFR_2 Mixture1              1 Mixture1_1
#> 90                                 DVPDFR_2 Mixture1              1 Mixture1_1
#> 91                  EAFLNTEDKGDSLDSVEALIK_3 Mixture1              1 Mixture1_1
#> 92                  EAFLNTEDKGDSLDSVEALIK_3 Mixture1              1 Mixture1_1
#> 93                  EAFLNTEDKGDSLDSVEALIK_3 Mixture1              1 Mixture1_1
#> 94                  EAFLNTEDKGDSLDSVEALIK_3 Mixture1              1 Mixture1_1
#> 95                  EAFLNTEDKGDSLDSVEALIK_3 Mixture1              1 Mixture1_1
#> 96                  EAFLNTEDKGDSLDSVEALIK_3 Mixture1              1 Mixture1_1
#> 97                  EAIVTSEELGQDLEHVEVLQK_3 Mixture1              1 Mixture1_1
#> 98                  EAIVTSEELGQDLEHVEVLQK_3 Mixture1              1 Mixture1_1
#> 99                  EAIVTSEELGQDLEHVEVLQK_3 Mixture1              1 Mixture1_1
#> 100                 EAIVTSEELGQDLEHVEVLQK_3 Mixture1              1 Mixture1_1
#> 101                 EAIVTSEELGQDLEHVEVLQK_3 Mixture1              1 Mixture1_1
#> 102                 EAIVTSEELGQDLEHVEVLQK_3 Mixture1              1 Mixture1_1
#> 103                 EAIVTSEELGQDLEHVEVLQK_4 Mixture1              1 Mixture1_1
#> 104                 EAIVTSEELGQDLEHVEVLQK_4 Mixture1              1 Mixture1_1
#> 105                 EAIVTSEELGQDLEHVEVLQK_4 Mixture1              1 Mixture1_1
#> 106                 EAIVTSEELGQDLEHVEVLQK_4 Mixture1              1 Mixture1_1
#> 107                 EAIVTSEELGQDLEHVEVLQK_4 Mixture1              1 Mixture1_1
#> 108                 EAIVTSEELGQDLEHVEVLQK_4 Mixture1              1 Mixture1_1
#> 109               EKEPIVGSTDYGKDEDSAEALLK_4 Mixture1              1 Mixture1_1
#> 110               EKEPIVGSTDYGKDEDSAEALLK_4 Mixture1              1 Mixture1_1
#> 111               EKEPIVGSTDYGKDEDSAEALLK_4 Mixture1              1 Mixture1_1
#> 112               EKEPIVGSTDYGKDEDSAEALLK_4 Mixture1              1 Mixture1_1
#> 113               EKEPIVGSTDYGKDEDSAEALLK_4 Mixture1              1 Mixture1_1
#> 114               EKEPIVGSTDYGKDEDSAEALLK_4 Mixture1              1 Mixture1_1
#> 115                     EQADYCVSHMKPYVDGK_4 Mixture1              1 Mixture1_1
#> 116                     EQADYCVSHMKPYVDGK_4 Mixture1              1 Mixture1_1
#> 117                     EQADYCVSHMKPYVDGK_4 Mixture1              1 Mixture1_1
#> 118                     EQADYCVSHMKPYVDGK_4 Mixture1              1 Mixture1_1
#> 119                     EQADYCVSHMKPYVDGK_4 Mixture1              1 Mixture1_1
#> 120                     EQADYCVSHMKPYVDGK_4 Mixture1              1 Mixture1_1
#> 121                              EVISPDTR_2 Mixture1              1 Mixture1_1
#> 122                              EVISPDTR_2 Mixture1              1 Mixture1_1
#> 123                              EVISPDTR_2 Mixture1              1 Mixture1_1
#> 124                              EVISPDTR_2 Mixture1              1 Mixture1_1
#> 125                              EVISPDTR_2 Mixture1              1 Mixture1_1
#> 126                              EVISPDTR_2 Mixture1              1 Mixture1_1
#> 127                        FEEFQTDLAAHEER_2 Mixture1              1 Mixture1_1
#> 128                        FEEFQTDLAAHEER_2 Mixture1              1 Mixture1_1
#> 129                        FEEFQTDLAAHEER_2 Mixture1              1 Mixture1_1
#> 130                        FEEFQTDLAAHEER_2 Mixture1              1 Mixture1_1
#> 131                        FEEFQTDLAAHEER_2 Mixture1              1 Mixture1_1
#> 132                        FEEFQTDLAAHEER_2 Mixture1              1 Mixture1_1
#> 133                        FEEFQTDLAAHEER_3 Mixture1              1 Mixture1_1
#> 134                        FEEFQTDLAAHEER_3 Mixture1              1 Mixture1_1
#> 135                        FEEFQTDLAAHEER_3 Mixture1              1 Mixture1_1
#> 136                        FEEFQTDLAAHEER_3 Mixture1              1 Mixture1_1
#> 137                        FEEFQTDLAAHEER_3 Mixture1              1 Mixture1_1
#> 138                        FEEFQTDLAAHEER_3 Mixture1              1 Mixture1_1
#> 139      GEIDAHEDSFKSADESGQALLAASHYASDEVR_5 Mixture1              1 Mixture1_1
#> 140      GEIDAHEDSFKSADESGQALLAASHYASDEVR_5 Mixture1              1 Mixture1_1
#> 141      GEIDAHEDSFKSADESGQALLAASHYASDEVR_5 Mixture1              1 Mixture1_1
#> 142      GEIDAHEDSFKSADESGQALLAASHYASDEVR_5 Mixture1              1 Mixture1_1
#> 143      GEIDAHEDSFKSADESGQALLAASHYASDEVR_5 Mixture1              1 Mixture1_1
#> 144      GEIDAHEDSFKSADESGQALLAASHYASDEVR_5 Mixture1              1 Mixture1_1
#> 145                              GFVDLVVK_2 Mixture1              1 Mixture1_1
#> 146                              GFVDLVVK_2 Mixture1              1 Mixture1_1
#> 147                              GFVDLVVK_2 Mixture1              1 Mixture1_1
#> 148                              GFVDLVVK_2 Mixture1              1 Mixture1_1
#> 149                              GFVDLVVK_2 Mixture1              1 Mixture1_1
#> 150                              GFVDLVVK_2 Mixture1              1 Mixture1_1
#> 151                      GNAMVEEGHFAAEDVK_2 Mixture1              1 Mixture1_1
#> 152                      GNAMVEEGHFAAEDVK_2 Mixture1              1 Mixture1_1
#> 153                      GNAMVEEGHFAAEDVK_2 Mixture1              1 Mixture1_1
#> 154                      GNAMVEEGHFAAEDVK_2 Mixture1              1 Mixture1_1
#> 155                      GNAMVEEGHFAAEDVK_2 Mixture1              1 Mixture1_1
#> 156                      GNAMVEEGHFAAEDVK_2 Mixture1              1 Mixture1_1
#> 157                      GNAMVEEGHFAAEDVK_3 Mixture1              1 Mixture1_1
#> 158                      GNAMVEEGHFAAEDVK_3 Mixture1              1 Mixture1_1
#> 159                      GNAMVEEGHFAAEDVK_3 Mixture1              1 Mixture1_1
#> 160                      GNAMVEEGHFAAEDVK_3 Mixture1              1 Mixture1_1
#> 161                      GNAMVEEGHFAAEDVK_3 Mixture1              1 Mixture1_1
#> 162                      GNAMVEEGHFAAEDVK_3 Mixture1              1 Mixture1_1
#> 163                          GRPEDKLEFMFR_4 Mixture1              1 Mixture1_1
#> 164                          GRPEDKLEFMFR_4 Mixture1              1 Mixture1_1
#> 165                          GRPEDKLEFMFR_4 Mixture1              1 Mixture1_1
#> 166                          GRPEDKLEFMFR_4 Mixture1              1 Mixture1_1
#> 167                          GRPEDKLEFMFR_4 Mixture1              1 Mixture1_1
#> 168                          GRPEDKLEFMFR_4 Mixture1              1 Mixture1_1
#> 169                          GSDKRPTLTDAK_4 Mixture1              1 Mixture1_1
#> 170                          GSDKRPTLTDAK_4 Mixture1              1 Mixture1_1
#> 171                          GSDKRPTLTDAK_4 Mixture1              1 Mixture1_1
#> 172                          GSDKRPTLTDAK_4 Mixture1              1 Mixture1_1
#> 173                          GSDKRPTLTDAK_4 Mixture1              1 Mixture1_1
#> 174                          GSDKRPTLTDAK_4 Mixture1              1 Mixture1_1
#> 175                     HPPVAILPLGTGNDLAR_3 Mixture1              1 Mixture1_1
#> 176                     HPPVAILPLGTGNDLAR_3 Mixture1              1 Mixture1_1
#> 177                     HPPVAILPLGTGNDLAR_3 Mixture1              1 Mixture1_1
#> 178                     HPPVAILPLGTGNDLAR_3 Mixture1              1 Mixture1_1
#> 179                     HPPVAILPLGTGNDLAR_3 Mixture1              1 Mixture1_1
#> 180                     HPPVAILPLGTGNDLAR_3 Mixture1              1 Mixture1_1
#> 181                        HQAFEAELHANADR_3 Mixture1              1 Mixture1_1
#> 182                        HQAFEAELHANADR_3 Mixture1              1 Mixture1_1
#> 183                        HQAFEAELHANADR_3 Mixture1              1 Mixture1_1
#> 184                        HQAFEAELHANADR_3 Mixture1              1 Mixture1_1
#> 185                        HQAFEAELHANADR_3 Mixture1              1 Mixture1_1
#> 186                        HQAFEAELHANADR_3 Mixture1              1 Mixture1_1
#> 187                        HQAFEAELHANADR_4 Mixture1              1 Mixture1_1
#> 188                        HQAFEAELHANADR_4 Mixture1              1 Mixture1_1
#> 189                        HQAFEAELHANADR_4 Mixture1              1 Mixture1_1
#> 190                        HQAFEAELHANADR_4 Mixture1              1 Mixture1_1
#> 191                        HQAFEAELHANADR_4 Mixture1              1 Mixture1_1
#> 192                        HQAFEAELHANADR_4 Mixture1              1 Mixture1_1
#> 193                   IAALQAFADQLIAVDHYAK_3 Mixture1              1 Mixture1_1
#> 194                   IAALQAFADQLIAVDHYAK_3 Mixture1              1 Mixture1_1
#> 195                   IAALQAFADQLIAVDHYAK_3 Mixture1              1 Mixture1_1
#> 196                   IAALQAFADQLIAVDHYAK_3 Mixture1              1 Mixture1_1
#> 197                   IAALQAFADQLIAVDHYAK_3 Mixture1              1 Mixture1_1
#> 198                   IAALQAFADQLIAVDHYAK_3 Mixture1              1 Mixture1_1
#> 199                   IAALQAFADQLIAVDHYAK_4 Mixture1              1 Mixture1_1
#> 200                   IAALQAFADQLIAVDHYAK_4 Mixture1              1 Mixture1_1
#> 201                   IAALQAFADQLIAVDHYAK_4 Mixture1              1 Mixture1_1
#> 202                   IAALQAFADQLIAVDHYAK_4 Mixture1              1 Mixture1_1
#> 203                   IAALQAFADQLIAVDHYAK_4 Mixture1              1 Mixture1_1
#> 204                   IAALQAFADQLIAVDHYAK_4 Mixture1              1 Mixture1_1
#> 205           IDGNLVIRPYTPVSSDDDKGFVDLVVK_4 Mixture1              1 Mixture1_1
#> 206           IDGNLVIRPYTPVSSDDDKGFVDLVVK_4 Mixture1              1 Mixture1_1
#> 207           IDGNLVIRPYTPVSSDDDKGFVDLVVK_4 Mixture1              1 Mixture1_1
#> 208           IDGNLVIRPYTPVSSDDDKGFVDLVVK_4 Mixture1              1 Mixture1_1
#> 209           IDGNLVIRPYTPVSSDDDKGFVDLVVK_4 Mixture1              1 Mixture1_1
#> 210           IDGNLVIRPYTPVSSDDDKGFVDLVVK_4 Mixture1              1 Mixture1_1
#> 211                   IDGNLVIRPYTPVSSDDDK_3 Mixture1              1 Mixture1_1
#> 212                   IDGNLVIRPYTPVSSDDDK_3 Mixture1              1 Mixture1_1
#> 213                   IDGNLVIRPYTPVSSDDDK_3 Mixture1              1 Mixture1_1
#> 214                   IDGNLVIRPYTPVSSDDDK_3 Mixture1              1 Mixture1_1
#> 215                   IDGNLVIRPYTPVSSDDDK_3 Mixture1              1 Mixture1_1
#> 216                   IDGNLVIRPYTPVSSDDDK_3 Mixture1              1 Mixture1_1
#> 217                              IGDTIEFR_2 Mixture1              1 Mixture1_1
#> 218                              IGDTIEFR_2 Mixture1              1 Mixture1_1
#> 219                              IGDTIEFR_2 Mixture1              1 Mixture1_1
#> 220                              IGDTIEFR_2 Mixture1              1 Mixture1_1
#> 221                              IGDTIEFR_2 Mixture1              1 Mixture1_1
#> 222                              IGDTIEFR_2 Mixture1              1 Mixture1_1
#> 223                       ILKDIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 224                       ILKDIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 225                       ILKDIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 226                       ILKDIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 227                       ILKDIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 228                       ILKDIESSTEIMLDR_3 Mixture1              1 Mixture1_1
#> 229                            IVVQPETQHR_3 Mixture1              1 Mixture1_1
#> 230                            IVVQPETQHR_3 Mixture1              1 Mixture1_1
#> 231                            IVVQPETQHR_3 Mixture1              1 Mixture1_1
#> 232                            IVVQPETQHR_3 Mixture1              1 Mixture1_1
#> 233                            IVVQPETQHR_3 Mixture1              1 Mixture1_1
#> 234                            IVVQPETQHR_3 Mixture1              1 Mixture1_1
#> 235                            LAQCSSVVIR_2 Mixture1              1 Mixture1_1
#> 236                            LAQCSSVVIR_2 Mixture1              1 Mixture1_1
#> 237                            LAQCSSVVIR_2 Mixture1              1 Mixture1_1
#> 238                            LAQCSSVVIR_2 Mixture1              1 Mixture1_1
#> 239                            LAQCSSVVIR_2 Mixture1              1 Mixture1_1
#> 240                            LAQCSSVVIR_2 Mixture1              1 Mixture1_1
#> 241                  LDETGNLMISEGHFASETIR_3 Mixture1              1 Mixture1_1
#> 242                  LDETGNLMISEGHFASETIR_3 Mixture1              1 Mixture1_1
#> 243                  LDETGNLMISEGHFASETIR_3 Mixture1              1 Mixture1_1
#> 244                  LDETGNLMISEGHFASETIR_3 Mixture1              1 Mixture1_1
#> 245                  LDETGNLMISEGHFASETIR_3 Mixture1              1 Mixture1_1
#> 246                  LDETGNLMISEGHFASETIR_3 Mixture1              1 Mixture1_1
#> 247                                LEFMFR_2 Mixture1              1 Mixture1_1
#> 248                                LEFMFR_2 Mixture1              1 Mixture1_1
#> 249                                LEFMFR_2 Mixture1              1 Mixture1_1
#> 250                                LEFMFR_2 Mixture1              1 Mixture1_1
#> 251                                LEFMFR_2 Mixture1              1 Mixture1_1
#> 252                                LEFMFR_2 Mixture1              1 Mixture1_1
#> 253                         LIDKEVISPDTRR_4 Mixture1              1 Mixture1_1
#> 254                         LIDKEVISPDTRR_4 Mixture1              1 Mixture1_1
#> 255                         LIDKEVISPDTRR_4 Mixture1              1 Mixture1_1
#> 256                         LIDKEVISPDTRR_4 Mixture1              1 Mixture1_1
#> 257                         LIDKEVISPDTRR_4 Mixture1              1 Mixture1_1
#> 258                         LIDKEVISPDTRR_4 Mixture1              1 Mixture1_1
#> 259                          LIDKEVISPDTR_3 Mixture1              1 Mixture1_1
#> 260                          LIDKEVISPDTR_3 Mixture1              1 Mixture1_1
#> 261                          LIDKEVISPDTR_3 Mixture1              1 Mixture1_1
#> 262                          LIDKEVISPDTR_3 Mixture1              1 Mixture1_1
#> 263                          LIDKEVISPDTR_3 Mixture1              1 Mixture1_1
#> 264                          LIDKEVISPDTR_3 Mixture1              1 Mixture1_1
#> 265             LWYTVDKAPDAWDYSQGFVNEEMIR_4 Mixture1              1 Mixture1_1
#> 266             LWYTVDKAPDAWDYSQGFVNEEMIR_4 Mixture1              1 Mixture1_1
#> 267             LWYTVDKAPDAWDYSQGFVNEEMIR_4 Mixture1              1 Mixture1_1
#> 268             LWYTVDKAPDAWDYSQGFVNEEMIR_4 Mixture1              1 Mixture1_1
#> 269             LWYTVDKAPDAWDYSQGFVNEEMIR_4 Mixture1              1 Mixture1_1
#> 270             LWYTVDKAPDAWDYSQGFVNEEMIR_4 Mixture1              1 Mixture1_1
#> 271                               LWYTVDK_2 Mixture1              1 Mixture1_1
#> 272                               LWYTVDK_2 Mixture1              1 Mixture1_1
#> 273                               LWYTVDK_2 Mixture1              1 Mixture1_1
#> 274                               LWYTVDK_2 Mixture1              1 Mixture1_1
#> 275                               LWYTVDK_2 Mixture1              1 Mixture1_1
#> 276                               LWYTVDK_2 Mixture1              1 Mixture1_1
#> 277                             MSQYLENMK_2 Mixture1              1 Mixture1_1
#> 278                             MSQYLENMK_2 Mixture1              1 Mixture1_1
#> 279                             MSQYLENMK_2 Mixture1              1 Mixture1_1
#> 280                             MSQYLENMK_2 Mixture1              1 Mixture1_1
#> 281                             MSQYLENMK_2 Mixture1              1 Mixture1_1
#> 282                             MSQYLENMK_2 Mixture1              1 Mixture1_1
#> 283                             MSQYLENMK_3 Mixture1              1 Mixture1_1
#> 284                             MSQYLENMK_3 Mixture1              1 Mixture1_1
#> 285                             MSQYLENMK_3 Mixture1              1 Mixture1_1
#> 286                             MSQYLENMK_3 Mixture1              1 Mixture1_1
#> 287                             MSQYLENMK_3 Mixture1              1 Mixture1_1
#> 288                             MSQYLENMK_3 Mixture1              1 Mixture1_1
#> 289                 NQALNTDNYGHDLASVQALQR_3 Mixture1              1 Mixture1_1
#> 290                 NQALNTDNYGHDLASVQALQR_3 Mixture1              1 Mixture1_1
#> 291                 NQALNTDNYGHDLASVQALQR_3 Mixture1              1 Mixture1_1
#> 292                 NQALNTDNYGHDLASVQALQR_3 Mixture1              1 Mixture1_1
#> 293                 NQALNTDNYGHDLASVQALQR_3 Mixture1              1 Mixture1_1
#> 294                 NQALNTDNYGHDLASVQALQR_3 Mixture1              1 Mixture1_1
#> 295                 NQALNTDNYGHDLASVQALQR_4 Mixture1              1 Mixture1_1
#> 296                 NQALNTDNYGHDLASVQALQR_4 Mixture1              1 Mixture1_1
#> 297                 NQALNTDNYGHDLASVQALQR_4 Mixture1              1 Mixture1_1
#> 298                 NQALNTDNYGHDLASVQALQR_4 Mixture1              1 Mixture1_1
#> 299                 NQALNTDNYGHDLASVQALQR_4 Mixture1              1 Mixture1_1
#> 300                 NQALNTDNYGHDLASVQALQR_4 Mixture1              1 Mixture1_1
#> 301                     NTDVMHHYWVEGNCPTK_4 Mixture1              1 Mixture1_1
#> 302                     NTDVMHHYWVEGNCPTK_4 Mixture1              1 Mixture1_1
#> 303                     NTDVMHHYWVEGNCPTK_4 Mixture1              1 Mixture1_1
#> 304                     NTDVMHHYWVEGNCPTK_4 Mixture1              1 Mixture1_1
#> 305                     NTDVMHHYWVEGNCPTK_4 Mixture1              1 Mixture1_1
#> 306                     NTDVMHHYWVEGNCPTK_4 Mixture1              1 Mixture1_1
#> 307                        QDILNQTIDFEGFK_2 Mixture1              1 Mixture1_1
#> 308                        QDILNQTIDFEGFK_2 Mixture1              1 Mixture1_1
#> 309                        QDILNQTIDFEGFK_2 Mixture1              1 Mixture1_1
#> 310                        QDILNQTIDFEGFK_2 Mixture1              1 Mixture1_1
#> 311                        QDILNQTIDFEGFK_2 Mixture1              1 Mixture1_1
#> 312                        QDILNQTIDFEGFK_2 Mixture1              1 Mixture1_1
#> 313                        QDILNQTIDFEGFK_3 Mixture1              1 Mixture1_1
#> 314                        QDILNQTIDFEGFK_3 Mixture1              1 Mixture1_1
#> 315                        QDILNQTIDFEGFK_3 Mixture1              1 Mixture1_1
#> 316                        QDILNQTIDFEGFK_3 Mixture1              1 Mixture1_1
#> 317                        QDILNQTIDFEGFK_3 Mixture1              1 Mixture1_1
#> 318                        QDILNQTIDFEGFK_3 Mixture1              1 Mixture1_1
#> 319              QDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 320              QDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 321              QDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 322              QDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 323              QDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 324              QDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 325                    QEFAQHANAFHQWIQETR_4 Mixture1              1 Mixture1_1
#> 326                    QEFAQHANAFHQWIQETR_4 Mixture1              1 Mixture1_1
#> 327                    QEFAQHANAFHQWIQETR_4 Mixture1              1 Mixture1_1
#> 328                    QEFAQHANAFHQWIQETR_4 Mixture1              1 Mixture1_1
#> 329                    QEFAQHANAFHQWIQETR_4 Mixture1              1 Mixture1_1
#> 330                    QEFAQHANAFHQWIQETR_4 Mixture1              1 Mixture1_1
#> 331        QETFDAGLQAFQQEGIANITALKDQLLAAK_4 Mixture1              1 Mixture1_1
#> 332        QETFDAGLQAFQQEGIANITALKDQLLAAK_4 Mixture1              1 Mixture1_1
#> 333        QETFDAGLQAFQQEGIANITALKDQLLAAK_4 Mixture1              1 Mixture1_1
#> 334        QETFDAGLQAFQQEGIANITALKDQLLAAK_4 Mixture1              1 Mixture1_1
#> 335        QETFDAGLQAFQQEGIANITALKDQLLAAK_4 Mixture1              1 Mixture1_1
#> 336        QETFDAGLQAFQQEGIANITALKDQLLAAK_4 Mixture1              1 Mixture1_1
#> 337               QETFDAGLQAFQQEGIANITALK_3 Mixture1              1 Mixture1_1
#> 338               QETFDAGLQAFQQEGIANITALK_3 Mixture1              1 Mixture1_1
#> 339               QETFDAGLQAFQQEGIANITALK_3 Mixture1              1 Mixture1_1
#> 340               QETFDAGLQAFQQEGIANITALK_3 Mixture1              1 Mixture1_1
#> 341               QETFDAGLQAFQQEGIANITALK_3 Mixture1              1 Mixture1_1
#> 342               QETFDAGLQAFQQEGIANITALK_3 Mixture1              1 Mixture1_1
#> 343                        QFQDAGHFDAENIK_3 Mixture1              1 Mixture1_1
#> 344                        QFQDAGHFDAENIK_3 Mixture1              1 Mixture1_1
#> 345                        QFQDAGHFDAENIK_3 Mixture1              1 Mixture1_1
#> 346                        QFQDAGHFDAENIK_3 Mixture1              1 Mixture1_1
#> 347                        QFQDAGHFDAENIK_3 Mixture1              1 Mixture1_1
#> 348                        QFQDAGHFDAENIK_3 Mixture1              1 Mixture1_1
#> 349                             QGLCCSFCK_2 Mixture1              1 Mixture1_1
#> 350                             QGLCCSFCK_2 Mixture1              1 Mixture1_1
#> 351                             QGLCCSFCK_2 Mixture1              1 Mixture1_1
#> 352                             QGLCCSFCK_2 Mixture1              1 Mixture1_1
#> 353                             QGLCCSFCK_2 Mixture1              1 Mixture1_1
#> 354                             QGLCCSFCK_2 Mixture1              1 Mixture1_1
#> 355                             QGLCCSFCK_3 Mixture1              1 Mixture1_1
#> 356                             QGLCCSFCK_3 Mixture1              1 Mixture1_1
#> 357                             QGLCCSFCK_3 Mixture1              1 Mixture1_1
#> 358                             QGLCCSFCK_3 Mixture1              1 Mixture1_1
#> 359                             QGLCCSFCK_3 Mixture1              1 Mixture1_1
#> 360                             QGLCCSFCK_3 Mixture1              1 Mixture1_1
#> 361                              QIIVEEAK_2 Mixture1              1 Mixture1_1
#> 362                              QIIVEEAK_2 Mixture1              1 Mixture1_1
#> 363                              QIIVEEAK_2 Mixture1              1 Mixture1_1
#> 364                              QIIVEEAK_2 Mixture1              1 Mixture1_1
#> 365                              QIIVEEAK_2 Mixture1              1 Mixture1_1
#> 366                              QIIVEEAK_2 Mixture1              1 Mixture1_1
#> 367              QQVAPMDDETGKELVLALYDYQEK_3 Mixture1              1 Mixture1_1
#> 368              QQVAPMDDETGKELVLALYDYQEK_3 Mixture1              1 Mixture1_1
#> 369              QQVAPMDDETGKELVLALYDYQEK_3 Mixture1              1 Mixture1_1
#> 370              QQVAPMDDETGKELVLALYDYQEK_3 Mixture1              1 Mixture1_1
#> 371              QQVAPMDDETGKELVLALYDYQEK_3 Mixture1              1 Mixture1_1
#> 372              QQVAPMDDETGKELVLALYDYQEK_3 Mixture1              1 Mixture1_1
#> 373                        QQYEQCMDLQLFYR_2 Mixture1              1 Mixture1_1
#> 374                        QQYEQCMDLQLFYR_2 Mixture1              1 Mixture1_1
#> 375                        QQYEQCMDLQLFYR_2 Mixture1              1 Mixture1_1
#> 376                        QQYEQCMDLQLFYR_2 Mixture1              1 Mixture1_1
#> 377                        QQYEQCMDLQLFYR_2 Mixture1              1 Mixture1_1
#> 378                        QQYEQCMDLQLFYR_2 Mixture1              1 Mixture1_1
#> 379             RQDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 380             RQDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 381             RQDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 382             RQDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 383             RQDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 384             RQDLEDSLQAQQYFADANEAESWMR_3 Mixture1              1 Mixture1_1
#> 385               SADESGQALLAASHYASDEVREK_4 Mixture1              1 Mixture1_1
#> 386               SADESGQALLAASHYASDEVREK_4 Mixture1              1 Mixture1_1
#> 387               SADESGQALLAASHYASDEVREK_4 Mixture1              1 Mixture1_1
#> 388               SADESGQALLAASHYASDEVREK_4 Mixture1              1 Mixture1_1
#> 389               SADESGQALLAASHYASDEVREK_4 Mixture1              1 Mixture1_1
#> 390               SADESGQALLAASHYASDEVREK_4 Mixture1              1 Mixture1_1
#> 391                 SADESGQALLAASHYASDEVR_2 Mixture1              1 Mixture1_1
#> 392                 SADESGQALLAASHYASDEVR_2 Mixture1              1 Mixture1_1
#> 393                 SADESGQALLAASHYASDEVR_2 Mixture1              1 Mixture1_1
#> 394                 SADESGQALLAASHYASDEVR_2 Mixture1              1 Mixture1_1
#> 395                 SADESGQALLAASHYASDEVR_2 Mixture1              1 Mixture1_1
#> 396                 SADESGQALLAASHYASDEVR_2 Mixture1              1 Mixture1_1
#> 397                 SADESGQALLAASHYASDEVR_3 Mixture1              1 Mixture1_1
#> 398                 SADESGQALLAASHYASDEVR_3 Mixture1              1 Mixture1_1
#> 399                 SADESGQALLAASHYASDEVR_3 Mixture1              1 Mixture1_1
#> 400                 SADESGQALLAASHYASDEVR_3 Mixture1              1 Mixture1_1
#> 401                 SADESGQALLAASHYASDEVR_3 Mixture1              1 Mixture1_1
#> 402                 SADESGQALLAASHYASDEVR_3 Mixture1              1 Mixture1_1
#> 403         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_3 Mixture1              1 Mixture1_1
#> 404         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_3 Mixture1              1 Mixture1_1
#> 405         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_3 Mixture1              1 Mixture1_1
#> 406         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_3 Mixture1              1 Mixture1_1
#> 407         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_3 Mixture1              1 Mixture1_1
#> 408         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_3 Mixture1              1 Mixture1_1
#> 409         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_4 Mixture1              1 Mixture1_1
#> 410         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_4 Mixture1              1 Mixture1_1
#> 411         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_4 Mixture1              1 Mixture1_1
#> 412         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_4 Mixture1              1 Mixture1_1
#> 413         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_4 Mixture1              1 Mixture1_1
#> 414         SLGYDLPMVEEGEPDPEFEAILDTVDPNR_4 Mixture1              1 Mixture1_1
#> 415                   SLPMQIDGEPWMQTPCTIK_3 Mixture1              1 Mixture1_1
#> 416                   SLPMQIDGEPWMQTPCTIK_3 Mixture1              1 Mixture1_1
#> 417                   SLPMQIDGEPWMQTPCTIK_3 Mixture1              1 Mixture1_1
#> 418                   SLPMQIDGEPWMQTPCTIK_3 Mixture1              1 Mixture1_1
#> 419                   SLPMQIDGEPWMQTPCTIK_3 Mixture1              1 Mixture1_1
#> 420                   SLPMQIDGEPWMQTPCTIK_3 Mixture1              1 Mixture1_1
#> 421                   SLPMQIDGEPWMQTPCTIK_4 Mixture1              1 Mixture1_1
#> 422                   SLPMQIDGEPWMQTPCTIK_4 Mixture1              1 Mixture1_1
#> 423                   SLPMQIDGEPWMQTPCTIK_4 Mixture1              1 Mixture1_1
#> 424                   SLPMQIDGEPWMQTPCTIK_4 Mixture1              1 Mixture1_1
#> 425                   SLPMQIDGEPWMQTPCTIK_4 Mixture1              1 Mixture1_1
#> 426                   SLPMQIDGEPWMQTPCTIK_4 Mixture1              1 Mixture1_1
#> 427                                SNPVVR_2 Mixture1              1 Mixture1_1
#> 428                                SNPVVR_2 Mixture1              1 Mixture1_1
#> 429                                SNPVVR_2 Mixture1              1 Mixture1_1
#> 430                                SNPVVR_2 Mixture1              1 Mixture1_1
#> 431                                SNPVVR_2 Mixture1              1 Mixture1_1
#> 432                                SNPVVR_2 Mixture1              1 Mixture1_1
#> 433                    SSLSSAQADFNQLAELDR_2 Mixture1              1 Mixture1_1
#> 434                    SSLSSAQADFNQLAELDR_2 Mixture1              1 Mixture1_1
#> 435                    SSLSSAQADFNQLAELDR_2 Mixture1              1 Mixture1_1
#> 436                    SSLSSAQADFNQLAELDR_2 Mixture1              1 Mixture1_1
#> 437                    SSLSSAQADFNQLAELDR_2 Mixture1              1 Mixture1_1
#> 438                    SSLSSAQADFNQLAELDR_2 Mixture1              1 Mixture1_1
#> 439                       SSPANTCSPEVIHLK_3 Mixture1              1 Mixture1_1
#> 440                       SSPANTCSPEVIHLK_3 Mixture1              1 Mixture1_1
#> 441                       SSPANTCSPEVIHLK_3 Mixture1              1 Mixture1_1
#> 442                       SSPANTCSPEVIHLK_3 Mixture1              1 Mixture1_1
#> 443                       SSPANTCSPEVIHLK_3 Mixture1              1 Mixture1_1
#> 444                       SSPANTCSPEVIHLK_3 Mixture1              1 Mixture1_1
#> 445                     STPAITLENPDIKYPLR_3 Mixture1              1 Mixture1_1
#> 446                     STPAITLENPDIKYPLR_3 Mixture1              1 Mixture1_1
#> 447                     STPAITLENPDIKYPLR_3 Mixture1              1 Mixture1_1
#> 448                     STPAITLENPDIKYPLR_3 Mixture1              1 Mixture1_1
#> 449                     STPAITLENPDIKYPLR_3 Mixture1              1 Mixture1_1
#> 450                     STPAITLENPDIKYPLR_3 Mixture1              1 Mixture1_1
#> 451                     STPAITLENPDIKYPLR_4 Mixture1              1 Mixture1_1
#> 452                     STPAITLENPDIKYPLR_4 Mixture1              1 Mixture1_1
#> 453                     STPAITLENPDIKYPLR_4 Mixture1              1 Mixture1_1
#> 454                     STPAITLENPDIKYPLR_4 Mixture1              1 Mixture1_1
#> 455                     STPAITLENPDIKYPLR_4 Mixture1              1 Mixture1_1
#> 456                     STPAITLENPDIKYPLR_4 Mixture1              1 Mixture1_1
#> 457                         STPAITLENPDIK_2 Mixture1              1 Mixture1_1
#> 458                         STPAITLENPDIK_2 Mixture1              1 Mixture1_1
#> 459                         STPAITLENPDIK_2 Mixture1              1 Mixture1_1
#> 460                         STPAITLENPDIK_2 Mixture1              1 Mixture1_1
#> 461                         STPAITLENPDIK_2 Mixture1              1 Mixture1_1
#> 462                         STPAITLENPDIK_2 Mixture1              1 Mixture1_1
#> 463                         STPAITLENPDIK_3 Mixture1              1 Mixture1_1
#> 464                         STPAITLENPDIK_3 Mixture1              1 Mixture1_1
#> 465                         STPAITLENPDIK_3 Mixture1              1 Mixture1_1
#> 466                         STPAITLENPDIK_3 Mixture1              1 Mixture1_1
#> 467                         STPAITLENPDIK_3 Mixture1              1 Mixture1_1
#> 468                         STPAITLENPDIK_3 Mixture1              1 Mixture1_1
#> 469                   SVGMIAGGTGITPMLQVIR_3 Mixture1              1 Mixture1_1
#> 470                   SVGMIAGGTGITPMLQVIR_3 Mixture1              1 Mixture1_1
#> 471                   SVGMIAGGTGITPMLQVIR_3 Mixture1              1 Mixture1_1
#> 472                   SVGMIAGGTGITPMLQVIR_3 Mixture1              1 Mixture1_1
#> 473                   SVGMIAGGTGITPMLQVIR_3 Mixture1              1 Mixture1_1
#> 474                   SVGMIAGGTGITPMLQVIR_3 Mixture1              1 Mixture1_1
#> 475             TYLLDGSCMVEESGTLESQLEATKR_3 Mixture1              1 Mixture1_1
#> 476             TYLLDGSCMVEESGTLESQLEATKR_3 Mixture1              1 Mixture1_1
#> 477             TYLLDGSCMVEESGTLESQLEATKR_3 Mixture1              1 Mixture1_1
#> 478             TYLLDGSCMVEESGTLESQLEATKR_3 Mixture1              1 Mixture1_1
#> 479             TYLLDGSCMVEESGTLESQLEATKR_3 Mixture1              1 Mixture1_1
#> 480             TYLLDGSCMVEESGTLESQLEATKR_3 Mixture1              1 Mixture1_1
#> 481              TYLLDGSCMVEESGTLESQLEATK_3 Mixture1              1 Mixture1_1
#> 482              TYLLDGSCMVEESGTLESQLEATK_3 Mixture1              1 Mixture1_1
#> 483              TYLLDGSCMVEESGTLESQLEATK_3 Mixture1              1 Mixture1_1
#> 484              TYLLDGSCMVEESGTLESQLEATK_3 Mixture1              1 Mixture1_1
#> 485              TYLLDGSCMVEESGTLESQLEATK_3 Mixture1              1 Mixture1_1
#> 486              TYLLDGSCMVEESGTLESQLEATK_3 Mixture1              1 Mixture1_1
#> 487   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK_4 Mixture1              1 Mixture1_1
#> 488   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK_4 Mixture1              1 Mixture1_1
#> 489   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK_4 Mixture1              1 Mixture1_1
#> 490   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK_4 Mixture1              1 Mixture1_1
#> 491   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK_4 Mixture1              1 Mixture1_1
#> 492   VAEDLESEGLMAEEVQAVQQQEVYGAMPRDEADSK_4 Mixture1              1 Mixture1_1
#> 493         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_3 Mixture1              1 Mixture1_1
#> 494         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_3 Mixture1              1 Mixture1_1
#> 495         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_3 Mixture1              1 Mixture1_1
#> 496         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_3 Mixture1              1 Mixture1_1
#> 497         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_3 Mixture1              1 Mixture1_1
#> 498         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_3 Mixture1              1 Mixture1_1
#> 499         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 500         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 501         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 502         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 503         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 504         VAEDLESEGLMAEEVQAVQQQEVYGAMPR_4 Mixture1              1 Mixture1_1
#> 505                          WGGGYEGENLMK_2 Mixture1              1 Mixture1_1
#> 506                          WGGGYEGENLMK_2 Mixture1              1 Mixture1_1
#> 507                          WGGGYEGENLMK_2 Mixture1              1 Mixture1_1
#> 508                          WGGGYEGENLMK_2 Mixture1              1 Mixture1_1
#> 509                          WGGGYEGENLMK_2 Mixture1              1 Mixture1_1
#> 510                          WGGGYEGENLMK_2 Mixture1              1 Mixture1_1
#> 511                          WGGGYEGENLMK_3 Mixture1              1 Mixture1_1
#> 512                          WGGGYEGENLMK_3 Mixture1              1 Mixture1_1
#> 513                          WGGGYEGENLMK_3 Mixture1              1 Mixture1_1
#> 514                          WGGGYEGENLMK_3 Mixture1              1 Mixture1_1
#> 515                          WGGGYEGENLMK_3 Mixture1              1 Mixture1_1
#> 516                          WGGGYEGENLMK_3 Mixture1              1 Mixture1_1
#> 517                               YAEYSTK_2 Mixture1              1 Mixture1_1
#> 518                               YAEYSTK_2 Mixture1              1 Mixture1_1
#> 519                               YAEYSTK_2 Mixture1              1 Mixture1_1
#> 520                               YAEYSTK_2 Mixture1              1 Mixture1_1
#> 521                               YAEYSTK_2 Mixture1              1 Mixture1_1
#> 522                               YAEYSTK_2 Mixture1              1 Mixture1_1
#> 523                YTEHSTVGLAQQWDQLDQLGMR_3 Mixture1              1 Mixture1_1
#> 524                YTEHSTVGLAQQWDQLDQLGMR_3 Mixture1              1 Mixture1_1
#> 525                YTEHSTVGLAQQWDQLDQLGMR_3 Mixture1              1 Mixture1_1
#> 526                YTEHSTVGLAQQWDQLDQLGMR_3 Mixture1              1 Mixture1_1
#> 527                YTEHSTVGLAQQWDQLDQLGMR_3 Mixture1              1 Mixture1_1
#> 528                YTEHSTVGLAQQWDQLDQLGMR_3 Mixture1              1 Mixture1_1
#>     Channel BioReplicate Condition  Intensity
#> 1    Int126           S1     Young    7861.00
#> 2    Int127           S2     Young    4463.00
#> 3    Int128           S3     Young    6274.00
#> 4    Int129           S4      Aged    7005.00
#> 5    Int130           S5      Aged    4565.00
#> 6    Int131           S6      Aged    5753.00
#> 7    Int126           S1     Young    4556.00
#> 8    Int127           S2     Young    4058.00
#> 9    Int128           S3     Young    5983.00
#> 10   Int129           S4      Aged    4357.00
#> 11   Int130           S5      Aged   15606.00
#> 12   Int131           S6      Aged    4101.00
#> 13   Int126           S1     Young    1195.00
#> 14   Int127           S2     Young     -92.62
#> 15   Int128           S3     Young    1751.00
#> 16   Int129           S4      Aged    -149.90
#> 17   Int130           S5      Aged    1792.00
#> 18   Int131           S6      Aged    1841.00
#> 19   Int126           S1     Young  127316.00
#> 20   Int127           S2     Young   87173.00
#> 21   Int128           S3     Young  140839.00
#> 22   Int129           S4      Aged   93095.00
#> 23   Int130           S5      Aged   93943.00
#> 24   Int131           S6      Aged  108463.00
#> 25   Int126           S1     Young    4484.00
#> 26   Int127           S2     Young    3381.00
#> 27   Int128           S3     Young    4392.00
#> 28   Int129           S4      Aged    2651.00
#> 29   Int130           S5      Aged    3069.00
#> 30   Int131           S6      Aged    3295.00
#> 31   Int126           S1     Young    4678.00
#> 32   Int127           S2     Young    5114.00
#> 33   Int128           S3     Young    5142.00
#> 34   Int129           S4      Aged    4154.00
#> 35   Int130           S5      Aged    2985.00
#> 36   Int131           S6      Aged    3170.00
#> 37   Int126           S1     Young   86175.00
#> 38   Int127           S2     Young   72498.00
#> 39   Int128           S3     Young  103458.00
#> 40   Int129           S4      Aged   77897.00
#> 41   Int130           S5      Aged   68510.00
#> 42   Int131           S6      Aged   78627.00
#> 43   Int126           S1     Young    7860.00
#> 44   Int127           S2     Young    8355.00
#> 45   Int128           S3     Young    9418.00
#> 46   Int129           S4      Aged    8931.00
#> 47   Int130           S5      Aged    6691.00
#> 48   Int131           S6      Aged    8421.00
#> 49   Int126           S1     Young   14490.00
#> 50   Int127           S2     Young   12342.00
#> 51   Int128           S3     Young   14130.00
#> 52   Int129           S4      Aged   17157.00
#> 53   Int130           S5      Aged   13235.00
#> 54   Int131           S6      Aged   13156.00
#> 55   Int126           S1     Young  103119.00
#> 56   Int127           S2     Young   82124.00
#> 57   Int128           S3     Young  115322.00
#> 58   Int129           S4      Aged  114894.00
#> 59   Int130           S5      Aged  102774.00
#> 60   Int131           S6      Aged  107063.00
#> 61   Int126           S1     Young         NA
#> 62   Int127           S2     Young    7542.00
#> 63   Int128           S3     Young    8204.00
#> 64   Int129           S4      Aged   10578.00
#> 65   Int130           S5      Aged    9828.00
#> 66   Int131           S6      Aged   11855.00
#> 67   Int126           S1     Young   37063.00
#> 68   Int127           S2     Young   32420.00
#> 69   Int128           S3     Young   50419.00
#> 70   Int129           S4      Aged   34120.00
#> 71   Int130           S5      Aged   39668.00
#> 72   Int131           S6      Aged   42109.00
#> 73   Int126           S1     Young    1431.00
#> 74   Int127           S2     Young    -112.40
#> 75   Int128           S3     Young    1797.00
#> 76   Int129           S4      Aged    2305.00
#> 77   Int130           S5      Aged    -179.30
#> 78   Int131           S6      Aged    2232.00
#> 79   Int126           S1     Young  114763.00
#> 80   Int127           S2     Young   94417.00
#> 81   Int128           S3     Young  123461.00
#> 82   Int129           S4      Aged  112301.00
#> 83   Int130           S5      Aged  102351.00
#> 84   Int131           S6      Aged  106529.00
#> 85   Int126           S1     Young   75058.00
#> 86   Int127           S2     Young   69123.00
#> 87   Int128           S3     Young   91346.00
#> 88   Int129           S4      Aged   75810.00
#> 89   Int130           S5      Aged   66094.00
#> 90   Int131           S6      Aged   80694.00
#> 91   Int126           S1     Young   58543.00
#> 92   Int127           S2     Young   46293.00
#> 93   Int128           S3     Young   62226.00
#> 94   Int129           S4      Aged   51159.00
#> 95   Int130           S5      Aged   55753.00
#> 96   Int131           S6      Aged   53177.00
#> 97   Int126           S1     Young  129452.00
#> 98   Int127           S2     Young   98333.00
#> 99   Int128           S3     Young  145714.00
#> 100  Int129           S4      Aged  104492.00
#> 101  Int130           S5      Aged  105819.00
#> 102  Int131           S6      Aged  105938.00
#> 103  Int126           S1     Young    9310.00
#> 104  Int127           S2     Young   10878.00
#> 105  Int128           S3     Young   10768.00
#> 106  Int129           S4      Aged    9343.00
#> 107  Int130           S5      Aged    8329.00
#> 108  Int131           S6      Aged    8114.00
#> 109  Int126           S1     Young   97907.00
#> 110  Int127           S2     Young   90016.00
#> 111  Int128           S3     Young  117288.00
#> 112  Int129           S4      Aged   93059.00
#> 113  Int130           S5      Aged  109068.00
#> 114  Int131           S6      Aged   98845.00
#> 115  Int126           S1     Young   74487.00
#> 116  Int127           S2     Young   65132.00
#> 117  Int128           S3     Young   91415.00
#> 118  Int129           S4      Aged   73564.00
#> 119  Int130           S5      Aged   81503.00
#> 120  Int131           S6      Aged   68471.00
#> 121  Int126           S1     Young  248198.00
#> 122  Int127           S2     Young  203349.00
#> 123  Int128           S3     Young  263747.00
#> 124  Int129           S4      Aged  245719.00
#> 125  Int130           S5      Aged  195150.00
#> 126  Int131           S6      Aged  214466.00
#> 127  Int126           S1     Young   57905.00
#> 128  Int127           S2     Young   59155.00
#> 129  Int128           S3     Young   80130.00
#> 130  Int129           S4      Aged   67666.00
#> 131  Int130           S5      Aged   67548.00
#> 132  Int131           S6      Aged   65014.00
#> 133  Int126           S1     Young   77053.00
#> 134  Int127           S2     Young   64489.00
#> 135  Int128           S3     Young   94444.00
#> 136  Int129           S4      Aged   78504.00
#> 137  Int130           S5      Aged   91308.00
#> 138  Int131           S6      Aged   76634.00
#> 139  Int126           S1     Young   39221.00
#> 140  Int127           S2     Young   46197.00
#> 141  Int128           S3     Young   41622.00
#> 142  Int129           S4      Aged   68307.00
#> 143  Int130           S5      Aged   68667.00
#> 144  Int131           S6      Aged   78506.00
#> 145  Int126           S1     Young  199082.00
#> 146  Int127           S2     Young  142694.00
#> 147  Int128           S3     Young  189549.00
#> 148  Int129           S4      Aged  177034.00
#> 149  Int130           S5      Aged  142824.00
#> 150  Int131           S6      Aged  146929.00
#> 151  Int126           S1     Young   72686.00
#> 152  Int127           S2     Young   72763.00
#> 153  Int128           S3     Young   97587.00
#> 154  Int129           S4      Aged   83471.00
#> 155  Int130           S5      Aged   83319.00
#> 156  Int131           S6      Aged   87863.00
#> 157  Int126           S1     Young  519086.00
#> 158  Int127           S2     Young  477185.00
#> 159  Int128           S3     Young  583524.00
#> 160  Int129           S4      Aged  431363.00
#> 161  Int130           S5      Aged  504550.00
#> 162  Int131           S6      Aged  449005.00
#> 163  Int126           S1     Young   33401.00
#> 164  Int127           S2     Young   27558.00
#> 165  Int128           S3     Young   39045.00
#> 166  Int129           S4      Aged   28675.00
#> 167  Int130           S5      Aged   33019.00
#> 168  Int131           S6      Aged   31281.00
#> 169  Int126           S1     Young  102499.00
#> 170  Int127           S2     Young   96299.00
#> 171  Int128           S3     Young  127058.00
#> 172  Int129           S4      Aged  105268.00
#> 173  Int130           S5      Aged  110667.00
#> 174  Int131           S6      Aged  119409.00
#> 175  Int126           S1     Young   12665.00
#> 176  Int127           S2     Young   13581.00
#> 177  Int128           S3     Young   16917.00
#> 178  Int129           S4      Aged   10973.00
#> 179  Int130           S5      Aged    7587.00
#> 180  Int131           S6      Aged   18806.00
#> 181  Int126           S1     Young  209946.00
#> 182  Int127           S2     Young  220134.00
#> 183  Int128           S3     Young  283347.00
#> 184  Int129           S4      Aged  240173.00
#> 185  Int130           S5      Aged  256718.00
#> 186  Int131           S6      Aged  263683.00
#> 187  Int126           S1     Young   51803.00
#> 188  Int127           S2     Young   47016.00
#> 189  Int128           S3     Young   46274.00
#> 190  Int129           S4      Aged   51613.00
#> 191  Int130           S5      Aged   44682.00
#> 192  Int131           S6      Aged   46237.00
#> 193  Int126           S1     Young   26874.00
#> 194  Int127           S2     Young   26781.00
#> 195  Int128           S3     Young   34611.00
#> 196  Int129           S4      Aged   25013.00
#> 197  Int130           S5      Aged   34561.00
#> 198  Int131           S6      Aged   31381.00
#> 199  Int126           S1     Young    3546.00
#> 200  Int127           S2     Young    4606.00
#> 201  Int128           S3     Young    4275.00
#> 202  Int129           S4      Aged    4190.00
#> 203  Int130           S5      Aged    2733.00
#> 204  Int131           S6      Aged    4590.00
#> 205  Int126           S1     Young   17042.00
#> 206  Int127           S2     Young   15018.00
#> 207  Int128           S3     Young   19853.00
#> 208  Int129           S4      Aged   18954.00
#> 209  Int130           S5      Aged   19135.00
#> 210  Int131           S6      Aged   20350.00
#> 211  Int126           S1     Young  156331.00
#> 212  Int127           S2     Young  114585.00
#> 213  Int128           S3     Young  144215.00
#> 214  Int129           S4      Aged  136545.00
#> 215  Int130           S5      Aged  125707.00
#> 216  Int131           S6      Aged  114838.00
#> 217  Int126           S1     Young  217046.00
#> 218  Int127           S2     Young  180073.00
#> 219  Int128           S3     Young  218970.00
#> 220  Int129           S4      Aged  219230.00
#> 221  Int130           S5      Aged  189245.00
#> 222  Int131           S6      Aged  208423.00
#> 223  Int126           S1     Young   11645.00
#> 224  Int127           S2     Young   15196.00
#> 225  Int128           S3     Young   19474.00
#> 226  Int129           S4      Aged   15814.00
#> 227  Int130           S5      Aged   19126.00
#> 228  Int131           S6      Aged   21208.00
#> 229  Int126           S1     Young   43518.00
#> 230  Int127           S2     Young   49149.00
#> 231  Int128           S3     Young   63325.00
#> 232  Int129           S4      Aged   59733.00
#> 233  Int130           S5      Aged   62849.00
#> 234  Int131           S6      Aged   66898.00
#> 235  Int126           S1     Young  169230.00
#> 236  Int127           S2     Young  155667.00
#> 237  Int128           S3     Young  195401.00
#> 238  Int129           S4      Aged  181479.00
#> 239  Int130           S5      Aged  163101.00
#> 240  Int131           S6      Aged  199154.00
#> 241  Int126           S1     Young   19753.00
#> 242  Int127           S2     Young   17551.00
#> 243  Int128           S3     Young   37157.00
#> 244  Int129           S4      Aged   63362.00
#> 245  Int130           S5      Aged   11638.00
#> 246  Int131           S6      Aged  102340.00
#> 247  Int126           S1     Young   90595.00
#> 248  Int127           S2     Young   57583.00
#> 249  Int128           S3     Young  100001.00
#> 250  Int129           S4      Aged   62687.00
#> 251  Int130           S5      Aged   57489.00
#> 252  Int131           S6      Aged   65730.00
#> 253  Int126           S1     Young   54026.00
#> 254  Int127           S2     Young   88645.00
#> 255  Int128           S3     Young   84440.00
#> 256  Int129           S4      Aged  115526.00
#> 257  Int130           S5      Aged  115833.00
#> 258  Int131           S6      Aged  121097.00
#> 259  Int126           S1     Young  353271.00
#> 260  Int127           S2     Young  359141.00
#> 261  Int128           S3     Young  410776.00
#> 262  Int129           S4      Aged  441094.00
#> 263  Int130           S5      Aged  456423.00
#> 264  Int131           S6      Aged  454753.00
#> 265  Int126           S1     Young    3951.00
#> 266  Int127           S2     Young    4348.00
#> 267  Int128           S3     Young    5596.00
#> 268  Int129           S4      Aged    4875.00
#> 269  Int130           S5      Aged    4991.00
#> 270  Int131           S6      Aged    5230.00
#> 271  Int126           S1     Young  230262.00
#> 272  Int127           S2     Young  120946.00
#> 273  Int128           S3     Young  202121.00
#> 274  Int129           S4      Aged  140565.00
#> 275  Int130           S5      Aged  126462.00
#> 276  Int131           S6      Aged  131735.00
#> 277  Int126           S1     Young  274369.00
#> 278  Int127           S2     Young  260641.00
#> 279  Int128           S3     Young  333608.00
#> 280  Int129           S4      Aged  292828.00
#> 281  Int130           S5      Aged  283466.00
#> 282  Int131           S6      Aged  292923.00
#> 283  Int126           S1     Young   17422.00
#> 284  Int127           S2     Young   14461.00
#> 285  Int128           S3     Young   19527.00
#> 286  Int129           S4      Aged   20570.00
#> 287  Int130           S5      Aged   18280.00
#> 288  Int131           S6      Aged   20299.00
#> 289  Int126           S1     Young  309173.00
#> 290  Int127           S2     Young  329409.00
#> 291  Int128           S3     Young  332731.00
#> 292  Int129           S4      Aged  350231.00
#> 293  Int130           S5      Aged  307258.00
#> 294  Int131           S6      Aged  305197.00
#> 295  Int126           S1     Young   12807.00
#> 296  Int127           S2     Young   12774.00
#> 297  Int128           S3     Young   14551.00
#> 298  Int129           S4      Aged   15281.00
#> 299  Int130           S5      Aged   12697.00
#> 300  Int131           S6      Aged   16294.00
#> 301  Int126           S1     Young   25715.00
#> 302  Int127           S2     Young   25052.00
#> 303  Int128           S3     Young   35943.00
#> 304  Int129           S4      Aged   24711.00
#> 305  Int130           S5      Aged   32446.00
#> 306  Int131           S6      Aged   31171.00
#> 307  Int126           S1     Young   35830.00
#> 308  Int127           S2     Young   34670.00
#> 309  Int128           S3     Young   49325.00
#> 310  Int129           S4      Aged   43464.00
#> 311  Int130           S5      Aged   45100.00
#> 312  Int131           S6      Aged   49575.00
#> 313  Int126           S1     Young   59681.00
#> 314  Int127           S2     Young   60937.00
#> 315  Int128           S3     Young   81080.00
#> 316  Int129           S4      Aged   53598.00
#> 317  Int130           S5      Aged   64320.00
#> 318  Int131           S6      Aged   68079.00
#> 319  Int126           S1     Young   10653.00
#> 320  Int127           S2     Young    8954.00
#> 321  Int128           S3     Young   14287.00
#> 322  Int129           S4      Aged   11624.00
#> 323  Int130           S5      Aged    7735.00
#> 324  Int131           S6      Aged    9432.00
#> 325  Int126           S1     Young  261339.00
#> 326  Int127           S2     Young  314666.00
#> 327  Int128           S3     Young  287928.00
#> 328  Int129           S4      Aged  344389.00
#> 329  Int130           S5      Aged  284725.00
#> 330  Int131           S6      Aged  341280.00
#> 331  Int126           S1     Young    6164.00
#> 332  Int127           S2     Young    6510.00
#> 333  Int128           S3     Young    9680.00
#> 334  Int129           S4      Aged    8892.00
#> 335  Int130           S5      Aged    8541.00
#> 336  Int131           S6      Aged    8108.00
#> 337  Int126           S1     Young   29667.00
#> 338  Int127           S2     Young   21042.00
#> 339  Int128           S3     Young   28738.00
#> 340  Int129           S4      Aged   25489.00
#> 341  Int130           S5      Aged   24756.00
#> 342  Int131           S6      Aged   20230.00
#> 343  Int126           S1     Young  405388.00
#> 344  Int127           S2     Young  341579.00
#> 345  Int128           S3     Young  467366.00
#> 346  Int129           S4      Aged  378380.00
#> 347  Int130           S5      Aged  380307.00
#> 348  Int131           S6      Aged  369580.00
#> 349  Int126           S1     Young   95483.00
#> 350  Int127           S2     Young   89353.00
#> 351  Int128           S3     Young  128521.00
#> 352  Int129           S4      Aged   98056.00
#> 353  Int130           S5      Aged   85049.00
#> 354  Int131           S6      Aged  104589.00
#> 355  Int126           S1     Young   36685.00
#> 356  Int127           S2     Young   28080.00
#> 357  Int128           S3     Young   41395.00
#> 358  Int129           S4      Aged   31342.00
#> 359  Int130           S5      Aged   25779.00
#> 360  Int131           S6      Aged   29508.00
#> 361  Int126           S1     Young  397149.00
#> 362  Int127           S2     Young  391438.00
#> 363  Int128           S3     Young  457870.00
#> 364  Int129           S4      Aged  394444.00
#> 365  Int130           S5      Aged  392321.00
#> 366  Int131           S6      Aged  391704.00
#> 367  Int126           S1     Young   18222.00
#> 368  Int127           S2     Young   19487.00
#> 369  Int128           S3     Young   28596.00
#> 370  Int129           S4      Aged   24249.00
#> 371  Int130           S5      Aged   25909.00
#> 372  Int131           S6      Aged   24806.00
#> 373  Int126           S1     Young   25394.00
#> 374  Int127           S2     Young   28637.00
#> 375  Int128           S3     Young   35828.00
#> 376  Int129           S4      Aged   35047.00
#> 377  Int130           S5      Aged   29376.00
#> 378  Int131           S6      Aged   35531.00
#> 379  Int126           S1     Young    9186.00
#> 380  Int127           S2     Young    7492.00
#> 381  Int128           S3     Young    9863.00
#> 382  Int129           S4      Aged    8548.00
#> 383  Int130           S5      Aged   10265.00
#> 384  Int131           S6      Aged    9995.00
#> 385  Int126           S1     Young   41637.00
#> 386  Int127           S2     Young   45788.00
#> 387  Int128           S3     Young   41690.00
#> 388  Int129           S4      Aged   52230.00
#> 389  Int130           S5      Aged   49071.00
#> 390  Int131           S6      Aged   55075.00
#> 391  Int126           S1     Young   16573.00
#> 392  Int127           S2     Young   12783.00
#> 393  Int128           S3     Young   21513.00
#> 394  Int129           S4      Aged   16301.00
#> 395  Int130           S5      Aged   19374.00
#> 396  Int131           S6      Aged   20353.00
#> 397  Int126           S1     Young  129168.00
#> 398  Int127           S2     Young  123658.00
#> 399  Int128           S3     Young  167076.00
#> 400  Int129           S4      Aged  128464.00
#> 401  Int130           S5      Aged  139097.00
#> 402  Int131           S6      Aged  129680.00
#> 403  Int126           S1     Young  890653.00
#> 404  Int127           S2     Young  776233.00
#> 405  Int128           S3     Young 1127000.00
#> 406  Int129           S4      Aged  935587.00
#> 407  Int130           S5      Aged  887078.00
#> 408  Int131           S6      Aged  922503.00
#> 409  Int126           S1     Young  154247.00
#> 410  Int127           S2     Young  191712.00
#> 411  Int128           S3     Young  285989.00
#> 412  Int129           S4      Aged  170749.00
#> 413  Int130           S5      Aged  226149.00
#> 414  Int131           S6      Aged  178835.00
#> 415  Int126           S1     Young   30556.00
#> 416  Int127           S2     Young   25329.00
#> 417  Int128           S3     Young   36435.00
#> 418  Int129           S4      Aged   25681.00
#> 419  Int130           S5      Aged   30169.00
#> 420  Int131           S6      Aged   32793.00
#> 421  Int126           S1     Young    6807.00
#> 422  Int127           S2     Young    6477.00
#> 423  Int128           S3     Young    6491.00
#> 424  Int129           S4      Aged   10069.00
#> 425  Int130           S5      Aged    6389.00
#> 426  Int131           S6      Aged    7822.00
#> 427  Int126           S1     Young  503037.00
#> 428  Int127           S2     Young  440529.00
#> 429  Int128           S3     Young  538371.00
#> 430  Int129           S4      Aged  605475.00
#> 431  Int130           S5      Aged  499790.00
#> 432  Int131           S6      Aged  563186.00
#> 433  Int126           S1     Young   12705.00
#> 434  Int127           S2     Young   10385.00
#> 435  Int128           S3     Young   16191.00
#> 436  Int129           S4      Aged   13560.00
#> 437  Int130           S5      Aged   11680.00
#> 438  Int131           S6      Aged   14529.00
#> 439  Int126           S1     Young   65328.00
#> 440  Int127           S2     Young   53024.00
#> 441  Int128           S3     Young   74512.00
#> 442  Int129           S4      Aged   59279.00
#> 443  Int130           S5      Aged   55524.00
#> 444  Int131           S6      Aged   66125.00
#> 445  Int126           S1     Young  110116.00
#> 446  Int127           S2     Young   61744.00
#> 447  Int128           S3     Young  103338.00
#> 448  Int129           S4      Aged   92655.00
#> 449  Int130           S5      Aged   84067.00
#> 450  Int131           S6      Aged   88576.00
#> 451  Int126           S1     Young   30588.00
#> 452  Int127           S2     Young   19488.00
#> 453  Int128           S3     Young   25335.00
#> 454  Int129           S4      Aged   21228.00
#> 455  Int130           S5      Aged   20312.00
#> 456  Int131           S6      Aged   21478.00
#> 457  Int126           S1     Young   27268.00
#> 458  Int127           S2     Young   16975.00
#> 459  Int128           S3     Young   26186.00
#> 460  Int129           S4      Aged   16114.00
#> 461  Int130           S5      Aged   14777.00
#> 462  Int131           S6      Aged   17144.00
#> 463  Int126           S1     Young   45615.00
#> 464  Int127           S2     Young   21317.00
#> 465  Int128           S3     Young   31804.00
#> 466  Int129           S4      Aged   26127.00
#> 467  Int130           S5      Aged   22500.00
#> 468  Int131           S6      Aged   24651.00
#> 469  Int126           S1     Young   21599.00
#> 470  Int127           S2     Young   20871.00
#> 471  Int128           S3     Young   29493.00
#> 472  Int129           S4      Aged   26153.00
#> 473  Int130           S5      Aged   21770.00
#> 474  Int131           S6      Aged   23449.00
#> 475  Int126           S1     Young   48096.00
#> 476  Int127           S2     Young   44829.00
#> 477  Int128           S3     Young   59058.00
#> 478  Int129           S4      Aged   51355.00
#> 479  Int130           S5      Aged   52614.00
#> 480  Int131           S6      Aged   55591.00
#> 481  Int126           S1     Young   16264.00
#> 482  Int127           S2     Young   14004.00
#> 483  Int128           S3     Young   19997.00
#> 484  Int129           S4      Aged   15074.00
#> 485  Int130           S5      Aged   15489.00
#> 486  Int131           S6      Aged   18055.00
#> 487  Int126           S1     Young    9194.00
#> 488  Int127           S2     Young    8008.00
#> 489  Int128           S3     Young   13329.00
#> 490  Int129           S4      Aged   11390.00
#> 491  Int130           S5      Aged   11518.00
#> 492  Int131           S6      Aged   12827.00
#> 493  Int126           S1     Young    7612.00
#> 494  Int127           S2     Young    4088.00
#> 495  Int128           S3     Young    6987.00
#> 496  Int129           S4      Aged    5719.00
#> 497  Int130           S5      Aged    5293.00
#> 498  Int131           S6      Aged    6796.00
#> 499  Int126           S1     Young    4256.00
#> 500  Int127           S2     Young    9113.00
#> 501  Int128           S3     Young    7752.00
#> 502  Int129           S4      Aged    7307.00
#> 503  Int130           S5      Aged    3149.00
#> 504  Int131           S6      Aged    5315.00
#> 505  Int126           S1     Young   15931.00
#> 506  Int127           S2     Young   18498.00
#> 507  Int128           S3     Young   20685.00
#> 508  Int129           S4      Aged   18930.00
#> 509  Int130           S5      Aged   12392.00
#> 510  Int131           S6      Aged   17782.00
#> 511  Int126           S1     Young  123151.00
#> 512  Int127           S2     Young   98105.00
#> 513  Int128           S3     Young  154558.00
#> 514  Int129           S4      Aged  114114.00
#> 515  Int130           S5      Aged  105281.00
#> 516  Int131           S6      Aged  137472.00
#> 517  Int126           S1     Young  201913.00
#> 518  Int127           S2     Young  199496.00
#> 519  Int128           S3     Young  259596.00
#> 520  Int129           S4      Aged  227527.00
#> 521  Int130           S5      Aged  227203.00
#> 522  Int131           S6      Aged  245297.00
#> 523  Int126           S1     Young    5972.00
#> 524  Int127           S2     Young    3731.00
#> 525  Int128           S3     Young    5337.00
#> 526  Int129           S4      Aged    4200.00
#> 527  Int130           S5      Aged    3611.00
#> 528  Int131           S6      Aged    5444.00
#> 
```
