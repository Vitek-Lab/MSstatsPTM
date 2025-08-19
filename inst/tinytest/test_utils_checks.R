# MSstatsPTM:::.checkAnnotation tests
## Test missing Run/Raw.file columns
annotation_no_run <- data.frame(
    Condition = c("A", "B"),
    BioReplicate = c(1, 2)
)

expect_error(
    MSstatsPTM:::.checkAnnotation(annotation_no_run, "LF"),
    pattern = "Run column missing in annotation file"
)

expect_error(
    MSstatsPTM:::.checkAnnotation(annotation_no_run, "TMT"),
    pattern = "Run column missing in annotation file"
)

## Test valid Run column
annotation_with_run <- data.frame(
    Run = c("run1", "run2"),
    Condition = c("A", "B")
)

expect_silent(MSstatsPTM:::.checkAnnotation(annotation_with_run, "LF"))

## Test valid Raw.file column
annotation_with_rawfile <- data.frame(
    Raw.file = c("file1.raw", "file2.raw"),
    Condition = c("A", "B")
)

expect_silent(MSstatsPTM:::.checkAnnotation(annotation_with_rawfile, "LF"))

## Test both Run and Raw.file columns
annotation_with_both <- data.frame(
    Run = c("run1", "run2"),
    Raw.file = c("file1.raw", "file2.raw"),
    Condition = c("A", "B")
)

expect_silent(MSstatsPTM:::.checkAnnotation(annotation_with_both, "LF"))

## Test LF label_type with all valid columns
annotation_lf_valid <- data.frame(
    Run = c("run1", "run2"),
    Raw.file = c("file1.raw", "file2.raw"),
    Condition = c("A", "B"),
    BioReplicate = c(1, 2),
    IsotopeLabelType = c("L", "H"),
    Fraction = c(1, 1)
)

expect_silent(MSstatsPTM:::.checkAnnotation(annotation_lf_valid, "LF"))

## Test LF label_type with subset of valid columns
annotation_lf_subset <- data.frame(
    Run = c("run1", "run2"),
    Condition = c("A", "B"),
    BioReplicate = c(1, 2)
)

expect_silent(MSstatsPTM:::.checkAnnotation(annotation_lf_subset, "LF"))

## Test TMT label_type with all valid columns
annotation_tmt_valid <- data.frame(
    Run = c("run1", "run2"),
    Raw.file = c("file1.raw", "file2.raw"),
    Fraction = c(1, 1),
    TechRepMixture = c(1, 1),
    Channel = c("126", "127"),
    Condition = c("A", "B"),
    Mixture = c(1, 1),
    BioReplicate = c(1, 2)
)

expect_silent(MSstatsPTM:::.checkAnnotation(annotation_tmt_valid, "TMT"))

## Test TMT label_type with subset of valid columns
annotation_tmt_subset <- data.frame(
    Run = c("run1", "run2"),
    Channel = c("126", "127"),
    Condition = c("A", "B")
)

expect_silent(MSstatsPTM:::.checkAnnotation(annotation_tmt_subset, "TMT"))

## Test extra columns for LF
annotation_lf_extra <- data.frame(
    Run = c("run1", "run2"),
    Condition = c("A", "B"),
    ExtraColumn = c("extra1", "extra2")  # This should not be allowed
)

expect_error(
    MSstatsPTM:::.checkAnnotation(annotation_lf_extra, "LF"),
    pattern = "Extra columns included in the annotation file"
)

## Test that error message contains allowed columns for LF
expect_error(
    MSstatsPTM:::.checkAnnotation(annotation_lf_extra, "LF"),
    pattern = "Run, Raw.file, Condition, BioReplicate, IsotopeLabelType, Fraction"
)

## Test extra columns for TMT
annotation_tmt_extra <- data.frame(
    Run = c("run1", "run2"),
    Channel = c("126", "127"),
    Condition = c("A", "B"),
    InvalidColumn = c("invalid1", "invalid2")  # This should not be allowed
)

expect_error(
    MSstatsPTM:::.checkAnnotation(annotation_tmt_extra, "TMT"),
    pattern = "Extra columns included in the annotation file"
)

## Test that error message contains allowed columns for TMT
expect_error(
    MSstatsPTM:::.checkAnnotation(annotation_tmt_extra, "TMT"),
    pattern = "Run, Raw.file, Fraction, TechRepMixture, Channel, Condition, Mixture, BioReplicate"
)

## Test empty data frame
annotation_empty <- data.frame()

expect_error(
    MSstatsPTM:::.checkAnnotation(annotation_empty, "LF"),
    pattern = "Run column missing in annotation file"
)

## Test only Run column
annotation_only_run <- data.frame(
    Run = c("run1", "run2")
)

expect_silent(MSstatsPTM:::.checkAnnotation(annotation_only_run, "LF"))
expect_silent(MSstatsPTM:::.checkAnnotation(annotation_only_run, "TMT"))

## Test only Raw.file column
annotation_only_rawfile <- data.frame(
    Raw.file = c("file1.raw", "file2.raw")
)

expect_silent(MSstatsPTM:::.checkAnnotation(annotation_only_rawfile, "LF"))
expect_silent(MSstatsPTM:::.checkAnnotation(annotation_only_rawfile, "TMT"))

## Test different label_type values
annotation_basic <- data.frame(
    Run = c("run1", "run2"),
    Condition = c("A", "B")
)

expect_silent(MSstatsPTM:::.checkAnnotation(annotation_basic, "LF"))
expect_silent(MSstatsPTM:::.checkAnnotation(annotation_basic, "TMT"))

## Test unknown label_type
expect_error(
    MSstatsPTM:::.checkAnnotation(annotation_basic, "UNKNOWN"),
    pattern = "Labeling type must be either LF or TMT"
)