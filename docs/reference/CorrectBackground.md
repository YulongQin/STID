# Correct background expression in spatial transcriptomics data

This function corrects background expression in spatial transcriptomics
data using specified background samples and features. It calculates
background gene expression statistics and subtracts them from all
samples.

## Usage

``` r
CorrectBackground(
  STID_obj = NULL,
  loop_id = "LoopAllSamp",
  bg_samp_id = NULL,
  bg_features = NULL,
  PosThres_prob = 0.95,
  adjust_UMI = TRUE,
  assay_id = "Spatial",
  layer_id = "counts",
  grp_nm = NULL,
  dir_nm = "M1_CorrectBackground"
)
```

## Arguments

- STID_obj:

  An STID object containing spatial transcriptomics data

- loop_id:

  Must be "LoopAllSamp" (default: "LoopAllSamp")

- bg_samp_id:

  Character vector specifying background sample IDs

- bg_features:

  Character vector specifying background features (genes) for correction

- PosThres_prob:

  Numeric, probability threshold (0-1) for determining positive
  expression (default: 0.95)

- adjust_UMI:

  Logical, whether to adjust correction values by mean UMI (default:
  TRUE)

- assay_id:

  Character, name of the assay to use (default: "Spatial")

- layer_id:

  Character, name of the layer/data slot to use (default: "counts")

- grp_nm:

  Character, group name for output organization (default: NULL, uses
  timestamp)

- dir_nm:

  Character, directory name for output (default: "M1_CorrectBackground")

## Value

Returns the modified STID object with corrected counts in the Spatial
assay

## Examples

``` r
if (FALSE) { # \dontrun{
# Correct background using specified background samples and features
STID_obj <- CorrectBackground(
  STID_obj = STID_object,
  bg_samp_id = c("sample1", "sample2"),
  bg_features = c("gene1", "gene2", "gene3"),
  PosThres_prob = 0.95,
  adjust_UMI = TRUE
)
} # }
```
