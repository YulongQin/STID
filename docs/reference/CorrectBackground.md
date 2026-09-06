# Correct background expression

Correct background expression in spatial transcriptomics data using
control samples, predefined background regions, or interactively
selected regions.

## Usage

``` r
CorrectBackground(
  STID_obj = NULL,
  loop_id = "LoopAllSamp",
  meta_key = "raw",
  bg_features = NULL,
  ctrl_samp_id = NULL,
  bg_region_cell = NULL,
  bg_region_lasso = FALSE,
  group_by = NULL,
  col = NULL,
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

  An `STID` object.

- loop_id:

  Sample loop name. Currently should be `"LoopAllSamp"`.

- meta_key:

  Metadata key used by
  [`GetMetaData()`](https://yulongqin.github.io/STID/reference/GetMetaData.md).
  Default is `"raw"`.

- bg_features:

  Character vector of features to correct.

- ctrl_samp_id:

  Optional character vector of control sample IDs.

- bg_region_cell:

  Optional data frame containing columns `sample_id` and `bg_cell`.

- bg_region_lasso:

  Logical. Whether to select background regions interactively using
  lasso selection.

- group_by:

  Metadata column used to color spots during lasso selection.

- col:

  Optional colors used for lasso visualization.

- PosThres_prob:

  Quantile of nonzero counts used to estimate background expression.
  Must be in `(0, 1]`. Default is `0.95`.

- adjust_UMI:

  Logical. Whether to scale control-sample correction values according
  to sample mean UMI. Default is `TRUE`.

- assay_id:

  Assay name. Default is `"Spatial"`.

- layer_id:

  Assay layer containing counts. Default is `"counts"`.

- grp_nm:

  Output group name. If `NULL`, a name is generated automatically.

- dir_nm:

  Output directory name. Default is `"M1_CorrectBackground"`.

## Value

The input `STID` object with corrected counts stored in the specified
assay and layer.

## Details

Background expression can be estimated from control samples,
sample-specific background regions, or both. When both are provided, the
larger correction value is used. Corrected counts below zero are set to
zero and rounded to integers.

## Examples

``` r
if (FALSE) { # \dontrun{
# Correct background using specified background samples and features
STID_obj <- CorrectBackground(
  STID_obj = STID_object,
  ctrl_samp_id = c("sample1", "sample2"),
  bg_features = c("gene1", "gene2", "gene3"),
  PosThres_prob = 0.95,
  adjust_UMI = TRUE
)
} # }
```
