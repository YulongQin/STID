# Scan pathogen count thresholds before SpotDetect_Geneset

This function calculates total pathogen counts and detected pathogen
gene numbers for each spot using a specified pathogen gene vector. It
generates pathogen count distribution curves, pathogen gene number
distribution curves, and positive spot number curves under different
count thresholds.

## Usage

``` r
ScanThreshold(
  STID_obj = NULL,
  pathogen_features = NULL,
  loop_id = "LoopAllSamp",
  meta_key = "raw",
  assay_id = "Spatial",
  layer_id = "counts",
  threshold_range = NULL,
  threshold_quantile = 0.99,
  threshold_n = 100,
  detect_turning_point = TRUE,
  turning_min_points = 3,
  detect_curvature = TRUE,
  curvature_min_points = 3,
  curvature_spar = 0.6,
  turning_transform = c("none", "log2", "zscore"),
  grp_nm = NULL,
  dir_nm = "M1_ScanThreshold"
)
```

## Arguments

- STID_obj:

  An STID object containing spatial transcriptomics data

- pathogen_features:

  Character vector of pathogen gene names

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- meta_key:

  Character, metadata key used by GetMetaData (default: "raw")

- assay_id:

  Character, name of the assay to use (default: "Spatial")

- layer_id:

  Character, name of the layer/data slot to use (default: "counts")

- threshold_range:

  Optional numeric vector of thresholds to scan. Default is NULL, in
  which case thresholds are generated automatically from positive-spot
  pathogen counts.

- threshold_quantile:

  Numeric value in (0, 1\]. The corresponding quantile of Pathogen_Count
  \> 0 spots is used as the automatic upper scan limit (default: 0.99).
  If this value exceeds max(threshold_range), the quantile-based
  automatic range replaces threshold_range.

- threshold_n:

  Integer, target number of evenly spaced threshold candidates when a
  quantile-based automatic range is generated (default: 100). The final
  scan is always limited to at most 200 points.

- detect_turning_point:

  Logical, whether to identify and mark the main turning point
  separately in each facet of the three line plots (default: TRUE). The
  turning point separates the curve into an early rapidly changing
  segment and a later slowly changing segment.

- turning_min_points:

  Integer, minimum number of observations required in each side of the
  two-segment regression, including the shared turning point (default:
  3).

- detect_curvature:

  Logical, whether to identify and mark a second turning point using
  normalized smoothing-spline curvature (default: TRUE).

- curvature_min_points:

  Integer, minimum number of observations retained on each side when
  excluding unstable curve endpoints (default: 3).

- curvature_spar:

  Numeric smoothing parameter passed to smooth.spline; larger values
  produce a smoother curve (default: 0.6).

- turning_transform:

  Character transformation used only for internal turning-point and
  curvature calculations. One of "none", "log2", or "zscore" (default:
  "none"). "log2" applies log2(x + 1) to the x-axis expression/threshold
  values before analysis; "zscore" applies (x - mean(x)) / sd(x). The
  selected point is mapped back by its original row index, so plots,
  reported x/y values and vertical-line positions remain in the original
  coordinate system.

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M1_ScanThreshold")

## Value

A list containing:

- valid_features:

  Valid pathogen genes found in STID_obj

- scan_threshold_range:

  Final threshold values used by both turning-point methods

- threshold_scan_info:

  Threshold source, quantile and point-count metadata

- threshold_quantile_value:

  Calculated positive-spot quantile threshold

- threshold_source:

  Whether the final range came from user input or the quantile rule

- spot_count:

  Pathogen counts and detected gene numbers for each spot

- count_distribution:

  Spot distribution for pathogen counts

- gene_distribution:

  Spot distribution for detected pathogen gene numbers

- threshold_stat:

  Positive spot statistics under each threshold

- turning_point:

  Turning points from two-segment linear regression

- curvature_turning_point:

  Turning points from curvature recognition

- turning_point_all:

  Combined turning points from both methods

- turning_transform:

  Transformation used only during turning-point calculations

- plot_count:

  Pathogen count distribution plot

- plot_count_hist:

  Pathogen count histogram with positive-spot quantile markers

- plot_gene:

  Pathogen gene number distribution plot

- plot_threshold:

  Positive spot number threshold scan plot

## Examples

``` r
if (FALSE) { # \dontrun{
threshold_res <- ScanThreshold(
  STID_obj = STID_object,
  pathogen_features = pathogen_gene_vec,
  threshold_range = NULL,
  threshold_quantile = 0.99,
  threshold_n = 100,
  turning_transform = "log2",
  assay_id = "Spatial",
  layer_id = "counts"
)
} # }
```
