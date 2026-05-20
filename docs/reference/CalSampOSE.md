# Calculate Spatial Organization Entropy (OSE) for samples

Performs spatial organization entropy analysis to quantify tissue
spatial heterogeneity based on cell type distribution patterns. This
function implements the OSE algorithm which partitions tissue into
regions and calculates entropy based on local cell type diversity.

## Usage

``` r
CalSampOSE(
  STID_obj = NULL,
  loop_id = NULL,
  meta_key = NULL,
  group_by = NULL,
  OSE_dist_m = 0.2,
  OSE_PCA_nPC = 15,
  OSE_window = NULL,
  OSE_minSpotNum = NULL,
  only_plot = FALSE,
  plot_params = list(p3_size = 2, p4_size = 0.3, p2_merge_r = 15, p3_merge_size = 9),
  col = COLOR_LIST[["PALETTE_WHITE_BG"]],
  return_data = FALSE,
  grp_nm = NULL,
  dir_nm = "M4_CalSampOSE"
)
```

## Arguments

- STID_obj:

  A STID object containing spatial transcriptomics data

- loop_id:

  Character, multi-sample analysis identifier

- meta_key:

  Character, metadata key for retrieving cell data

- group_by:

  Character, column name for cell type grouping

- OSE_dist_m:

  Numeric, distance weight (default: 0.2). Smaller values give more
  weight to expression distance in region segmentation

- OSE_PCA_nPC:

  Integer, number of PCs for PCA dimensionality reduction (default: 15)

- OSE_window:

  Numeric, window size for superpixel generation. If NULL, automatically
  set to max(x_range, y_range)/10

- OSE_minSpotNum:

  Integer, minimum spots per region. If NULL, automatically set to
  nrow(meta_data)/1000

- only_plot:

  Logical, whether to only regenerate plots from existing results
  (default: FALSE)

- plot_params:

  List of plotting parameters:

  - p3_size: Point size for partition plot

  - p4_size: Point size for group plot

  - p2_merge_r: Radius for pie chart in merged plot

  - p3_merge_size: Point size for entropy in merged plot

- col:

  Character vector, color palette for cell types

- return_data:

  Logical, whether to return results as a list (default: FALSE)

- grp_nm:

  Character, group name for output directory (default: NULL, uses
  timestamp)

- dir_nm:

  Character, directory name for output (default: "M4_CalSampOSE")

## Value

If return_data = TRUE, returns a list containing:

- center_df: Data frame of region centers

- data: Cell-level OSE data with region assignments

- entropy_results: Region entropy calculations

- celltype_percent: Cell type proportions per region

- segment_df: Region boundary segments

## Examples

``` r
if (FALSE) { # \dontrun{
results <- CalSampOSE(
  STID_obj = stid_object,
  loop_id = "comparison_1",
  meta_key = "M2_NicheDetect_STS_20240101",
  group_by = "cell_type",
  OSE_window = 50,
  OSE_minSpotNum = 100
)
} # }
```
