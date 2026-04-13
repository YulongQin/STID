# Calculate Niche Aggregation Index

Computes spatial aggregation metrics for cell types within niches using
k-nearest neighbors and graph-based clustering. Calculates aggregation
index and aggregation region ratios.

## Usage

``` r
CalNicheAggIndex(
  STID_obj = NULL,
  samp_mode = "SS",
  loop_id = "LoopAllSamp",
  samp_grp_index = FALSE,
  meta_key = NULL,
  niche_key = NULL,
  group_by = NULL,
  dist_thres = NULL,
  k_neighbors = 8,
  min_agg_size = 2,
  col = COLOR_List[["PALETTE_WHITE_BG"]][-1],
  return_data = FALSE
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- samp_mode:

  Character, sample type - "SS" (single sample) or "MS" (multi-sample)
  (default: "SS")

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- samp_grp_index:

  Logical, whether to group by sample groups in MS mode (default: FALSE)

- meta_key:

  Character, metadata key for MS mode when niche_key is NULL

- niche_key:

  Character, niche key to analyze

- group_by:

  Character, column name for cell type grouping

- dist_thres:

  Numeric, distance threshold for graph construction (default: NULL,
  uses median of 1st neighbor distances)

- k_neighbors:

  Integer, number of nearest neighbors to consider (default: 8)

- min_agg_size:

  Integer, minimum aggregation cluster size (default: 2)

- col:

  Color palette for visualization (default:
  `COLOR_List$PALETTE_WHITE_BG[-1]`)

- return_data:

  Logical, whether to return the results list (default: FALSE)

## Value

If return_data = TRUE, returns a list of results per sample; otherwise
NULL

## Examples

``` r
if (FALSE) { # \dontrun{
# Calculate aggregation index for niche cells
CalNicheAggIndex(
  STID_obj = STID_obj,
  samp_mode  = "SS",
  niche_key = "niche_virulence",
  group_by = "cell_type",
  k_neighbors = 8,
  min_agg_size = 3
)
} # }
```
