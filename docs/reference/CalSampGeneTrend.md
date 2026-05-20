# Calculate gene expression trends across time points

Analyzes gene expression patterns across temporal samples using linear
regression or Mfuzz clustering to identify genes with different trend
patterns.

## Usage

``` r
CalSampGeneTrend(
  STID_obj = NULL,
  loop_id = NULL,
  samp_grp_index = FALSE,
  meta_key = NULL,
  niche_key = NULL,
  group_by = NULL,
  assay_id = "Spatial",
  layer_id = "data",
  method = "fitting",
  fitting_paras = list(lm_min_slope = 2, lm_pval_cutoff = 0.2),
  mfuzz_paras = list(cluster_num = 10, min_acore = 0.3),
  gene_list = NULL,
  remove_genes = NULL,
  col = COLOR_LIST[["PALETTE_WHITE_BG"]],
  return_data = FALSE,
  grp_nm = NULL,
  dir_nm = "M4_CalSampGeneTrend"
)
```

## Arguments

- STID_obj:

  A STID object

- loop_id:

  Character, multi-sample analysis identifier

- samp_grp_index:

  Logical, whether to group by sample groups instead of individual
  sample IDs (default: FALSE)

- meta_key:

  Character, metadata key for retrieving cell data

- niche_key:

  Character, niche key for filtering niche cells (optional)

- group_by:

  Character, column name for cell type grouping

- assay_id:

  Character, assay name (default: "Spatial")

- layer_id:

  Character, layer/slot name (default: "data")

- method:

  Character, trend analysis method - "fitting" (linear/non-linear) or
  "mfuzz" (fuzzy clustering) (default: "fitting")

- fitting_paras:

  List of parameters for fitting method:

  - lm_min_slope: Minimum slope threshold for significant linear trend

  - lm_pval_cutoff: P-value cutoff for linear regression significance

- mfuzz_paras:

  List of parameters for Mfuzz method:

  - cluster_num: Number of clusters for Mfuzz (default: 10)

  - min_acore: Minimum membership score for gene inclusion (default:
    0.3)

- gene_list:

  Character vector, specific genes to analyze (default: NULL, uses all
  variable features)

- remove_genes:

  Character vector, genes to exclude from analysis

- col:

  Character vector, color palette for visualization

- return_data:

  Logical, whether to return results (default: FALSE)

- grp_nm:

  Character, group name for output (default: NULL, uses timestamp)

- dir_nm:

  Character, directory name for output (default: "M4_CalSampGeneTrend")

## Value

If return_data = TRUE, returns a list of gene trend results per cell
type

## Details

The fitting method uses linear regression to identify genes with
significant monotonic trends (increasing/decreasing) and quadratic
regression to classify concave/convex patterns. The Mfuzz method
performs fuzzy clustering to group genes with similar temporal
expression patterns.

## Examples

``` r
if (FALSE) { # \dontrun{
# Linear regression based trend analysis
results <- CalSampGeneTrend(
  STID_obj = stid_object,
  loop_id = "time_series",
  group_by = "cell_type",
  method = "fitting",
  fitting_paras = list(lm_min_slope = 1.5, lm_pval_cutoff = 0.05)
)

# Mfuzz clustering based trend analysis
results <- CalSampGeneTrend(
  STID_obj = stid_object,
  loop_id = "time_series",
  group_by = "cell_type",
  method = "mfuzz",
  mfuzz_paras = list(cluster_num = 8, min_acore = 0.4)
)
} # }
```
