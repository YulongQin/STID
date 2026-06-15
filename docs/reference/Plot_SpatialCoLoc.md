# Plot Spatial Co-localization

Visualizes spatial co-localization of two cell groups, genes, or
metadata features in selected samples or niches.

## Usage

``` r
Plot_SpatialCoLoc(
  STID_obj = NULL,
  loop_id = "LoopAllSamp",
  meta_key = NULL,
  niche_key = NULL,
  group_by = NULL,
  group_use = NULL,
  features = NULL,
  feature_colnm = NULL,
  exp_thres = 1,
  col = COLOR_LIST$PALETTE_WHITE_BG,
  pt_size = 1.5
)
```

## Arguments

- STID_obj:

  An STID object.

- loop_id:

  Character, sample loop identifier (default: `"LoopAllSamp"`).

- meta_key:

  Character, metadata key used when `niche_key` is `NULL`.

- niche_key:

  Character, niche key used to retrieve niche cells.

- group_by:

  Character, metadata column for group-based co-localization.

- group_use:

  Character vector, groups to keep from `group_by`.

- features:

  Character vector, feature or gene names for co-localization.

- feature_colnm:

  Character vector, metadata columns used as features.

- exp_thres:

  Numeric, threshold below which feature values are set to zero
  (default: `1`).

- col:

  Character vector, color palette (default:
  `COLOR_LIST$PALETTE_WHITE_BG`).

- pt_size:

  Numeric, point size for spatial plots (default: `1.5`).

## Value

Invisibly returns `NULL`; spatial plots are printed.

## Examples

``` r
if (FALSE) { # \dontrun{
Plot_SpatialCoLoc(
STID_obj = STID_obj,
niche_key = "niche_virulence",
group_by = "cell_type",
group_use = c("Macrophage", "T cell")
)

Plot_SpatialCoLoc(
STID_obj = STID_obj,
niche_key = "niche_virulence",
features = c("Cxcl10", "Isg15")
)
} # }
```
