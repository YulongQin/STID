# Plot Gene Expression Along Distance from Niche Center

Visualizes gene expression patterns as a function of distance from niche
centers, with optional smoothing and scaling.

## Usage

``` r
Plot_DistLine_Exp(
  STID_obj = NULL,
  features = NULL,
  feature_colnm = NULL,
  facet_grpnm = NULL,
  meta_key = NULL,
  loop_id = "LoopAllSamp",
  smooth_method = "gam",
  exp_scale = TRUE,
  distance_scale = TRUE,
  col = COLOR_List[["PALETTE_WHITE_BG"]],
  linewidth = 1,
  ncol = 4,
  assay_id = "Spatial",
  layer_id = "data"
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- features:

  Character vector, gene names to plot

- feature_colnm:

  Character vector, metadata column names to plot

- facet_grpnm:

  Character, column name for faceting

- meta_key:

  Character, metadata key containing distance information

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- smooth_method:

  Character, smoothing method for geom_smooth (default: "gam")

- exp_scale:

  Logical, whether to scale expression values (default: TRUE)

- distance_scale:

  Logical, whether to scale distance values (default: TRUE)

- col:

  Color palette (default: COLOR_List\$PALETTE_WHITE_BG)

- linewidth:

  Numeric, line width (default: 1)

- ncol:

  Integer, number of facet columns (default: 4)

- assay_id:

  Character, assay name (default: "Spatial")

- layer_id:

  Character, layer name (default: "data")

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot gene expression along distance from niche center
Plot_DistLine_Exp(
  STID_obj = STID_obj,
  features = c("gene1", "gene2"),
  meta_key = "M2_NicheDetect_STS_20240101"
)
} # }
```
