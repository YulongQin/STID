# Plot Gene Expression Along Distance from Niche Center

Visualizes gene expression patterns as a function of distance from niche
centers, with optional smoothing and scaling.

## Usage

``` r
Plot_DistLine_Exp(
  STID_obj = NULL,
  loop_id = "LoopAllSamp",
  meta_key = NULL,
  facet_grpnm = NULL,
  assay_id = "Spatial",
  layer_id = "data",
  features = NULL,
  feature_colnm = NULL,
  smooth_method = "gam",
  exp_scale = TRUE,
  distance_scale = TRUE,
  linewidth = 1,
  ncol = 4,
  col = COLOR_LIST[["PALETTE_WHITE_BG"]]
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- meta_key:

  Character, metadata key containing distance information

- facet_grpnm:

  Character, column name for faceting

- assay_id:

  Character, assay name (default: "Spatial")

- layer_id:

  Character, layer name (default: "data")

- features:

  Character vector, gene names to plot

- feature_colnm:

  Character vector, metadata column names to plot

- smooth_method:

  Character, smoothing method for geom_smooth (default: "gam")

- exp_scale:

  Logical, whether to scale expression values (default: TRUE)

- distance_scale:

  Logical, whether to scale distance values (default: TRUE)

- linewidth:

  Numeric, line width (default: 1)

- ncol:

  Integer, number of facet columns (default: 4)

- col:

  Color palette (default: COLOR_LIST\$PALETTE_WHITE_BG)

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
