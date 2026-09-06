# Plot Cell Type Ratio Along Distance from Niche Center

Visualizes the fraction of cell types as a function of distance from
niche centers, with distance binning based on spatial coord_intervals.

## Usage

``` r
Plot_DistLine_Ratio(
  STID_obj = NULL,
  loop_id = "LoopAllSamp",
  meta_key = NULL,
  group_by = NULL,
  facet_grpnm = NULL,
  celltypes = NULL,
  coord_interval_ratio = NULL,
  linewidth = 1,
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

- group_by:

  Character, column name for cell type grouping

- facet_grpnm:

  Character, column name for faceting

- celltypes:

  Character vector, cell types to plot

- coord_interval_ratio:

  Numeric, multiplier for coord_interval to set bin width

- linewidth:

  Numeric, line width (default: 1)

- col:

  Color palette (default: COLOR_LIST\$PALETTE_WHITE_BG)

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot cell type ratios along distance from niche center
Plot_DistLine_Ratio(
  STID_obj = STID_obj,
  celltypes = c("Tcell", "Bcell", "Macrophage"),
  group_by = "cell_type",
  meta_key = "M2_NicheDetect_STS_20240101",
  coord_interval_ratio = 5
)
} # }
```
