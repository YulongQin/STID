# Plot Cell Type Ratio Along Distance from Niche Center

Visualizes the fraction of cell types as a function of distance from
niche centers, with distance binning based on spatial intervals.

## Usage

``` r
Plot_DistLine_Ratio(
  STID_obj = NULL,
  celltypes = NULL,
  group_by = NULL,
  facet_grpnm = NULL,
  meta_key = NULL,
  loop_id = "LoopAllSamp",
  interval_ratio = NULL,
  col = COLOR_List[["PALETTE_WHITE_BG"]],
  linewidth = 1,
  ncol = 4
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- celltypes:

  Character vector, cell types to plot

- group_by:

  Character, column name for cell type grouping

- facet_grpnm:

  Character, column name for faceting

- meta_key:

  Character, metadata key containing distance information

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- interval_ratio:

  Numeric, multiplier for interval to set bin width

- col:

  Color palette (default: COLOR_List\$PALETTE_WHITE_BG)

- linewidth:

  Numeric, line width (default: 1)

- ncol:

  Integer, number of facet columns (default: 4)

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
  interval_ratio = 5
)
} # }
```
