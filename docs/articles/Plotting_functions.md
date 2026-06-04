# Plotting functions

## Overview

This vignette introduces the public plotting functions whose names start
with `Plot_` in the STID analysis scripts:

- [`Plot_Spatial()`](https://yulongqin.github.io/STID/reference/Plot_Spatial.md)
  for spatial visualization of discrete labels or continuous scores;
- `Plot_SpatialCoLoc()` for visualizing two-variable spatial
  co-localization;
- [`Plot_NicheCellComm()`](https://yulongqin.github.io/STID/reference/Plot_NicheCellComm.md)
  for visualizing niche cell-cell communication results;
- [`Plot_DistLine_Exp()`](https://yulongqin.github.io/STID/reference/Plot_DistLine_Exp.md)
  for plotting expression trends along distance from ROI/niche centers;
- [`Plot_DistLine_Ratio()`](https://yulongqin.github.io/STID/reference/Plot_DistLine_Ratio.md)
  for plotting cell-type ratio trends along distance from ROI/niche
  centers.

> **Note:** The examples are shown with `eval = FALSE` because STID
> objects, metadata keys, niche keys, CellChat results, gene names, and
> output folders are dataset-specific. Replace the placeholder values
> before running the examples.

## Prerequisites

``` r
library(tidyverse)
library(Seurat)
library(STID)
```

Depending on which plotting function is used, the following packages may
also be required.

``` r
# Common plotting/data manipulation dependencies
library(dplyr)
library(ggh4x)
library(scales)

# Cell communication visualization
library(CellChat)
library(nichenetr)
```

## Plotting-function map

| Function | Main purpose | Typical use |
|----|----|----|
| [`Plot_Spatial()`](https://yulongqin.github.io/STID/reference/Plot_Spatial.md) | Create spatial scatter plots from an STID object or a custom data frame | Visualize cell types, spot classes, gene scores, or continuous metadata on tissue coordinates |
| `Plot_SpatialCoLoc()` | Visualize spatial co-localization of two groups, genes, or metadata features | Check whether two cell groups or two expression features overlap in niche/ROI regions |
| [`Plot_NicheCellComm()`](https://yulongqin.github.io/STID/reference/Plot_NicheCellComm.md) | Plot CellChat-based niche cell communication results | Inspect interaction counts, interaction weights, signaling roles, pathway bubbles, and spatial signaling maps |
| [`Plot_DistLine_Exp()`](https://yulongqin.github.io/STID/reference/Plot_DistLine_Exp.md) | Plot gene or metadata-feature expression along distance from ROI/niche centers | Examine how expression changes from ROI center toward the edge and outside regions |
| [`Plot_DistLine_Ratio()`](https://yulongqin.github.io/STID/reference/Plot_DistLine_Ratio.md) | Plot cell-type proportions along distance from ROI/niche centers | Examine how cell-type composition changes across binned distance layers |

## Spatial plots

[`Plot_Spatial()`](https://yulongqin.github.io/STID/reference/Plot_Spatial.md)
creates spatial coordinate plots with either discrete or continuous
coloring. It can read metadata directly from an STID object, or it can
plot a custom data frame if coordinate columns and a grouping column are
supplied.

``` r
# Plot default metadata from an STID object
p <- Plot_Spatial(
  STID_obj = STID_obj,
  meta_key = "coord",
  group_by = "cell_type",
  datatype = "discrete",
  pt_size = 1.1,
  title = "Cell type distribution"
)

p
```

``` r
# Plot a continuous score with percentile clipping
p <- Plot_Spatial(
  STID_obj = STID_obj,
  meta_key = "coord",
  group_by = "Pathogen_Score",
  datatype = "continuous",
  vmin = "p1",
  vmax = "p99",
  black_bg = TRUE,
  title = "Pathogen score"
)

p
```

``` r
# Plot from a custom data frame
p <- Plot_Spatial(
  plot_data = my_spatial_df,
  x_colnm = "x",
  y_colnm = "y",
  group_by = "ROI_label",
  facet_grpnm = "sample",
  datatype = "discrete",
  pt_size = 1.2
)

p
```

## Spatial co-localization plots

`Plot_SpatialCoLoc()` visualizes co-localization for exactly two
variables. The variables can be two groups from a `group_by` column, two
genes from `features`, two metadata features from `feature_colnm`, or a
combination of gene and metadata features. When `niche_key` is supplied,
the function focuses on niche cells; otherwise, it can use metadata
through `meta_key`.

``` r
# Co-localization of selected cell groups inside a niche
Plot_SpatialCoLoc(
  STID_obj = STID_obj,
  loop_id = "LoopAllSamp",
  niche_key = "niche_virulence",
  group_by = "cell_type",
  group_use = c("Macrophage", "Epithelial"),
  col = COLOR_LIST$PALETTE_WHITE_BG,
  pt_size = 1.5
)
```

``` r
# Co-localization of two gene-expression features
Plot_SpatialCoLoc(
  STID_obj = STID_obj,
  loop_id = "LoopAllSamp",
  meta_key = "coord",
  features = c("GeneA", "GeneB"),
  exp_thres = 1,
  col = COLOR_LIST$PALETTE_WHITE_BG,
  pt_size = 1.5
)
```

``` r
# Co-localization using two numeric metadata columns
Plot_SpatialCoLoc(
  STID_obj = STID_obj,
  meta_key = "M2_NicheDetect_Lasso_20240101",
  feature_colnm = c("Score_A", "Score_B"),
  exp_thres = 1,
  pt_size = 1.5
)
```

> **Note:** This function is designed for two-variable co-localization.
> If more than two variables are supplied, the function stops with an
> error.

## Niche cell communication plots

[`Plot_NicheCellComm()`](https://yulongqin.github.io/STID/reference/Plot_NicheCellComm.md)
visualizes CellChat results generated from niche cell communication
analysis. In single-sample mode, it plots interaction circle plots,
heatmaps, signaling-role heatmaps, spatial signaling maps when
available, and ligand-receptor bubble plots. In multi-sample mode, it
compares CellChat objects across samples in a multi-sample loop.

``` r
# Single-sample CellChat visualization
Plot_NicheCellComm(
  STID_obj = STID_obj,
  samp_mode = "SS",
  loop_id = "LoopAllSamp",
  CellComm_data = cellcomm_results,
  signaling = c("TGFb", "WNT"),
  sources.use = c("Epithelial"),
  targets.use = c("Immune", "Macrophage")
)
```

``` r
# Multi-sample CellChat comparison
Plot_NicheCellComm(
  STID_obj = STID_obj,
  samp_mode = "MS",
  loop_id = "LoopAllMulti",
  CellComm_data = cellcomm_results,
  signaling = c("TGFb", "WNT"),
  sources.use = NULL,
  targets.use = NULL
)
```

If a custom color vector is supplied through `col`, it should be a named
vector whose names match the cell-type names in the CellChat object.

``` r
celltype_cols <- c(
  Epithelial = "#4DBBD5FF",
  Macrophage = "#E64B35FF",
  Immune = "#00A087FF"
)

Plot_NicheCellComm(
  STID_obj = STID_obj,
  samp_mode = "SS",
  CellComm_data = cellcomm_results,
  signaling = "TGFb",
  col = celltype_cols
)
```

## Expression trends along distance

[`Plot_DistLine_Exp()`](https://yulongqin.github.io/STID/reference/Plot_DistLine_Exp.md)
plots expression values along distance from ROI/niche centers. It
expects a metadata table containing distance-related columns such as
`ROI_edge`, `All_Dist2ROIcenter`, and `All_Dist2ROIedge`, which are
typically produced by niche/ROI detection workflows.

The function can plot genes from the assay layer, numeric columns from
metadata, or both.

``` r
# Plot gene-expression trends along distance from ROI center
Plot_DistLine_Exp(
  STID_obj = STID_obj,
  loop_id = "LoopAllSamp",
  meta_key = "M2_NicheDetect_Lasso_20240101",
  facet_grpnm = "sample",
  assay_id = "Spatial",
  layer_id = "data",
  features = c("GeneA", "GeneB"),
  smooth_method = "gam",
  exp_scale = TRUE,
  distance_scale = TRUE,
  linewidth = 1,
  col = COLOR_LIST$PALETTE_WHITE_BG
)
```

``` r
# Plot numeric metadata features along distance from ROI center
Plot_DistLine_Exp(
  STID_obj = STID_obj,
  meta_key = "M2_NicheDetect_Lasso_20240101",
  facet_grpnm = "sample",
  feature_colnm = c("Pathogen_Score", "Inflammation_Score"),
  smooth_method = "loess",
  exp_scale = TRUE,
  distance_scale = TRUE
)
```

Use `distance_scale = TRUE` to scale distance by the 95th percentile
within each sample. Use `exp_scale = TRUE` to standardize expression or
score values within each sample and feature.

## Cell-type ratio trends along distance

[`Plot_DistLine_Ratio()`](https://yulongqin.github.io/STID/reference/Plot_DistLine_Ratio.md)
summarizes spot/cell composition by distance bins and plots the fraction
of selected cell types across distance layers. The bin width is
controlled by `coord_interval_ratio`, which multiplies the
sample-specific spatial coordinate interval.

``` r
# Plot cell-type ratio trends across distance bins
Plot_DistLine_Ratio(
  STID_obj = STID_obj,
  loop_id = "LoopAllSamp",
  meta_key = "M2_NicheDetect_Lasso_20240101",
  group_by = "cell_type",
  facet_grpnm = "sample",
  celltypes = c("Tcell", "Bcell", "Macrophage"),
  coord_interval_ratio = 5,
  linewidth = 1,
  col = COLOR_LIST$PALETTE_WHITE_BG
)
```

``` r
# Plot all cell types in the grouping column
Plot_DistLine_Ratio(
  STID_obj = STID_obj,
  meta_key = "M2_NicheDetect_Lasso_20240101",
  group_by = "cell_type",
  facet_grpnm = "sample",
  celltypes = NULL,
  coord_interval_ratio = 5
)
```

The shaded region and dashed line represent the median ROI-edge distance
bin for each sample, helping compare composition inside and outside the
ROI boundary.

## Suggested workflow

The plotting functions are often used after spot detection, ROI/niche
detection, and downstream niche analysis.

``` r
# 1. Visualize spatial metadata or scores
Plot_Spatial(
  STID_obj = STID_obj,
  meta_key = "coord",
  group_by = "cell_type",
  datatype = "discrete"
)

# 2. Visualize two-variable co-localization
Plot_SpatialCoLoc(
  STID_obj = STID_obj,
  niche_key = "niche_virulence",
  group_by = "cell_type",
  group_use = c("Macrophage", "Epithelial")
)

# 3. Plot communication results from CalSampCellComm()
Plot_NicheCellComm(
  STID_obj = STID_obj,
  CellComm_data = cellcomm_results,
  signaling = c("TGFb", "WNT")
)

# 4. Plot expression changes along ROI-center distance
Plot_DistLine_Exp(
  STID_obj = STID_obj,
  meta_key = "M2_NicheDetect_Lasso_20240101",
  facet_grpnm = "sample",
  features = c("GeneA", "GeneB")
)

# 5. Plot cell-type composition changes along ROI-center distance
Plot_DistLine_Ratio(
  STID_obj = STID_obj,
  meta_key = "M2_NicheDetect_Lasso_20240101",
  group_by = "cell_type",
  facet_grpnm = "sample",
  coord_interval_ratio = 5
)
```

## Session information

``` r
sessionInfo()
#> R version 4.2.0 (2022-04-22 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 22000)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
#> [2] LC_CTYPE=Chinese (Simplified)_China.utf8   
#> [3] LC_MONETARY=Chinese (Simplified)_China.utf8
#> [4] LC_NUMERIC=C                               
#> [5] LC_TIME=Chinese (Simplified)_China.utf8    
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.35     R6_2.6.1          jsonlite_1.8.8    lifecycle_1.0.5  
#>  [5] evaluate_1.0.1    cachem_1.1.0      rlang_1.1.7       cli_3.6.5        
#>  [9] rstudioapi_0.15.0 fs_1.6.3          jquerylib_0.1.4   bslib_0.8.0      
#> [13] ragg_1.3.0        rmarkdown_2.29    pkgdown_2.2.0     textshaping_0.3.6
#> [17] desc_1.4.3        tools_4.2.0       htmlwidgets_1.6.4 yaml_2.3.10      
#> [21] xfun_0.49         fastmap_1.2.0     compiler_4.2.0    systemfonts_1.0.4
#> [25] htmltools_0.5.8.1 knitr_1.49        sass_0.4.9
```
