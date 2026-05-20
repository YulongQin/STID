# Calculate Niche Cell Communication using CellChat

Analyzes cell-cell communication within niches using the CellChat
framework. Supports both spatial and non-spatial modes with
comprehensive visualization.

## Usage

``` r
CalSampCellComm(
  STID_obj = NULL,
  loop_id = "LoopAllSamp",
  niche_key = NULL,
  meta_key = NULL,
  group_by = NULL,
  assay_id = "Spatial",
  layer_id = "data",
  is_Spatial = TRUE,
  spatial.factors = NULL,
  interaction.range = 250,
  remove_genes = NULL,
  col = NULL,
  return_data = TRUE,
  grp_nm = NULL,
  dir_nm = "M3_CalSampCellComm"
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- niche_key:

  Character, niche key to analyze (default: NULL)

- meta_key:

  Character, metadata key for cell annotations (default: NULL)

- group_by:

  Character, column name for cell type grouping

- assay_id:

  Character, assay name (default: "Spatial")

- layer_id:

  Character, layer name (default: "data")

- is_Spatial:

  Logical, whether to use spatial mode (default: TRUE)

- spatial.factors:

  Data frame, spatial scaling factors (default: NULL)

- interaction.range:

  Numeric, interaction range for spatial mode (default: 250)

- remove_genes:

  Character vector, genes to exclude

- col:

  Named vector, colors for cell types

- return_data:

  Logical, whether to return results list (default: TRUE)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M3_CalSampCellComm")

## Value

If return_data = TRUE, returns a list of CellChat objects per sample

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze cell communication in niches
results <- CalSampCellComm(
  STID_obj = STID_obj,
  niche_key = "niche_virulence",
  group_by = "cell_type",
  is_Spatial = TRUE,
  interaction.range = 200
)
} # }
```
