# Calculate Niche Cell Communication using CellChat

Analyzes cell-cell communication within niches using the CellChat
framework. Supports both spatial and non-spatial modes with
comprehensive visualization.

## Usage

``` r
CalNicheCellComm(
  STID_obj = NULL,
  niche_key = NULL,
  group_by = NULL,
  assay_id = "Spatial",
  layer_id = "data",
  loop_id = "LoopAllSamp",
  col = NULL,
  is_Spatial = TRUE,
  spatial.factors = NULL,
  interaction.range = 250,
  remove_genes = NULL,
  return_data = TRUE,
  grp_nm = NULL,
  dir_nm = "M3_CalNicheCellComm"
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- niche_key:

  Character, niche key to analyze

- group_by:

  Character, column name for cell type grouping

- assay_id:

  Character, assay name (default: "Spatial")

- layer_id:

  Character, layer name (default: "data")

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- col:

  Named vector, colors for cell types

- is_Spatial:

  Logical, whether to use spatial mode (default: TRUE)

- spatial.factors:

  Data frame, spatial scaling factors (default: NULL)

- interaction.range:

  Numeric, interaction range for spatial mode (default: 250)

- remove_genes:

  Character vector, genes to exclude

- return_data:

  Logical, whether to return results list (default: TRUE)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M3_CalNicheCellComm")

## Value

If return_data = TRUE, returns a list of CellChat objects per sample

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze cell communication in niches
results <- CalNicheCellComm(
  STID_obj = STID_obj,
  niche_key = "niche_virulence",
  group_by = "cell_type",
  is_Spatial = TRUE,
  interaction.range = 200
)
} # }
```
