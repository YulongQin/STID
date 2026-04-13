# Interactive lasso selection for ROI (Region of Interest) detection

This function launches an interactive Shiny application that allows
users to manually draw lasso selections to define regions of interest
(ROIs) on spatial transcriptomics data. Supports multiple ROIs per
sample and saves the results for downstream analysis.

## Usage

``` r
NicheDetect_Lasso(
  STID_obj = NULL,
  loop_id = "LoopAllSamp",
  meta_key = "coord",
  group_by = NULL,
  col = COLOR_LIST$PALETTE_WHITE_BG,
  description = NULL,
  grp_nm = NULL,
  dir_nm = "M2_NicheDetect_Lasso"
)
```

## Arguments

- STID_obj:

  An STID object containing spatial transcriptomics data

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- meta_key:

  Character, metadata key containing spatial coordinates and grouping
  information (default: "coord")

- group_by:

  Character, column name for coloring points in the interactive plot

- col:

  Color palette for visualization (default:
  COLOR_LIST\$PALETTE_WHITE_BG)

- description:

  Character, description of the analysis (default: NULL)

- grp_nm:

  Character, group name for output organization (default: NULL, uses
  timestamp)

- dir_nm:

  Character, directory name for output (default: "M2_NicheDetect_Lasso")

## Value

Returns the modified STID object with added metadata containing ROI
labels, distances to ROI centers, and edge information

## Examples

``` r
if (FALSE) { # \dontrun{
# Launch interactive ROI selection
STID_obj <- NicheDetect_Lasso(
  STID_obj = ist_object,
  group_by = "cell_type",
  meta_key = "coord"
)
} # }
```
