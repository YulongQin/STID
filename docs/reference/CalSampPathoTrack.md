# Calculate pathogen tracking across samples

Analyzes pathogen distribution across temporal or comparative samples,
infers dissemination networks, and visualizes pathogen spread patterns.

## Usage

``` r
CalSampPathoTrack(
  STID_obj = NULL,
  loop_id = NULL,
  samp_grp_index = FALSE,
  meta_key = NULL,
  niche_key = NULL,
  group_by = NULL,
  pos_colnm = NULL,
  neg_value = "neg",
  col = COLOR_LIST[["PALETTE_WHITE_BG"]],
  return_data = FALSE,
  grp_nm = NULL,
  dir_nm = "M4_CalSampPathoTrack"
)
```

## Arguments

- STID_obj:

  A STID object

- loop_id:

  Character, multi-sample analysis identifier

- samp_grp_index:

  Logical, whether to group by sample groups (default: FALSE)

- meta_key:

  Character, metadata key for retrieving cell data

- niche_key:

  Character, niche key for filtering niche cells (optional)

- group_by:

  Character, column name for cell type grouping

- pos_colnm:

  Character, column name containing positive spot labels

- neg_value:

  Character, value indicating negative spots (default: "neg")

- col:

  Character vector, color palette for cell types

- return_data:

  Logical, whether to return results (default: FALSE)

- grp_nm:

  Character, group name for output (default: NULL, uses timestamp)

- dir_nm:

  Character, directory name for output (default: "M4_CalSampPathoTrack")

## Value

If return_data = TRUE, returns a list containing:

- pathogen_data: Pathogen-positive cell counts and ratios per time point

- pathogen_wide: Wide format pathogen ratio matrix

- region_summary: Peak time and total burden per cell type

- edges_df: Inferred dissemination network edges

- spatiotemporal_edges: Time-resolved transition edges

## Details

The function performs three types of analyses:

1.  Temporal trend analysis of pathogen-positive cell type composition

2.  Peak-based dissemination network inference (which cell types seed
    others)

3.  Transition-based spatiotemporal mapping of pathogen spread

## Examples

``` r
if (FALSE) { # \dontrun{
# Track pathogen spread across time points
results <- CalSampPathoTrack(
  STID_obj = stid_object,
  loop_id = "infection_series",
  group_by = "cell_type",
  pos_colnm = "Label_pathogen"
)
} # }
```
