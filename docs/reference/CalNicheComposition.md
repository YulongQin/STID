# Calculate Niche Cell Type Composition

Analyzes and visualizes the cell type composition within and outside
niches. Generates stacked barplots and faceted barplots showing cell
type proportions for niche vs. bystander regions.

## Usage

``` r
CalNicheComposition(
  STID_obj = NULL,
  samp_mode = "SS",
  loop_id = "LoopAllSamp",
  samp_grp_index = FALSE,
  meta_key = NULL,
  niche_key = NULL,
  group_by = NULL,
  col = COLOR_List[["PALETTE_WHITE_BG"]],
  return_data = FALSE
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- samp_mode:

  Character, sample type - "SS" (single sample) or "MS" (multi-sample)
  (default: "SS")

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- samp_grp_index:

  Logical, whether to group by sample groups in MS mode (default: FALSE)

- meta_key:

  Character, metadata key for MS mode when niche_key is NULL

- niche_key:

  Character, niche key to analyze (only one value supported)

- group_by:

  Character, column name for cell type grouping

- col:

  Color palette for visualization (default:
  COLOR_List\$PALETTE_WHITE_BG)

- return_data:

  Logical, whether to return the plot data (default: FALSE)

## Value

If return_data = TRUE, returns a list of plots per sample; otherwise
NULL

## Examples

``` r
if (FALSE) { # \dontrun{
# Single-sample niche composition
CalNicheComposition(
  STID_obj = STID_obj,
  samp_mode  = "SS",
  niche_key = "niche_virulence",
  group_by = "cell_type"
)

# Multi-sample niche composition with sample grouping
CalNicheComposition(
  STID_obj = STID_obj,
  samp_mode  = "MS",
  loop_id = "LoopAllMulti",
  niche_key = "niche_virulence",
  samp_grp_index = TRUE,
  group_by = "cell_type"
)
} # }
```
