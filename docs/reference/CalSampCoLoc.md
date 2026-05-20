# Calculate Niche Co-localization using MSTIDyR

Performs spatial co-localization analysis using the MSTIDyR framework to
identify intra-view, juxtaview, and paraview relationships between cell
types or gene expression patterns within niches.

## Usage

``` r
CalSampCoLoc(
  STID_obj = NULL,
  loop_id = "LoopAllSamp",
  meta_key = NULL,
  niche_key = NULL,
  group_by = NULL,
  group_use = NULL,
  features = NULL,
  feature_colnm = NULL,
  method = "mistyR",
  mistyR_params = c(juxtaview_radius = 15, paraview_radius = 10),
  return_data = TRUE,
  grp_nm = NULL,
  dir_nm = "M3_CalSampCoLoc"
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- meta_key:

  Character, metadata key for analysis when niche_key is NULL

- niche_key:

  Character, niche key to analyze

- group_by:

  Character, column name for cell type grouping

- group_use:

  Character vector, specific groups to include

- features:

  Character vector, gene names to analyze

- feature_colnm:

  Character vector, metadata column names to analyze

- method:

  Character, method for co-localization analysis - "mistyR" or "squidpy"
  (default: "mistyR")

- mistyR_params:

  Numeric vector, parameters for MSTIDyR analysis (default:
  c(juxtaview_radius = 15, paraview_radius = 10)) - juxtaview_radius:
  radius for juxtaview in spatial units - paraview_radius: radius for
  paraview in spatial units

- return_data:

  Logical, whether to return the results list (default: TRUE)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M3_CalSampCoLoc")

## Value

If return_data = TRUE, returns a list of MSTIDyR results per sample

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze co-localization of cell types within niches
results <- CalSampCoLoc(
  STID_obj = STID_obj,
  niche_key = "niche_virulence",
  group_by = "cell_type",
  mistyR_params = c(15,10)
)
} # }
```
