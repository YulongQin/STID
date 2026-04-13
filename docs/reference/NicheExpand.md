# Expand niche regions based on distance threshold

Expands existing niche regions by including neighboring spots within a
specified distance threshold. Useful for defining broader influence
zones around detected niches.

## Usage

``` r
NicheExpand(
  STID_obj = NULL,
  meta_key = NULL,
  loop_id = "LoopAllSamp",
  pos_colnm = NULL,
  neg_value = "neg",
  center_colnm = NULL,
  expand_dist = 1,
  description = NULL,
  grp_nm = NULL,
  dir_nm = "M2_NicheExpand"
)
```

## Arguments

- STID_obj:

  An STID object containing niche detection results

- meta_key:

  Character, metadata key containing niche information

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- pos_colnm:

  Character, column name containing positive spot labels

- neg_value:

  Character, value indicating negative spots (default: "neg")

- center_colnm:

  Character, column name indicating niche center spots (default: NULL)

- expand_dist:

  Numeric, distance threshold for expansion (default: 1)

- description:

  Character, description of the analysis

- grp_nm:

  Character, group name for output organization (default: NULL, uses
  timestamp)

- dir_nm:

  Character, directory name for output (default: "M2_NicheExpand")

## Value

Returns the modified STID object with expanded ROI labels and updated
distance calculations

## Examples

``` r
if (FALSE) { # \dontrun{
# Expand niches by including spots within distance 5
STID_obj <- NicheExpand(
  STID_obj = ist_object,
  meta_key = "M2_NicheDetect_STS_20240101",
  pos_colnm = "ROI_label",
  expand_dist = 5
)
} # }
```
