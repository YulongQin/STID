# Automated ROI detection using density-based spatial clustering

This function automatically detects regions of interest (ROIs) using
density-based clustering (DBSCAN) on spatial transcriptomics data. It
includes preprocessing steps like density filtering and iterative spot
updating.

## Usage

``` r
NicheDetect_STS(
  STID_obj = NULL,
  meta_key = NULL,
  loop_id = "LoopAllSamp",
  pos_colnm = NULL,
  neg_value = "neg",
  spatial_scale_method = "region",
  region_detect_method = "convex",
  concavity = 2,
  update_spots = TRUE,
  density_thres = 0.9,
  ROI_size = 10,
  minPts = NULL,
  k_kNNdist = NULL,
  description = NULL,
  grp_nm = NULL,
  dir_nm = "M2_NicheDetect_STS"
)
```

## Arguments

- STID_obj:

  An STID object containing spatial transcriptomics data

- meta_key:

  Character, metadata key containing positive spot information

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- pos_colnm:

  Character, column name containing positive spot labels

- neg_value:

  Character, value indicating negative/non-ROI spots (default: "neg")

- spatial_scale_method:

  Character, detection level - "region" or "spot" (default: "region")

- region_detect_method:

  Character, hull method - "convex" or "concave" (default: "convex")

- concavity:

  Numeric, concavity parameter for concave hull (default: 2)

- update_spots:

  Logical, whether to iteratively update positive spots (default: TRUE)

- density_thres:

  Numeric, density threshold for filtering (0-1, default: 0.9)

- ROI_size:

  Integer, minimum number of spots to form an ROI (default: 10)

- minPts:

  Integer, minimum points parameter for DBSCAN (default: NULL,
  automatically determined based on data format)

- k_kNNdist:

  Integer, k value for kNN distance calculation (default: NULL, uses
  minPts if not specified)

- description:

  Character, description of the analysis (default: NULL)

- grp_nm:

  Character, group name for output organization (default: NULL, uses
  timestamp)

- dir_nm:

  Character, directory name for output (default: "M2_NicheDetect_STS")

## Value

Returns the modified STID object with added metadata containing ROI
labels, distances to ROI centers, edge information, and region
classifications

## Examples

``` r
if (FALSE) { # \dontrun{
# Automatically detect ROIs based on positive spots
STID_obj <- NicheDetect_STS(
  STID_obj = ist_object,
  meta_key = "M1_SpotDetect_Gene_20240101",
  pos_colnm = "Label_geneA",
  density_thres = 0.9,
  ROI_size = 10
)
} # }
```
