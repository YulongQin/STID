# Convert AnnData h5ad file to Seurat rds file

Converts a Python AnnData h5ad file to a Seurat rds file, preserving
spatial coordinates, dimensional reductions, and optionally
SCT-transformed data.

## Usage

``` r
h5ad2rds(
  file_path = NULL,
  data_type = "stRNA",
  convert_mode = "scanpy",
  assay_id = NULL,
  X_index = "rawX",
  binsize = 1,
  SCT_index = FALSE,
  reduction_index = FALSE,
  image_index = NULL,
  return_object = TRUE,
  grp_nm = NULL,
  dir_nm = "M0_h5ad2rds"
)
```

## Arguments

- file_path:

  Character, path to the input h5ad file

- data_type:

  Character, data type - "stRNA" (spatial) or "scRNA" (single-cell)
  (default: "stRNA")

- convert_mode:

  Character, conversion method - "scanpy" or "seurat" (default:
  "scanpy")

- assay_id:

  Character, assay name for the Seurat object (default: NULL,
  auto-detected based on data_type)

- X_index:

  Character, which matrix to use as counts - "X" or "rawX" (default:
  "rawX")

- binsize:

  Numeric, bin size for spatial coordinate scaling (default: 1)

- SCT_index:

  Logical, whether to convert SCT results from stereopy (default: FALSE)

- reduction_index:

  Logical, whether to preserve dimensional reductions (default: FALSE)

- image_index:

  Logical, whether to create spatial image object (default: NULL,
  auto-detected based on data_type)

- return_object:

  Logical, whether to return the Seurat object (default: TRUE)

- grp_nm:

  Character, group name for output organization (default: "sample1")

- dir_nm:

  Character, directory name for output (default: "M0_h5ad2rds")

## Value

If return_object = TRUE, returns a Seurat object; otherwise returns NULL

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert spatial transcriptomics h5ad to Seurat
configure_conda()
seurat_obj <- h5ad2rds(
  file_path = "data/spatial_data.h5ad",
  data_type = "stRNA",
  binsize = 100,
  reduction_index = TRUE,
  image_index = TRUE
)
} # }
```
