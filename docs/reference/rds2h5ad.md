# Convert Seurat rds file to AnnData h5ad file

Converts a Seurat rds file to a Python AnnData h5ad file for use with
scanpy and other Python-based tools.

## Usage

``` r
rds2h5ad(
  seurat_obj = NULL,
  data_type = "stRNA",
  convert_mode = "seurat",
  assay_id = NULL,
  grp_nm = NULL,
  dir_nm = "M0_rds2h5ad"
)
```

## Arguments

- seurat_obj:

  A Seurat object to convert

- data_type:

  Character, data type - "stRNA" (spatial) or "scRNA" (single-cell)
  (default: "stRNA")

- convert_mode:

  Character, conversion method - "seurat" (using SeuratDisk) (default:
  "seurat")

- assay_id:

  Character, assay name to convert (default: NULL, auto-detected)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M0_rds2h5ad")

## Value

NULL (invisible), saves h5ad file to disk

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert Seurat object to h5ad
rds2h5ad(
  seurat_obj = seurat_object,
  python_path = "~/anaconda3/envs/scanpy/bin/python",
  data_type = "stRNA",
  grp_nm = "sample1"
)
} # }
```
