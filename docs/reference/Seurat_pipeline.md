# Standard Seurat Analysis Pipeline for Spatial Transcriptomics

Performs a complete Seurat analysis pipeline for spatial transcriptomics
data, including normalization, feature selection, scaling, PCA,
clustering, and UMAP/t-SNE dimensionality reduction. Also generates
spatial feature plots.

## Usage

``` r
Seurat_pipeline(
  seurat_obj = NULL,
  data_format = "stereo",
  data_type = "stRNA",
  resolution_index = seq(0.1, 1.3, 0.2),
  runTSNE_index = FALSE,
  assay_nm = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object containing spatial transcriptomics data

- data_format:

  Character, data format - "stereo" (Stereo-seq) or "visium" (default:
  "stereo")

- data_type:

  Character, data type - "stRNA" (spatial) or "scRNA" (single-cell)
  (default: "stRNA")

- resolution_index:

  Numeric vector, clustering resolutions to test (default: seq(0.1, 1.3,
  0.2))

- runTSNE_index:

  Logical, whether to perform t-SNE (default: FALSE)

- assay_nm:

  Character, assay name (default: NULL, auto-detected)

## Value

A Seurat object with processed data and cluster annotations

## Examples

``` r
if (FALSE) { # \dontrun{
# Run standard pipeline on spatial data
seurat_obj <- Seurat_pipeline(
  seurat_obj = my_spatial_data,
  data_format = "stereo",
  resolution_index = seq(0.2, 1.0, 0.2)
)
} # }
```
