# Standard Seurat Analysis Pipeline for Spatial Transcriptomics

Performs a complete Seurat analysis pipeline for spatial transcriptomics
data, including normalization, feature selection, scaling, PCA,
clustering, and UMAP/t-SNE dimensionality reduction. Also generates
spatial feature plots.

## Usage

``` r
Seurat_pipeline(
  seurat_obj = NULL,
  data_type = "stRNA",
  resolution_index = seq(0.1, 1.3, 0.2),
  runTSNE_index = FALSE,
  assay_nm = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object containing spatial transcriptomics data

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
stRNA <- Seurat_pipeline(
  seurat_obj = stRNA,
  resolution_index = seq(0.2, 1.0, 0.2)
)
} # }
```
