# Annotate Cell Types Using SingleR

Performs automated cell type annotation using the SingleR package. Uses
a reference dataset (e.g., from celldex) to annotate cells based on gene
expression profiles.

## Usage

``` r
anno_SingleR(
  seurat_obj = NULL,
  ref_obj = NULL,
  seurat_colnm = NULL,
  ref_colnm = NULL,
  species_index = NULL,
  assay_nm = "Spatial",
  layer_nm = "data"
)
```

## Arguments

- seurat_obj:

  A Seurat object to annotate

- ref_obj:

  Reference dataset (matrix or SummarizedExperiment) with cell type
  labels

- seurat_colnm:

  Character, column name in Seurat metadata containing cluster IDs for
  aggregation

- ref_colnm:

  Character, column name in reference object containing cell type labels

- species_index:

  Character, species for built-in references (currently unused)

- assay_nm:

  Character, assay to use (default: "Spatial")

- layer_nm:

  Character, layer/slot to use (default: "data")

## Value

A Seurat object with added 'anno_SingleR' column in metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# Annotate using built-in reference
library(SingleR)
ref <- HumanPrimaryCellAtlasData()

seurat_obj <- anno_SingleR(
  seurat_obj = my_seurat,
  ref_obj = ref,
  seurat_colnm = "seurat_clusters",
  ref_colnm = "label.main"
)
} # }
```
