# Identify top expressed genes in the dataset

This function identifies the top expressed genes either globally or
grouped by sample and/or cell type. Useful for finding marker genes or
highly expressed pathogen genes.

## Usage

``` r
GetTopGenes(
  STID_obj,
  top_n = 10,
  pattern = NULL,
  grp_by_samp = FALSE,
  grp_by_celltype = FALSE,
  assay_id = "Spatial",
  layer_id = "counts"
)
```

## Arguments

- STID_obj:

  A Seurat or STID object

- top_n:

  Numeric, number of top genes to return (default: 10)

- pattern:

  Character, regular expression pattern to filter genes (default: NULL)

- grp_by_samp:

  Logical, whether to group by sample (default: FALSE)

- grp_by_celltype:

  Logical, whether to group by cell type (default: FALSE)

- assay_id:

  Character, name of the assay to use (default: "Spatial")

- layer_id:

  Character, name of the layer/data slot to use (default: "counts")

## Value

If no grouping, returns a character vector of top genes. If grouping,
returns a named list of top genes per group.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get top 20 genes globally
top_genes <- GetTopGenes(STID_obj = STID_object, top_n = 20)

# Get top 10 pathogen genes per sample
pathogen_genes <- GetTopGenes(
  STID_obj = STID_object,
  top_n = 10,
  pattern = "^VP",
  grp_by_samp = TRUE
)
} # }
```
