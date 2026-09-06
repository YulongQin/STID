# Convert STID object to Seurat object

Converts an STID object back to a standard Seurat object, preserving all
assay data and basic metadata but removing STID-specific analysis
structures.

## Usage

``` r
as.Seurat.STID(STID_obj)
```

## Arguments

- STID_obj:

  An STID object to be converted

## Value

A Seurat object containing the same assay data and basic metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert STID object back to Seurat
seurat_obj <- as.Seurat.STID(ist_obj)
} # }
```
