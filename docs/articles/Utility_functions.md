# Utility functions

## Overview

This vignette introduces four utility functions used for environment
configuration, cross-format data conversion, and annotation parsing in
STID workflows. `configure_conda()` prepares the Python environment
required by Python-backed conversion steps.
[`h5ad2rds()`](https://yulongqin.github.io/STID/reference/h5ad2rds.md)
and
[`rds2h5ad()`](https://yulongqin.github.io/STID/reference/rds2h5ad.md)
convert data between AnnData `h5ad` files and Seurat-compatible R
objects.
[`parse_gtf()`](https://yulongqin.github.io/STID/reference/parse_gtf.md)
extracts structured gene annotation information from GTF files.

## Function summary

| Function | Purpose | Typical input | Typical output | Notes |
|----|----|----|----|----|
| `configure_conda()` | Configure the Python environment through `reticulate` | Conda environment name, or default Miniconda setup | Activated Python/Conda environment | Run before Python-backed conversion workflows |
| [`h5ad2rds()`](https://yulongqin.github.io/STID/reference/h5ad2rds.md) | Convert AnnData `h5ad` files to Seurat objects | `.h5ad` file path | Seurat object and/or `.rds` file | Supports spatial and single-cell modes, reductions, coordinates, and optional SCT conversion |
| [`rds2h5ad()`](https://yulongqin.github.io/STID/reference/rds2h5ad.md) | Convert Seurat objects to AnnData-compatible `h5ad` files | Seurat object | `.h5ad` file | Current helper supports the SeuratDisk conversion path |
| [`parse_gtf()`](https://yulongqin.github.io/STID/reference/parse_gtf.md) | Parse GTF annotation files into tidy tables | `.gtf` file path | Data frame containing genomic features and attributes | Use `fil_label` to filter feature types such as `gene` |

## Prerequisites

``` r
library(tidyverse)
library(Seurat)
library(STID)
```

## `configure_conda()`

`configure_conda()` configures Python through `reticulate`. If no Conda
environment name is supplied, it attempts to use or install the default
Miniconda environment. If a Conda environment name is supplied, it
checks whether that environment exists before activating it.

``` r
# Use default Miniconda
configure_conda()

# Use an existing Conda environment
configure_conda(conda_nm = "scanpy_env")
```

Run this configuration before using Python-backed conversion workflows
such as
[`h5ad2rds()`](https://yulongqin.github.io/STID/reference/h5ad2rds.md)
and
[`rds2h5ad()`](https://yulongqin.github.io/STID/reference/rds2h5ad.md).

## `h5ad2rds()`

[`h5ad2rds()`](https://yulongqin.github.io/STID/reference/h5ad2rds.md)
converts an AnnData `h5ad` file into a Seurat object. It supports
spatial and single-cell modes, optional reductions, optional spatial
coordinates, and optional SCT conversion.

``` r
stRNA <- h5ad2rds(
  file_path = "data/sample_spatial.h5ad",
  data_type = "stRNA",
  convert_mode = "scanpy",
  assay_id = "Spatial",
  X_index = "rawX",
  binsize = 100,
  SCT_index = FALSE,
  reduction_index = TRUE,
  image_index = TRUE,
  return_object = TRUE,
  grp_nm = "sample_spatial",
  dir_nm = "M0_h5ad2rds"
)
```

## `rds2h5ad()`

[`rds2h5ad()`](https://yulongqin.github.io/STID/reference/rds2h5ad.md)
converts a Seurat object into an AnnData-compatible `h5ad` file using
the SeuratDisk path.

``` r
rds2h5ad(
  seurat_obj = stRNA,
  data_type = "stRNA",
  convert_mode = "seurat",
  assay_id = "Spatial",
  grp_nm = "sample_spatial",
  dir_nm = "M0_rds2h5ad"
)
```

> **Note:** The Scanpy mode of
> [`rds2h5ad()`](https://yulongqin.github.io/STID/reference/rds2h5ad.md)
> is not implemented in the current helper. Use
> `convert_mode = "seurat"` for this direction.

## `parse_gtf()`

[`parse_gtf()`](https://yulongqin.github.io/STID/reference/parse_gtf.md)
reads a GTF file and returns a tidy data frame. The `fil_label` argument
filters the feature column, for example to keep only genes.

``` r
gtf_gene <- parse_gtf(
  gtf_file = "reference/annotation.gtf",
  fil_label = "gene"
)

head(gtf_gene)
```

``` r
gene_table <- parse_gtf("reference/annotation.gtf", fil_label = "gene") %>%
  dplyr::select(seqname, start, end, strand, gene_id, gene_name, gene_biotype)
```

## Session information

``` r
sessionInfo()
#> R version 4.2.0 (2022-04-22 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 22000)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
#> [2] LC_CTYPE=Chinese (Simplified)_China.utf8   
#> [3] LC_MONETARY=Chinese (Simplified)_China.utf8
#> [4] LC_NUMERIC=C                               
#> [5] LC_TIME=Chinese (Simplified)_China.utf8    
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.35     R6_2.6.1          jsonlite_1.8.8    lifecycle_1.0.5  
#>  [5] evaluate_1.0.1    cachem_1.1.0      rlang_1.1.7       cli_3.6.5        
#>  [9] rstudioapi_0.15.0 fs_1.6.3          jquerylib_0.1.4   bslib_0.8.0      
#> [13] ragg_1.3.0        rmarkdown_2.29    pkgdown_2.2.0     textshaping_0.3.6
#> [17] desc_1.4.3        tools_4.2.0       htmlwidgets_1.6.4 yaml_2.3.10      
#> [21] xfun_0.49         fastmap_1.2.0     compiler_4.2.0    systemfonts_1.0.4
#> [25] htmltools_0.5.8.1 knitr_1.49        sass_0.4.9
```
