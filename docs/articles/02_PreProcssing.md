# 02_PreProcssing

## Overview

Preprocessing standardizes spatial transcriptomics data prior to
infection analysis. In ISTtools, this step ensures that input data are
normalized, dimensionally reduced, and annotated in a consistent format
compatible with downstream modules.

This vignette describes:

- Seurat-based preprocessing workflow
- Cell-type annotation using reference mapping

------------------------------------------------------------------------

## Setup

``` r
suppressMessages({
  library(tidyverse)
  library(ISTools)
})
```

------------------------------------------------------------------------

## Seurat preprocessing pipeline

ISTtools provides a wrapper function to streamline Seurat-based
preprocessing, including normalization, feature selection,
dimensionality reduction, and clustering.

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

### Description

The preprocessing pipeline performs:

- Data normalization and scaling
- Identification of highly variable features
- Principal component analysis (PCA)
- Graph-based clustering across multiple resolutions

The parameter `resolution_index` allows exploration of clustering
granularity, which may influence downstream spatial interpretation.

------------------------------------------------------------------------

## Cell-type annotation using SingleR

To facilitate biological interpretation, spatial spots can be annotated
using reference-based methods.

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

### Description

This function maps query spots to reference cell types using the
**SingleR** framework:

- `ref_obj`: reference dataset with known annotations
- `seurat_colnm`: clustering labels in query data
- `ref_colnm`: annotation labels in reference

The resulting annotations are stored in the Seurat object and can be
used for:

- Interpreting infection niches
- Comparing host cell-type composition
- Downstream niche-level analysis

------------------------------------------------------------------------

## Integration with downstream analysis

Preprocessed data serve as input for:

- Background correction
- Infection spot detection
- Spatial niche reconstruction

Consistent preprocessing ensures comparability across samples and
conditions.

------------------------------------------------------------------------

## Remarks

- Preprocessing parameters should be adapted to sequencing depth and
  platform
- Multiple clustering resolutions can be retained for flexibility
- Annotation quality depends on the relevance of the reference dataset

------------------------------------------------------------------------

## Session information

``` r
sessionInfo()
```
