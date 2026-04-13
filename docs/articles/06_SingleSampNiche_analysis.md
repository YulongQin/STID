# 06_SingleSampNiche_analysis

## Overview

Within a single sample, infection niches define local microenvironments
where host–pathogen interactions occur. ISTtools enables quantitative
characterization of these niches by modeling **spatial gradients
relative to infection regions**.

------------------------------------------------------------------------

## Setup

``` r
suppressMessages({
  library(tidyverse)
  library(ISTools)
})
```

------------------------------------------------------------------------

## Distance-dependent expression profiling

To quantify spatial organization, gene expression is modeled as a
function of distance to infection niches.

``` r
Plot_DistLine_Exp(
  IST_obj = IST_obj_detect,
  features = c("NS5","Ccl2"),
  feature_colnm = "all_gene_nFeature(sum)",
  loop_id = "D5_1",
  col = c("#F81B02FF","#3B95C4FF","#F81B02FF"),
  meta_key = list(c(
    "M1_SpotDetect_Gene_JEV_correct_before_all_gene_white",
    "M2_NicheDetect_STS_STS_JEV_microbe_region"
  ))
)
```

This analysis captures continuous spatial transitions from infected to
non-infected regions.

------------------------------------------------------------------------

## Niche expansion

To further model local microenvironments, infection niches can be
expanded radially.

``` r
IST_obj_expand <- NicheExpand(
  IST_obj_detect,
  meta_key = "M2_NicheDetect_STS_STS_JEV_microbe_region",
  pos_colnm = "ROI_label",
  center_colnm = "ROI_center",
  expand_dist = 30
)
```

Expanded regions approximate the spatial influence range of infection.

------------------------------------------------------------------------

## Remarks

- Distance-based modeling reveals gradients in host response
- Niche expansion provides a controllable spatial context definition
- These analyses enable inference of local microenvironment structure

------------------------------------------------------------------------

## Session information

``` r
sessionInfo()
```
