# Calculate Niche Gene Expression Correlations

Computes correlation between pathogen and host gene expression within
niches. Supports both feature-based and gene-based correlation analysis.

## Usage

``` r
CalNicheGeneCor(
  STID_obj = NULL,
  loop_id = "LoopAllSamp",
  niche_key = NULL,
  meta_key = NULL,
  p.features = NULL,
  h.features = NULL,
  p.feature_colnm = NULL,
  h.feature_colnm = NULL,
  method = "spearman",
  assay_id = "Spatial",
  layer_id = "data",
  return_data = TRUE,
  grp_nm = NULL,
  dir_nm = "M3_CalNicheGeneCor"
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- niche_key:

  Character, niche key to analyze

- meta_key:

  Character, metadata key for analysis when niche_key is NULL

- p.features:

  Character vector, pathogen gene names

- h.features:

  Character vector, host gene names

- p.feature_colnm:

  Character vector, pathogen feature column names in metadata

- h.feature_colnm:

  Character vector, host feature column names in metadata

- method:

  Character, correlation method - "spearman" or "pearson" (default:
  "spearman")

- assay_id:

  Character, assay name (default: "Spatial")

- layer_id:

  Character, layer name (default: "data")

- return_data:

  Logical, whether to return results list (default: TRUE)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M3_CalNicheGeneCor")

## Value

If return_data = TRUE, returns a list of correlation results per sample

## Examples

``` r
if (FALSE) { # \dontrun{
# Correlate pathogen and host gene expression
results <- CalNicheGeneCor(
  STID_obj = STID_obj,
  niche_key = "niche_virulence",
  p.features = c("pathogen_gene1", "pathogen_gene2"),
  h.features = c("host_gene1", "host_gene2"),
  method = "spearman"
)
} # }
```
