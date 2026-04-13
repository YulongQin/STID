# Calculate summary statistics for selected genes

This function computes summary statistics (nCount and nFeature) for a
set of selected genes, useful for creating pathogen-specific metrics.

## Usage

``` r
GetGeneStat(
  STID_obj = NULL,
  pattern = NULL,
  features = NULL,
  prefix = NULL,
  func = "sum",
  assay_id = "Spatial",
  layer_id = "counts"
)
```

## Arguments

- STID_obj:

  A Seurat or STID object

- pattern:

  Character, regular expression pattern to select genes (default: NULL)

- features:

  Character vector of specific gene names to include (default: NULL)

- prefix:

  Character, prefix for output column names

- func:

  Character, function to apply for nCount calculation. Options: "sum",
  "mean", "median", "max", "min" (default: "sum")

- assay_id:

  Character, name of the assay to use (default: "Spatial")

- layer_id:

  Character, name of the layer/data slot to use (default: "counts")

## Value

A data frame with two columns:

- prefix_nCount(func):

  Summary statistic for total counts

- prefix_nFeature(sum):

  Number of detected features

## Examples

``` r
if (FALSE) { # \dontrun{
# Calculate pathogen gene statistics
pathogen_stats <- GetGeneStat(
  STID_obj = STID_object,
  pattern = "^MTB",
  prefix = "Pathogen",
  func = "sum"
)

# Add to metadata
STID_object <- AddMetaData(STID_object, pathogen_stats)
} # }
```
