# Retrieve sample-specific thresholds from previous analyses

This function retrieves stored threshold information from previous
SpotDetect analyses (both gene and gene set based)

## Usage

``` r
GetSapThreshold(STID_obj, meta_key, ...)

# S3 method for class 'STID'
GetSapThreshold(STID_obj = NULL, meta_key = NULL, ...)
```

## Arguments

- STID_obj:

  An STID object

- meta_key:

  Character vector of metadata keys to retrieve thresholds for

- ...:

  Additional arguments passed to methods

## Value

A list containing threshold information for each requested metadata key,
with components 'stat' and 'threshold'

## Examples

``` r
if (FALSE) { # \dontrun{
# Retrieve thresholds from previous analyses
thresholds <- GetSapThreshold(
  STID_obj = STID_object,
  meta_key = c("M1_SpotDetect_Gene_20240101", "M1_SpotDetect_Geneset_20240101")
)
} # }
```
