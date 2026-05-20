# Compare two niche detection results

Compares two sets of niche/ROI detection results using Venn diagrams and
distance analysis. Identifies overlapping spots and calculates distances
between niche centers.

## Usage

``` r
CompareNiche(STID_obj = NULL, meta_key1 = NULL, meta_key2 = NULL, bins = 15)
```

## Arguments

- STID_obj:

  An STID object containing niche detection results

- meta_key1:

  Character, first metadata key for comparison

- meta_key2:

  Character, second metadata key for comparison

- bins:

  Integer, number of bins for distance distribution plots (default: 15)

## Value

NULL (invisible), generates plots and prints comparison statistics

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare two niche detection methods
CompareNiche(
  STID_obj = ist_object,
  meta_key1 = "M2_NicheDetect_Lasso_20240101",
  meta_key2 = "M2_NicheDetect_STS_20240101"
)
} # }
```
