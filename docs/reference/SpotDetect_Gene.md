# Detect infection spots based on individual gene expression

This function identifies positive spots based on expression thresholds
for individual genes. It can apply spatial smoothing and generates
visualization of gene expression patterns.

## Usage

``` r
SpotDetect_Gene(
  STID_obj = NULL,
  features = NULL,
  feature_colnm = NULL,
  loop_id = "LoopAllSamp",
  assay_id = "Spatial",
  layer_id = "counts",
  PosThres_prob = 0,
  PosThres_count = 1,
  PosThres_operator = c(">=", ">"),
  col = COLOR_DIS_CON,
  pt_size = 0.5,
  vmin = NULL,
  vmax = "p99",
  black_bg = FALSE,
  plot_method = "merge",
  blur_method = NULL,
  blur_n = 1,
  blur_sigma = 0.25,
  description = NULL,
  grp_nm = NULL,
  dir_nm = "M1_SpotDetect_Gene"
)
```

## Arguments

- STID_obj:

  An STID object containing spatial transcriptomics data

- features:

  Character vector of gene names to analyze

- feature_colnm:

  Character, column name in metadata containing pre-computed features

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- assay_id:

  Character, name of the assay to use (default: "Spatial")

- layer_id:

  Character, name of the layer/data slot to use (default: "counts")

- PosThres_prob:

  Numeric, probability threshold (0-1) for positive detection (default:
  0)

- PosThres_count:

  Numeric, absolute count threshold for positive detection (default: 1)

- PosThres_operator:

  Character, operator for threshold comparison - "\>=" or "\>" (default:
  "\>=")

- col:

  Color palette for visualization (default: COLOR_DIS_CON)

- pt_size:

  Numeric, point size in spatial plots (default: 0.5)

- vmin:

  Numeric or character, minimum value for color scale (default: NULL)

- vmax:

  Numeric or character, maximum value for color scale (default: "p99")

- black_bg:

  Logical, whether to use black background in plots (default: FALSE)

- plot_method:

  Character, plotting mode - "merge" or "single" (default: "merge")

- blur_method:

  Character, spatial smoothing method - "isoblur" or "medianblur"
  (default: NULL)

- blur_n:

  Numeric, number of iterations for median blur (default: 1)

- blur_sigma:

  Numeric, sigma parameter for isoblur (default: 0.25)

- description:

  Character, description of the analysis (default: NULL)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M1_SpotDetect_Gene")

## Value

Returns the modified STID object with added metadata containing gene
expression labels and thresholds

## Examples

``` r
if (FALSE) { # \dontrun{
# Detect spots based on pathogen gene expression
STID_obj <- SpotDetect_Gene(
  STID_obj = STID_object,
  features = c("geneA", "geneB", "geneC"),
  PosThres_prob = 0.95,
  PosThres_count = 1,
  plot_method = "merge"
)
} # }
```
