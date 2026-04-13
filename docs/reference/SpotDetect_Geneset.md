# Detect infection spots based on gene set scores

This function calculates scores for gene sets using various methods and
identifies positive spots based on score thresholds. Supports multiple
scoring algorithms including AddModuleScore, AUCell, UCell, and basic
statistics.

## Usage

``` r
SpotDetect_Geneset(
  STID_obj = NULL,
  geneset_list = NULL,
  score_method = "AddModuleScore",
  n_iter = 5,
  nbin = 10,
  seed = 10,
  loop_id = "LoopAllSamp",
  assay_id = "Spatial",
  layer_id = "data",
  PosThres_prob = 0,
  PosThres_score = 0,
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
  dir_nm = "M1_SpotDetect_Geneset"
)
```

## Arguments

- STID_obj:

  An STID object containing spatial transcriptomics data

- geneset_list:

  A named list of gene sets for scoring

- score_method:

  Character, scoring method to use. Options: "AddModuleScore", "AUCell",
  "UCell", "MeanExp", "SumExp" (default: "AddModuleScore")

- n_iter:

  Numeric, number of iterations for AddModuleScore (default: 5)

- nbin:

  Numeric, number of bins for AddModuleScore (default: 10)

- seed:

  Numeric, random seed for reproducibility (default: 10)

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- assay_id:

  Character, name of the assay to use (default: "Spatial")

- layer_id:

  Character, name of the layer/data slot to use (default: "data")

- PosThres_prob:

  Numeric, probability threshold (0-1) for positive detection (default:
  0)

- PosThres_score:

  Numeric, absolute score threshold for positive detection (default: 0)

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

  Character, directory name for output (default:
  "M1_SpotDetect_Geneset")

## Value

Returns the modified STID object with added metadata containing gene set
scores and labels

## Examples

``` r
if (FALSE) { # \dontrun{
# Define gene sets
pathogen_genes <- list(
  "Virulence" = c("gene1", "gene2", "gene3"),
  "Toxin" = c("gene4", "gene5", "gene6")
)

# Detect spots based on gene set scores
STID_obj <- SpotDetect_Geneset(
  STID_obj = STID_object,
  geneset_list = pathogen_genes,
  score_method = "AddModuleScore",
  PosThres_prob = 0.95
)
} # }
```
