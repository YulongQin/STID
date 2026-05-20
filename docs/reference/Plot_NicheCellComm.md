# Plot Niche Cell Communication Results

Visualizes CellChat analysis results with various plot types including
circle plots, heatmaps, bubble plots, and spatial signaling plots.

## Usage

``` r
Plot_NicheCellComm(
  STID_obj = NULL,
  samp_mode = "SS",
  loop_id = "LoopAllSamp",
  CellComm_data = NULL,
  sources.use = NULL,
  targets.use = NULL,
  signaling = NULL,
  pairLR.use = NULL,
  col = NULL
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- samp_mode:

  Character, sample type - "SS" or "MS" (default: "SS")

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- CellComm_data:

  List, CellChat results from CalSampCellComm

- sources.use:

  Character vector, source cell types to include

- targets.use:

  Character vector, target cell types to include

- signaling:

  Character vector, signaling pathways to plot

- pairLR.use:

  Data frame, specific ligand-receptor pairs to plot

- col:

  Named vector, colors for cell types

## Value

NULL (invisible), generates plots

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot CellChat results
Plot_NicheCellComm(
  STID_obj = STID_obj,
  CellComm_data = cellcomm_results,
  signaling = c("TGFb", "WNT"),
  sources.use = "Epithelial",
  targets.use = "Immune"
)
} # }
```
