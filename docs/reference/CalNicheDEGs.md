# Calculate Niche Differential Expression Genes

Identifies differentially expressed genes between niche and bystander
cells, both globally and per cell type. Generates volcano plots and
faceted volcano plots for visualization.

## Usage

``` r
CalNicheDEGs(
  STID_obj = NULL,
  samp_mode = "SS",
  loop_id = "LoopAllSamp",
  samp_grp_index = FALSE,
  meta_key = NULL,
  niche_key = NULL,
  group_by = NULL,
  group_value = NULL,
  assay_id = "Spatial",
  test_method = "wilcox",
  logfc_thres = 1,
  min_pct = 0.01,
  padj_thres = 0.05,
  adjust_method = "BH",
  topGeneN = 3,
  col = COLOR_List[["PALETTE_WHITE_BG"]],
  remove_genes = NULL,
  return_data = TRUE,
  grp_nm = NULL,
  dir_nm = "M3_CalNicheDEGs",
  ...
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- samp_mode:

  Character, sample type - "SS" (single sample) or "MS" (multi-sample)
  (default: "SS")

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- samp_grp_index:

  Logical, whether to group by sample groups in MS mode (default: FALSE)

- meta_key:

  Character, metadata key for MS mode when niche_key is NULL

- niche_key:

  Character, niche key to analyze

- group_by:

  Character, column name for cell type grouping

- group_value:

  Character vector, specific cell types to analyze

- assay_id:

  Character, assay name (default: "Spatial")

- test_method:

  Character, statistical test to use (default: "wilcox")

- logfc_thres:

  Numeric, log2 fold change threshold (default: 1)

- min_pct:

  Numeric, minimum percentage of cells expressing gene (default: 0.01)

- padj_thres:

  Numeric, adjusted p-value threshold (default: 0.05)

- adjust_method:

  Character, p-value adjustment method (default: "BH")

- topGeneN:

  Integer, number of top genes to label (default: 3)

- col:

  Color palette for visualization (default:
  COLOR_List\$PALETTE_WHITE_BG)

- remove_genes:

  Character vector, genes to exclude from analysis

- return_data:

  Logical, whether to return results list (default: TRUE)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M3_CalNicheDEGs")

- ...:

  Additional parameters for FindAllMarkers and FindMarkers functions

## Value

If return_data = TRUE, returns a list of DEG results per sample

## Examples

``` r
if (FALSE) { # \dontrun{
# Find DEGs between niche and bystander cells
results <- CalNicheDEGs(
  STID_obj = STID_obj,
  niche_key = "niche_virulence",
  group_by = "cell_type",
  logfc_thres = 0.5,
  padj_thres = 0.05,
  topGeneN = 5
)
} # }
```
