# Perform Gene Enrichment Analysis (GO/KEGG or GSEA)

Performs gene enrichment analysis using either traditional
over-representation analysis (GO/KEGG) or Gene Set Enrichment Analysis
(GSEA). Supports both human and mouse organisms.

## Usage

``` r
GeneEnrichment(
  STID_obj = NULL,
  up_gene = NULL,
  down_gene = NULL,
  DEGs = NULL,
  host_org = NULL,
  method = "GO_KEGG",
  plot_pathway_num = 12,
  internal = FALSE,
  col = viridis(100, option = "H")[15:85] %>% rev(),
  return_data = FALSE,
  grp_nm = NULL,
  dir_nm = "M3_GeneEnrichment",
  ...
)
```

## Arguments

- STID_obj:

  STID object (optional, used to extract host organism)

- up_gene:

  Character vector, upregulated genes for GO/KEGG analysis

- down_gene:

  Character vector, downregulated genes for GO/KEGG analysis

- DEGs:

  Data frame, differential expression results with columns 'SYMBOL' and
  'LogFC' for GSEA analysis

- host_org:

  Character, host organism - "human" or "mouse" (if NULL, extracted from
  STID_obj)

- method:

  Character, method to use - "GO_KEGG" or "GSEA_GO_KEGG" (default:
  "GO_KEGG")

- plot_pathway_num:

  Integer, number of top pathways to display in plots (default: 12)

- internal:

  Logical, whether to use internal KEGG database (default: FALSE)

- col:

  Color palette for plots (default:
  `viridis(100, option = "H")[15:85] %>% rev()`)

- return_data:

  Logical, whether to return enrichment results (default: FALSE)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M3_GeneEnrichment")

- ...:

  Additional arguments passed to internal functions

## Value

If return_data = TRUE, returns list of enrichment results; otherwise
NULL

## Examples

``` r
if (FALSE) { # \dontrun{
# GO/KEGG enrichment with up and down regulated genes
GeneEnrichment(
  up_gene = upregulated_genes,
  down_gene = downregulated_genes,
  host_org = "human",
  plot_pathway_num = 10
)

# GSEA enrichment with DEGs
GeneEnrichment(
  DEGs = deg_results,
  host_org = "mouse",
  method = "GSEA_GO_KEGG"
)
} # }
```
