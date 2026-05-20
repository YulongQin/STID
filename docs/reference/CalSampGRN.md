# Calculate Niche Gene Regulatory Networks using NicheNet

Predicts ligand-target gene regulatory networks within niches using the
NicheNet framework. Identifies potential signaling from sender to
receiver cells.

## Usage

``` r
CalSampGRN(
  STID_obj = NULL,
  loop_id = "LoopAllSamp",
  niche_key = NULL,
  group_by = NULL,
  ref_data = NULL,
  sender_celltypes = NULL,
  receiver_celltypes = NULL,
  target_features = NULL,
  expression_pct = 0.05,
  top_ligand_num = 10,
  top_target_num = 5,
  remove_genes = NULL,
  return_data = TRUE,
  grp_nm = NULL,
  dir_nm = "M3_CalSampGRN"
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- niche_key:

  Character, niche key to analyze

- group_by:

  Character, column name for cell type grouping

- ref_data:

  List, NicheNet reference data (default: NULL, downloads automatically)

- sender_celltypes:

  Character vector, sender cell types

- receiver_celltypes:

  Character vector, receiver cell types

- target_features:

  Character vector, target genes of interest

- expression_pct:

  Numeric, expression percentage threshold (default: 0.05)

- top_ligand_num:

  Integer, number of top ligands to prioritize (default: 10)

- top_target_num:

  Integer, number of top targets per ligand (default: 5)

- remove_genes:

  Character vector, genes to exclude

- return_data:

  Logical, whether to return results list (default: TRUE)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M3_CalSampGRN")

## Value

If return_data = TRUE, returns a list of NicheNet results per sample

## Examples

``` r
if (FALSE) { # \dontrun{
# Predict ligand-target networks
results <- CalSampGRN(
  STID_obj = STID_obj,
  niche_key = "niche_virulence",
  sender_celltypes = "Epithelial",
  receiver_celltypes = "Immune",
  target_features = c("gene1", "gene2", "gene3")
)
} # }
```
