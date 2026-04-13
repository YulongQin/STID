# Add annotations to SingleSampNiche genes

Adds new columns to the gene data in SingleSampNiche objects.

## Usage

``` r
AddSSNicheGenes(STID_obj, ...)

# S3 method for class 'STID'
AddSSNicheGenes(
  STID_obj = NULL,
  gene = NULL,
  label = NULL,
  add_colnm = NULL,
  loop_id = "LoopAllSamp",
  niche_key = NULL,
  ...
)
```

## Arguments

- STID_obj:

  An STID object

- ...:

  Additional arguments passed to methods

- gene:

  Character vector, gene names to annotate

- label:

  Character vector, labels to assign to genes (optional)

- add_colnm:

  Character, name of new column to add

- loop_id:

  Character vector, sample identifiers

- niche_key:

  Character, niche key to modify

## Value

Modified STID object with updated niche gene data

## Examples

``` r
if (FALSE) { # \dontrun{
# Add gene family annotations
STID_obj <- AddSSNicheGenes(ist_obj,
                           gene = c("gene1", "gene2", "gene3"),
                           label = c("toxin", "toxin", "adhesin"),
                           add_colnm = "gene_family",
                           niche_key = "niche_virulence")
} # }
```
