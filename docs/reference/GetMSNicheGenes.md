# Get MultiSampNiche genes

Retrieves gene-level data from MultiSampNiche objects for specified
multi-sample analyses.

## Usage

``` r
GetMSNicheGenes(STID_obj, loop_id, niche_key, ...)

# S3 method for class 'STID'
GetMSNicheGenes(
  STID_obj = NULL,
  loop_id = "LoopAllMulti",
  niche_key = NULL,
  ...
)
```

## Arguments

- STID_obj:

  An STID object

- loop_id:

  Character vector, multi-sample analysis identifiers

- niche_key:

  Character, niche key to retrieve

- ...:

  Additional arguments passed to methods

## Value

List of gene data frames, one per multi-sample analysis

## Examples

``` r
if (FALSE) { # \dontrun{
multi_genes <- GetMSNicheGenes(ist_obj,
                               loop_id = c("comparison1", "comparison2"),
                               niche_key = "niche_virulence")
} # }
```
