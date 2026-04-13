# Add columns to SingleSampNiche cells

Adds new columns to the cell data in SingleSampNiche objects.

## Usage

``` r
AddSSNicheCells(STID_obj, loop_id, meta_key, select_colnm, niche_key, ...)

# S3 method for class 'STID'
AddSSNicheCells(
  STID_obj = NULL,
  loop_id = "LoopAllSamp",
  meta_key = NULL,
  select_colnm = NULL,
  niche_key = NULL,
  ...
)
```

## Arguments

- STID_obj:

  An STID object

- loop_id:

  Character vector, sample identifiers

- meta_key:

  Character, metadata key containing new data

- select_colnm:

  Character vector, columns to add

- niche_key:

  Character, niche key to modify

- ...:

  Additional arguments passed to methods

## Value

Modified STID object with updated niche cell data

## Examples

``` r
if (FALSE) { # \dontrun{
# Add new annotation columns to niche cells
STID_obj <- AddSSNicheCells(ist_obj,
                           loop_id = "LoopAllSamp",
                           meta_key = "custom_annotation",
                           select_colnm = c("cell_type", "confidence"),
                           niche_key = "niche_virulence")
} # }
```
