# Get SingleSampNiche information

Retrieves niche information from SingleSampNiche objects for specified
samples.

## Usage

``` r
GetSSNicheInfo(STID_obj, loop_id, niche_key, ...)

# S3 method for class 'STID'
GetSSNicheInfo(STID_obj = NULL, loop_id = "LoopAllSamp", niche_key = NULL, ...)
```

## Arguments

- STID_obj:

  An STID object

- loop_id:

  Character vector, sample identifiers

- niche_key:

  Character, niche key to retrieve

- ...:

  Additional arguments passed to methods

## Value

List of niche information data frames, one per sample

## Examples

``` r
if (FALSE) { # \dontrun{
niche_info <- GetSSNicheInfo(ist_obj,
                             loop_id = c("sample1", "sample2"),
                             niche_key = "niche_virulence")
} # }
```
