# Create a MultiSampNiche object for cross-sample comparison

Creates and stores a MultiSampNiche object within an STID object,
combining niche analysis results from multiple samples for comparative
analysis.

## Usage

``` r
CreateMultiSampNiche(
  STID_obj = NULL,
  multi_id = NULL,
  loop_id = "LoopAllSamp",
  compare_mode = NULL,
  niche_key = NULL,
  description = NULL
)

# S3 method for class 'STID'
CreateMultiSampNiche(
  STID_obj = NULL,
  multi_id = NULL,
  loop_id = "LoopAllSamp",
  compare_mode = NULL,
  niche_key = NULL,
  description = NULL
)

CreateMSNiche(...)
```

## Arguments

- STID_obj:

  An STID object containing SingleSampNiche objects

- multi_id:

  Character, unique identifier for the multi-sample analysis

- loop_id:

  Character, sample grouping identifier (default: "LoopAllSamp")

- compare_mode:

  Character, comparison mode - "Comparative" or "Temporal"

- niche_key:

  Character, niche key to combine across samples

- description:

  Character, description of the multi-sample analysis

- ...:

  Additional arguments passed to methods

## Value

Modified STID object with added MultiSampNiche information

## Examples

``` r
if (FALSE) { # \dontrun{
# Create MultiSampNiche object for comparative analysis
STID_obj <- CreateMultiSampNiche(
  STID_obj = ist_obj,
  multi_id = "comparison_infected_vs_control",
  loop_id = c("infected1", "infected2", "control1", "control2"),
  compare_mode = "Comparative",
  niche_key = "niche_virulence"
)
} # }
```
