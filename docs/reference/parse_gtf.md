# Parse GTF File into Tidy Data Frame

Parses a GTF (Gene Transfer Format) file and converts it into a tidy
data frame with key-value attributes properly separated into columns.

## Usage

``` r
parse_gtf(gtf_file = NULL, fil_label = "gene")
```

## Arguments

- gtf_file:

  Character, path to the GTF file

- fil_label:

  Character, feature type to filter (e.g., "gene", "exon") (default:
  "gene")

## Value

A data frame with GTF columns (seqname, source, feature, start, end,
score, strand, frame) plus additional columns for each attribute

## Examples

``` r
if (FALSE) { # \dontrun{
# Parse a GTF file and keep only gene features
gtf_data <- parse_gtf(
  gtf_file = "path/to/annotation.gtf",
  fil_label = "gene"
)

# View the parsed data
head(gtf_data)
} # }
```
