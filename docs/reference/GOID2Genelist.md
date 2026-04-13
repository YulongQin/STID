# Convert GO IDs to Gene Lists

Retrieves gene lists associated with specified Gene Ontology (GO) terms.
This function extracts genes annotated to given GO IDs from the
appropriate organism database.

## Usage

``` r
GOID2Genelist(
  GOID = NULL,
  STID_obj = NULL,
  host_org = NULL,
  ont = "BP",
  keyType = "SYMBOL"
)
```

## Arguments

- GOID:

  Character vector, GO term identifiers to query

- STID_obj:

  STID object (optional, used to extract host organism)

- host_org:

  Character, host organism - "human" or "mouse" (if NULL, extracted from
  STID_obj)

- ont:

  Character, GO ontology type - "BP", "CC", "MF", or "ALL" (default:
  "BP")

- keyType:

  Character, gene identifier type - "SYMBOL" or "ENTREZID" (default:
  "SYMBOL")

## Value

A data frame where columns are GO terms (with names as descriptions) and
rows are genes. Contains gene lists for each requested GO term.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get gene list for specific GO terms
go_genes <- GOID2Genelist(
  GOID = c("GO:0006955", "GO:0002376"),
  host_org = "human",
  ont = "BP"
)

# View the gene list for the first GO term
head(go_genes[,1])
} # }
```
