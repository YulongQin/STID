# Remove columns from existing metadata

Removes specified columns from an existing metadata data frame in the
STID object.

## Usage

``` r
RemoveMetaColumn(STID_obj, ...)

# S3 method for class 'STID'
RemoveMetaColumn(STID_obj = NULL, meta_key = NULL, remove_colnm = NULL, ...)
```

## Arguments

- STID_obj:

  An STID object

- ...:

  Additional arguments passed to methods

- meta_key:

  Character, metadata key to modify

- remove_colnm:

  Character vector, column names to remove

## Value

Modified STID object with removed metadata columns

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove temporary columns
STID_obj <- RemoveMetaColumn(ist_obj,
                            meta_key = "raw",
                            remove_colnm = c("temp_score", "temp_label"))
} # }
```
