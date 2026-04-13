# Remove metadata from STID object

Removes specified metadata from the STID object's meta_data_record.

## Usage

``` r
RemoveMetaData(STID_obj, meta_key, ...)

# S3 method for class 'STID'
RemoveMetaData(STID_obj = NULL, meta_key = NULL, ...)
```

## Arguments

- STID_obj:

  An STID object

- meta_key:

  Character vector, metadata keys to remove

- ...:

  Additional arguments passed to methods

## Value

Modified STID object with removed metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove temporary metadata
STID_obj <- RemoveMetaData(ist_obj, meta_key = "temp_analysis")
} # }
```
