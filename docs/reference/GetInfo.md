# Get information from STID object

Retrieves specific information from the STID object's STID_info slot.

## Usage

``` r
GetInfo(STID_obj, info_key, ...)

# S3 method for class 'STID'
GetInfo(STID_obj = NULL, info_key = NULL, sub_key = NULL, ...)
```

## Arguments

- STID_obj:

  An STID object

- info_key:

  Character, main information category (e.g., "samp_info", "data_info")

- ...:

  Additional arguments passed to methods

- sub_key:

  Character vector, specific sub-keys to retrieve (optional)

## Value

Requested information from the STID object

## Examples

``` r
if (FALSE) { # \dontrun{
# Get sample column name
samp_colnm <- GetInfo(ist_obj, info_key = "data_info", sub_key = "samp_colnm")

# Get all data information
data_info <- GetInfo(ist_obj, info_key = "data_info")
} # }
```
