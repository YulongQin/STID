# Get coordinate-related column names

Convenience function to retrieve all coordinate-related column names
from an STID object.

## Usage

``` r
GetCoordInfo(STID_obj)

# S3 method for class 'STID'
GetCoordInfo(STID_obj = NULL)
```

## Arguments

- STID_obj:

  An STID object

## Value

Named vector of coordinate column names

## Examples

``` r
if (FALSE) { # \dontrun{
coord_info <- GetCoordInfo(ist_obj)
} # }
```
