# Add information to STID object

Appends information to existing slots in the STID_info, primarily for
adding comments and notes.

## Usage

``` r
AddInfo(STID_obj, info_key, sub_key, info_value, ...)

# S3 method for class 'STID'
AddInfo(
  STID_obj = NULL,
  info_key = "comment_info",
  sub_key = NULL,
  info_value = NULL,
  ...
)
```

## Arguments

- STID_obj:

  An STID object

- info_key:

  Character, information category (default: "comment_info")

- sub_key:

  Character, specific sub-key to add to

- info_value:

  Value to append

- ...:

  Additional arguments passed to methods

## Value

Modified STID object

## Examples

``` r
if (FALSE) { # \dontrun{
# Add a comment about analysis parameters
STID_obj <- AddInfo(ist_obj,
                   info_key = "comment_info",
                   sub_key = "parameters",
                   info_value = "Used DBSCAN with eps=5, minPts=10")
} # }
```
