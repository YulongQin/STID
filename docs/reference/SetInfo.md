# Set information in STID object

Modifies specific information in the STID object's STID_info slot.

## Usage

``` r
SetInfo(STID_obj, info_key, sub_key, info_value, ...)

# S3 method for class 'STID'
SetInfo(
  STID_obj = NULL,
  info_key = NULL,
  sub_key = NULL,
  info_value = NULL,
  ...
)
```

## Arguments

- STID_obj:

  An STID object

- info_key:

  Character, main information category

- sub_key:

  Character, specific sub-key to set (optional)

- info_value:

  Value to set

- ...:

  Additional arguments passed to methods

## Value

Modified STID object

## Examples

``` r
if (FALSE) { # \dontrun{
# Update project description
STID_obj <- SetInfo(ist_obj,
                   info_key = "project_info",
                   sub_key = "description",
                   info_value = "Updated project description")
} # }
```
