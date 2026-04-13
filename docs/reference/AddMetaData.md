# Add metadata to STID object

Adds a new metadata data frame to the STID object's meta_data_record and
updates the tracking information.

## Usage

``` r
AddMetaData(STID_obj, ...)

# S3 method for class 'STID'
AddMetaData(
  STID_obj = NULL,
  meta_key = NULL,
  add_data = NULL,
  dir_nm = NA,
  grp_nm = NA,
  asso_key = NULL,
  description = NULL,
  ...
)
```

## Arguments

- STID_obj:

  An STID object

- ...:

  Additional arguments passed to methods

- meta_key:

  Character, unique identifier for the metadata

- add_data:

  Data frame, metadata to add

- dir_nm:

  Character, directory name for output (optional)

- grp_nm:

  Character, group name for output (optional)

- asso_key:

  Character, associated metadata key (optional)

- description:

  Character, description of the metadata (optional)

## Value

Modified STID object with added metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# Add custom metadata
STID_obj <- AddMetaData(ist_obj,
                       meta_key = "custom_annotation",
                       add_data = custom_metadata,
                       description = "Manual cell annotations")
} # }
```
