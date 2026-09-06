# Convert Seurat object to STID object

Creates an STID object from a Seurat object with infection-specific
metadata.

## Usage

``` r
CreateSTIDObject(
  seurat_obj = NULL,
  host_org = NULL,
  pathogen_grp = NULL,
  pathogen_org = NULL,
  samp_colnm = NULL,
  samp_grp_colnm = NULL,
  celltype_colnm = NULL,
  x_colnm = NULL,
  y_colnm = NULL,
  pathogen_genes = NULL,
  data_format = NULL,
  data_platform = NULL,
  binsize = NULL,
  coord_interval = NULL,
  base_unit = NULL,
  project_id = NULL,
  description = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object to be converted

- host_org:

  Character, host organism - "human", "mouse", or "unknown"

- pathogen_grp:

  Character, pathogen group - "virus", "bacteria", "parasite", or
  "unknown"

- pathogen_org:

  Character, specific pathogen organism name

- samp_colnm:

  Character, column name in metadata containing sample IDs

- samp_grp_colnm:

  Character, column name in metadata containing sample groups

- celltype_colnm:

  Character, column name in metadata containing cell types

- x_colnm:

  Character, column name for x coordinates

- y_colnm:

  Character, column name for y coordinates

- pathogen_genes:

  Character vector, names of pathogen genes

- data_format:

  Character, spatial data format. One of "square_grid", "hex_grid", or
  "single_cell". For "single_cell", the default spatial unit, binsize,
  and coord_interval are 1 unless otherwise specified.

- data_platform:

  Character, data platform - "StereoSeq", "Visium", "VisiumHD",
  "SlideSeq", or "unknown"

- binsize:

  Numeric, bin size for spatial data

- coord_interval:

  Numeric, coordinate coord_interval

- base_unit:

  Numeric, base unit for spatial coordinates, and the unit is mm

- project_id:

  Character, project identifier

- description:

  Character, project description

- ...:

  Additional arguments passed to
  [`as.STID`](https://yulongqin.github.io/STID/reference/as.STID.md)

## Value

An STID object

## Examples

``` r
if (FALSE) { # \dontrun{
ist_obj <- CreateSTIDObject(
  seurat_obj = seurat_object,
  samp_colnm = "sample_id",
  celltype_colnm = "cell_type",
  host_org = "human",
  pathogen_grp = "virus",
  pathogen_org = "SARS-CoV-2",
  data_format = "square_grid",
  data_platform = "StereoSeq"
)
} # }
```
