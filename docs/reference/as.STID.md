# Convert Seurat object to STID object

Converts a Seurat object to an STID (Integrated Single-cell
Transcriptomics) object, adding infection-specific metadata and analysis
structures.

## Usage

``` r
as.STID(seurat_obj, ...)

# S4 method for class 'Seurat'
as.STID(
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
  binsize = NULL,
  interval = NULL,
  project_id = NULL,
  description = NULL,
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object to be converted

- ...:

  Additional arguments passed to methods

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

  Character, column name for x coordinates (optional)

- y_colnm:

  Character, column name for y coordinates (optional)

- pathogen_genes:

  Character vector, names of pathogen genes (optional)

- data_format:

  Character, data format - "StereoSeq" or "Visium"

- binsize:

  Numeric, bin size for spatial data (optional)

- interval:

  Numeric, coordinate interval (optional)

- project_id:

  Character, project identifier (optional)

- description:

  Character, project description (optional)

## Value

An STID object containing the converted data with infection-specific
metadata and analysis structures

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert Seurat object to STID object
ist_obj <- as.STID(
  seurat_obj = seurat_object,
  samp_colnm = "sample_id",
  celltype_colnm = "cell_type",
  host_org = "human",
  pathogen_grp = "virus",
  pathogen_org = "SARS-CoV-2",
  data_format = "StereoSeq"
)
} # }
```
