# The STID_analysis Class

A container for all analysis results and metadata within an STID object.
Stores project information, metadata records, and both single-sample and
multi-sample niche analysis results.

## Slots

- `STID_info`:

  List containing comprehensive project information:

  - samp_info: Sample-level metadata

  - data_info: Data structure information

  - project_info: Project-level metadata

  - comment_info: Additional comments and notes

- `meta_data_record`:

  List containing metadata tracking information:

  - meta_data_info: Data frame tracking all metadata additions

  - meta_data_list: List of all stored metadata data frames

- `SingleSampNiche`:

  List of SingleSampNiche objects, one per sample

- `MultiSampNiche`:

  List of MultiSampNiche objects for cross-sample analyses

## See also

Other STID-classes:
[`STID-class`](https://yulongqin.github.io/STID/reference/STID-class.md)
