# The SingleSampNiche Class

A container for niche analysis results from a single sample. Stores
information about detected niches, including their spatial coordinates,
associated cells, and gene expression patterns.

## Slots

- `samp_info`:

  List containing sample information:

  - samp_id: Character, sample identifier

- `niche_info`:

  Data frame containing metadata for each niche:

  - niche_key: Character, unique identifier for the niche

  - meta_key: Character, reference to metadata source

  - ROI_type: Character, type of ROI ("ROI", "Spot", or "Unknown")

  - pos_colnm: Character, column name for positive spot labels

  - center_colnm: Character, column name for center spot labels

  - edge_colnm: Character, column name for edge spot labels

  - all_label_colnm: Character, column name for all spot labels

  - all_dist_colnm: Character, column name for distance to center

  - niche_cells_num: Integer, number of cells in the niche

  - niche_genes_num: Integer, number of genes analyzed

  - description: Character, description of the niche

- `niche_cells`:

  List containing data frames of cell-level information for each niche
  key

- `niche_genes`:

  List containing data frames of gene-level information for each niche
  key

## See also

Other niche-classes:
[`MultiSampNiche-class`](https://yulongqin.github.io/STID/reference/MultiSampNiche-class.md)
