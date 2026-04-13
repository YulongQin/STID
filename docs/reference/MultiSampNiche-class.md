# The MultiSampNiche Class

A container for comparative niche analysis across multiple samples.
Stores integrated information from multiple SingleSampNiche objects for
cross-sample comparisons.

## Slots

- `samp_info`:

  List containing multi-sample information:

  - multi_id: Character, unique identifier for the multi-sample analysis

  - samp_id: Character vector, sample identifiers included

  - samp_grp: Character, sample group assignments

  - compare_mode: Character, comparison mode ("Comparative" or
    "Temporal")

- `niche_info`:

  List of data frames containing niche metadata for each niche key
  across samples

- `niche_cells`:

  List of data frames containing cell-level information for each niche
  key across samples

- `niche_genes`:

  List of data frames containing gene-level information for each niche
  key across samples

## See also

Other niche-classes:
[`SingleSampNiche-class`](https://yulongqin.github.io/STID/reference/SingleSampNiche-class.md)
