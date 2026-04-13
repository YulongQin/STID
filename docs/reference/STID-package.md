# STID: Spatial Transcriptomics toolkit for Infectious Diseases

A comprehensive toolkit for quality control, analysis, and exploration
of spatial transcriptomics data. 'STID' aims to enable users to identify
spatial niches, interpret cell-cell communication contexts, and
integrate diverse types of spatial omics measurements (e.g., Visium,
Slide-seq, MERFISH).

The package provides a unified workflow extending the Seurat ecosystem,
featuring:

- **Niche Detection**: Robust algorithms to define spatial
  microenvironments based on cellular composition.

- **Spatial Aggregation**: Metrics to quantify cell clustering and
  spatial organization.

- **Contextual Interaction**: Enhanced ligand-receptor analysis
  incorporating spatial distance constraints.

- **Advanced Visualization**: Specialized plotting functions for spatial
  maps, niche heatmaps, and interaction networks.

See Yulong Qin et al. (Year) doi:10.xxxx/xxxxx for methodological
details on niche detection,

## Package options

STID uses the following [`options`](https://rdrr.io/r/base/options.html)
to configure behaviour:

- `STID.memsafe`:

  Global option to call [`gc()`](https://rdrr.io/r/base/gc.html) after
  memory-intensive operations (e.g., large matrix permutations). This
  helps clean up the R session memory and prevents swap space usage.
  Defaults to `TRUE`. Setting to `FALSE` can speed up computation in
  high-RAM environments.

- `STID.parallel.threads`:

  Controls the default number of threads for parallel operations (using
  foreach and parallel). Defaults to half of the available physical
  cores. Set to `1` to disable parallelization.

- `STID.warn.seurat.v5`:

  Show warning if the input Seurat object is v5 format and specific
  v4-compatible layers are missing.

- `STID.check.dots`:

  For functions that have `...` as a parameter, this controls the
  behavior when an argument isn't used. Can be one of `"warn"`,
  `"stop"`, or `"silent"`. Defaults to `"warn"`.

- `STID.msg.nichenet`:

  Show message about using cached NicheNet networks to speed up
  initialization.

- `STID.msg.cellchat`:

  Show message about the version of CellChatDB being loaded.

## References

If you use STID in your publication, please cite:

- YulongQin, et al. "STID: Spatial Transcriptomics toolkit for
  Infectious Diseases." *Journal Name*, Year. DOI: 10.xxxx/xxxxx

- Stuart T, Butler A, et al. (2019) "Comprehensive Integration of
  Single-Cell Data." *Cell*. doi:10.1016/j.cell.2019.05.031 (For Seurat
  foundation)

## See also

Useful links:

- Homepage: <https://github.com/YulongQin/STID>

- Documentation: <https://YulongQin.github.io/STID/>

- Report bugs at <https://github.com/YulongQin/STID/issues>

## Author

**Maintainer**: Yulong Qin <qyl3700@foxmail.com> (ORCID:
0009-0009-2761-0750)

**Contributors**:

**Funders**:
