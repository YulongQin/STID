# STID: Spatial Transcriptomics toolkit for Infectious Diseases

Provides a standardized and extensible framework for
infection-associated spatial transcriptomic analysis. The package
integrates with the Seurat ecosystem and incorporates Python-based
modules to provide a full workflow for infection-related spatial
transcriptomic data.

## Details

STID introduces an infection-specific data structure and supports:
pathogen background correction, infection-associated spot and niche
identification, single-sample and multi-sample analysis, as well as
temporal analysis across infection stages.

The framework enables systematic characterization of: structural
features of infection niches, cellular composition, molecular functions,
host–pathogen interactions, and pathogen-infected or host-responsive
niches.

## Main functionalities

- Data conversion: `h5ad2rds`, `rds2h5ad`

- Preprocessing: `Seurat_pipeline`, `anno_SingleR`

- Background correction: `CorrectBackground`

- Spot detection: `SpotDetect_Gene`, `SpotDetect_Geneset`

- Niche identification: `NicheDetect_Lasso`, `NicheDetect_STS`

- Single-sample niche analysis

- Multi-sample niche analysis

- Customized visualization

## Package options

STID initializes several package-level options during loading. Current
options can be viewed with
[`get_STID_options()`](https://yulongqin.github.io/STID/reference/get_STID_options.md).

- `parallel_workers`:

  Default number of workers for parallel computation:
  `min(4, parallel::detectCores()-2)` to avoid excessive CPU usage.

- `default_plan`:

  Default parallel backend used by future:

  - Windows: `"multisession"`

  - Linux/macOS: `"multicore"`

- `future.globals.maxSize`:

  Maximum size of global objects exported in future-based computation:
  default is 1 GB (`1 * 1024^3`), can be increased for large datasets.

## Global constants

Access internal constants via
[`get_STID_globals()`](https://yulongqin.github.io/STID/reference/get_STID_globals.md),
including: supported platforms, hosts, pathogen types, scoring methods,
visualization modes, and sample analysis modes.

## References

Qin Y et al. STID: Spatial Transcriptomics toolkit for Infectious
Diseases. Stuart T, Butler A, et al. Comprehensive Integration of
Single-Cell Data. *Cell* 2019.
[doi:10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031)

## See also

Useful links:

- Homepage: <https://github.com/YulongQin/STID>

- Documentation: <https://YulongQin.github.io/STID>

- Tutorial: <https://YulongQin.github.io/STID/articles>

- Report bugs: <https://github.com/YulongQin/STID/issues>

## Author

**Maintainer:** Yulong Qin Email: <qyl3700@foxmail.com> ORCID:
<https://orcid.org/0009-0009-2761-0750>
