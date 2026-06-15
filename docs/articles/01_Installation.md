# Installation

## Overview

This vignette describes a reproducible installation workflow for `STID`.
Because `STID` integrates spatial transcriptomics, niche detection,
cell-cell communication, regulatory-network inference, and host-pathogen
analysis, the package depends on a relatively broad R, Bioconductor,
GitHub, and system-software ecosystem.

This vignette covers:

- preparing a clean R library for `STID`;
- configuring CRAN and Bioconductor mirrors;
- installing core R, Bioconductor, and GitHub dependencies;
- installing `STID` from GitHub;
- checking whether required and suggested packages are available;

> **Note:** Before running the commands below, update the CRAN mirror,
> Bioconductor mirror, library path, R version, compiler settings, and
> optional dependency list according to the operating system and network
> environment used for the analysis.

> **Reference code:** For more runnable code examples, please refer to
> the [article
> code](https://github.com/YulongQin/STID/tree/main/article_code).

## Installation strategy

A controlled installation environment is recommended. The general
workflow is:

1.  confirm that R and system build tools are available;
2.  configure a dedicated library path for `STID`;
3.  configure CRAN and Bioconductor mirrors;
4.  install package-management utilities;
5.  install heavy or failure-prone dependencies first;
6.  install `STID` from GitHub;
7.  verify core and optional dependencies.

> **Note:** If an older package version is already loaded from the
> default library, simply calling `.libPaths("./library_STID/")` may not
> fully solve version conflicts. In that case, restart R and update or
> remove the conflicting package from the default library before
> reinstalling `STID`.

## Configure the R installation environment

### Configure mirrors

Use a CRAN mirror and Bioconductor mirror that are stable in the current
network environment.

``` r
# Recommended for users in China
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# Recommended for users outside China
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(BioC_mirror = "https://bioconductor.org")
```

### Configure an isolated library path

A dedicated library path helps avoid conflicts with an existing R
environment and makes the installation easier to reproduce.

``` r
dir.create("./library_STID/", showWarnings = FALSE, recursive = TRUE)
.libPaths(normalizePath("./library_STID/"))
.libPaths()
```

> **Note:** If older dependency versions already exist in the default
> library (for example, `Matrix`, `rlang`, or `SeuratObject`), version
> conflicts may still occur during installation. It is recommended to
> restart R and prioritize upgrading the conflicting packages in the
> default library, or reinstall everything in a fresh R environment.

### Configure download and compilation options

``` r
# Increase download timeout for large packages or unstable networks
options(timeout = 300)

# Use parallel compilation when source installation is required
n_cores <- max(1, parallel::detectCores() - 2)
Sys.setenv(MAKEFLAGS = paste0("-j", n_cores))

# Reduce debug information and use optimized compilation flags
Sys.setenv(
  PKG_CFLAGS = "-O2 -g0 -DNDEBUG",
  PKG_CXXFLAGS = "-O2 -g0 -DNDEBUG"
)
```

## Install package-management utilities

Install `remotes`, `devtools`, and `BiocManager` before installing
`STID` and its dependencies.

``` r
install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs)
  }
  invisible(missing_pkgs)
}

install_if_missing(c("remotes", "devtools", "BiocManager"))
```

## Pre-installation requirements

We recommend using R 4.5 or later. Installing packages becomes
challenging on R versions older than 4.5, as numerous dependencies will
require manual compilation from source. We suggest installing the
packages listed below first, and source-based installation is advised.

``` r
# Core infrastructure package
install.packages("rlang") # >= 1.1.7


# These packages require network access (e.g., GitHub or external sources).
# Installation failures (e.g., "Error in download.file") are often caused by network issues.
# Retry installation or consider using a proxy/VPN if necessary.
# Note: this is not a complete list, but these are the most failure-prone packages.
install.packages("openssl") # >= 2.3.5
install.packages("curl") # >= 7.0.0 
install.packages("systemfonts")


# These packages are recommended to be installed as binary versions
# to avoid compilation errors and reduce installation time
install.packages("Matrix", type = "binary") # >= 1.6.5
install.packages("RcppArmadillo", type = "binary")
install.packages("uwot", type = "binary")
install.packages("units", type = "binary")
install.packages("ggforce", type = "binary")
```

> **Note:** These issues mainly occur in older versions of R (e.g., R
> 4.2). R versions 4.5 and above ship with binary packages compatible
> with the latest source packages, so you can skip these installation
> steps if you are using R ≥ 4.5.

## Install STID

Install `STID` directly from GitHub. When prompted to update
dependencies, selecting `3: None` is often safer in an already
configured environment because it avoids unexpectedly overwriting
packages such as `Matrix`, `SeuratObject`, or `Seurat`.

``` r
remotes::install_github(
  "YulongQin/STID",
  dependencies = TRUE
)
```

> **Note:** The full installation may take anywhere from 10 minutes to
> over 1 hour, depending on the operating system, R version, network
> speed, and whether packages are installed from source. With R 4.5 or
> later, installation typically takes around 10 minutes, whereas older R
> versions may require extensive source compilation and can take
> approximately 1 hour or longer.

``` r
library(STID)
packageVersion("STID")
```

## Check dependencies

You can check the dependencies of STID
[here](https://github.com/YulongQin/STID/blob/main/DESCRIPTION).

The `Imports` packages are required for core `STID` functions. Packages
in `Suggests` are not always installed automatically and should be
installed manually before running the corresponding downstream analysis.

To ensure all required packages are installed with correct versions, the
following utility function can be used.

``` r
pkg_imports <- c(
  "anndata", "AnnotationDbi", "assertthat", "AUCell", "Biobase", "Biostrings",
  "clusterProfiler", "concaveman", "data.table", "dbscan", "distances",
  "doParallel", "DOSE", "dplyr", "DT", "FNN", "foreach", "furrr",
  "future", "ggdensity", "ggforce", "ggh4x", "ggnewscale", "ggplot2",
  "ggprism", "ggraph", "ggrepel", "ggsignif", "ggvenn", "GO.db", "GOSemSim",
  "igraph", "imager", "magrittr", "Matrix", "methods", "Mfuzz", "paletteer",
  "patchwork", "proxy", "psych", "purrr", "rBLAST", "reticulate",
  "rjson", "rlang", "scales", "Seurat", "SeuratDisk", "SeuratObject",
  "shiny", "shinycssloaders", "sp", "STRINGdb", "stringr",
  "SummarizedExperiment", "tibble", "tidyr", "tidyselect", "UCell", "viridis", "zip"
)

pkg_suggests <- c(
  "CellChat", "knitr", "mistyR", "nichenetr", "org.Hs.eg.db",
  "org.Mm.eg.db", "rmarkdown", "SingleR", "testthat"
)

pkg_min_versions <- c(
  ggplot2 = "4.0.2",
  rlang = "1.1.7",
  Matrix = "1.6.5",
  SeuratObject = "5.0.2",
  Seurat = "5.1.0"
)
```

``` r
check_stid_packages <- function(pkg_list, min_versions = NULL) {
  installed_pkgs <- rownames(installed.packages())
  installed_ok <- intersect(pkg_list, installed_pkgs)
  missing_pkgs <- setdiff(pkg_list, installed_pkgs)

  version_wrong <- character()
  if (!is.null(min_versions)) {
    for (pkg in names(min_versions)) {
      if (pkg %in% installed_ok) {
        current_version <- as.character(packageVersion(pkg))
        required_version <- min_versions[[pkg]]
        if (utils::compareVersion(current_version, required_version) < 0) {
          version_wrong <- c(
            version_wrong,
            sprintf(
              "%s (requires >= %s, current = %s)",
              pkg, required_version, current_version
            )
          )
        }
      }
    }
  }

  cat("==================== Package Installation Check ====================\n")
  cat("Installed packages: ", length(installed_ok), "\n", sep = "")
  cat("Missing packages: ", length(missing_pkgs), "\n", sep = "")

  if (length(missing_pkgs) > 0) {
    cat("\nMissing list:\n")
    cat(paste0("  - ", missing_pkgs, "\n"), sep = "")
  } else {
    cat("\nAll packages in this list are installed.\n")
  }

  if (length(version_wrong) > 0) {
    cat("\nVersion mismatches:\n")
    cat(paste0("  - ", version_wrong, "\n"), sep = "")
  }
  cat("====================================================================\n")

  invisible(list(
    installed = installed_ok,
    missing = missing_pkgs,
    version_mismatch = version_wrong
  ))
}

imports_check <- check_stid_packages(pkg_imports, pkg_min_versions)
suggests_check <- check_stid_packages(pkg_suggests)
```

## Install missing dependencies from the check result

In some environments (e.g., restricted networks or incompatible
binaries), automatic installation may fail. The following steps provide
manual installation options.

### Seurat dependencies

``` r
remotes::install_version('rlang', version = '1.1.7')
remotes::install_version('Matrix', version = '1.6.5')
remotes::install_version('SeuratObject', version = '5.0.2')
remotes::install_version('Seurat', version = '5.1.0')
```

### Bioconductor dependencies

``` r
BiocManager::install(c(
  "AnnotationDbi",
  "AUCell",
  "Biostrings",
  "clusterProfiler",
  "DOSE",
  "GO.db",
  "GOSemSim",
  "SingleR",
  "STRINGdb",
  "SummarizedExperiment",
  "UCell",
  "org.Hs.eg.db",
  "org.Mm.eg.db"
), ask = FALSE, update = FALSE)
```

### GitHub dependencies

Some downstream analyses use packages that are usually installed from
GitHub. Install these packages only if the corresponding analysis
modules are needed.

``` r
devtools::install_github("jinworks/CellChat")
devtools::install_github("saezlab/mistyR")
devtools::install_github("saeyslab/nichenetr")
devtools::install_github("mhahsler/rBLAST")
devtools::install_github("mojaveazure/seurat-disk#198")
```

> **Note:** When GitHub installation fails with
> `Error in download.file`, the most common causes are network
> interruption, rate limits, or proxy configuration. Retry the command,
> configure a proxy if needed, or install the package from a downloaded
> source archive.

## Notes

- Use a clean library path and restart R before installation to avoid
  stale namespaces.
- Use R 4.5 or higher when possible to reduce source-compilation and
  dependency-resolution problems.
- Install heavy or failure-prone dependencies before installing `STID`
  from GitHub.
- Suggested packages are not always installed automatically; install
  them manually before running optional downstream analyses.
- For reproducible projects, record the R version, package versions,
  system libraries, and installation date.

## Next steps

After installing STID and the required dependencies, you are ready to
load spatial transcriptomic data into R. In the next vignette, we will
preprocess the input data, perform basic quality control and annotation,
and convert the processed object into the `STID` class.

## Session information

``` r
sessionInfo()
#> R version 4.2.0 (2022-04-22 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 22000)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
#> [2] LC_CTYPE=Chinese (Simplified)_China.utf8   
#> [3] LC_MONETARY=Chinese (Simplified)_China.utf8
#> [4] LC_NUMERIC=C                               
#> [5] LC_TIME=Chinese (Simplified)_China.utf8    
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.35     R6_2.6.1          jsonlite_1.8.8    lifecycle_1.0.5  
#>  [5] evaluate_1.0.1    cachem_1.1.0      rlang_1.1.7       cli_3.6.5        
#>  [9] rstudioapi_0.15.0 fs_1.6.3          jquerylib_0.1.4   bslib_0.8.0      
#> [13] ragg_1.3.0        rmarkdown_2.29    pkgdown_2.2.0     textshaping_0.3.6
#> [17] desc_1.4.3        tools_4.2.0       htmlwidgets_1.6.4 yaml_2.3.10      
#> [21] xfun_0.49         fastmap_1.2.0     compiler_4.2.0    systemfonts_1.0.4
#> [25] htmltools_0.5.8.1 knitr_1.49        sass_0.4.9
```
