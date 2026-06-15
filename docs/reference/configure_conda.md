# Configure Python Conda Environment

Configures the Python environment used by STID through reticulate. If
`conda_nm` is not provided, the function tries to use the default
Miniconda environment managed by reticulate; if it is unavailable,
Miniconda will be installed automatically. If `conda_nm` is provided,
the function checks whether the specified Conda environment exists and
then activates it.

## Usage

``` r
configure_conda(conda_nm = NULL, conda_path = NULL)
```

## Arguments

- conda_nm:

  Character, name of an existing Conda environment to use. If `NULL`,
  the default Miniconda path returned by
  [`miniconda_path`](https://rstudio.github.io/reticulate/reference/miniconda_path.html)
  will be used (default: `NULL`).

- conda_path:

  Character, path to a Conda installation or environment. Currently
  reserved for future use and not used directly in this function
  (default: `NULL`).

## Value

Invisibly returns `NULL`. The function is mainly called for its side
effects, including configuring the Python environment and printing
Python configuration information.

## Examples

``` r
if (FALSE) { # \dontrun{
# Use the default reticulate Miniconda environment
configure_conda()

# Use an existing Conda environment
configure_conda(conda_nm = "scanpy")
} # }
```
