# Get STID global constants

Returns a list of global constants used throughout the package.

## Usage

``` r
get_STID_globals()
```

## Value

A list containing supported formats, hosts, pathogens, etc.

## Examples

``` r
if (FALSE) { # \dontrun{
globals <- get_STID_globals()
print(globals$supported_formats)
} # }
```
