# Spatial plot theme

A ggplot2 theme optimized for spatial transcriptomics plots with minimal
grid lines and clean background.

## Usage

``` r
theme_spatial(base_size = 11)
```

## Arguments

- base_size:

  Numeric, base font size (default: 11)

## Value

A ggplot2 theme object

## Examples

``` r
if (FALSE) { # \dontrun{
ggplot(spatial_data, aes(x, y, color = expression)) +
  geom_point() +
  theme_spatial()
} # }
```
