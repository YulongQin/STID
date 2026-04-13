# Common analysis theme

A versatile ggplot2 theme for general analysis plots with options for
continuous or discrete x-axis.

## Usage

``` r
theme_common(base_size = 11, x_type = "continuous")
```

## Arguments

- base_size:

  Numeric, base font size (default: 11)

- x_type:

  Character, x-axis type - "continuous" or "discrete" (default:
  "continuous")

## Value

A ggplot2 theme object

## Examples

``` r
if (FALSE) { # \dontrun{
# For continuous x-axis
ggplot(continuous_data, aes(x, y)) +
  geom_point() +
  theme_common(x_type = "continuous")

# For discrete x-axis with rotated labels
ggplot(discrete_data, aes(category, value)) +
  geom_bar() +
  theme_common(x_type = "discrete")
} # }
```
