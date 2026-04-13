# Predefined color palettes

A collection of color palettes for consistent visualization across the
package. Includes palettes for discrete and continuous scales, with
options for different backgrounds.

## Usage

``` r
COLOR_DIS_CON

COLOR_LIST
```

## Format

- COLOR_DIS_CON:

  List with discrete (2 colors) and continuous (5 colors) palettes

- COLOR_LIST:

  Named list of various color palettes:

  - PALETTE_WHITE_BGExtended palette for white background

  - PALETTE_BLACK_BGExtended palette for black background

  - PALETTE_7_VECTOR7-color palette with variations

  - PALETTE_7_CLASSICClassic 7-color palette

  - PALETTE_7_LIGHLight 7-color palette

  - PALETTE_9_CLASSICClassic 9-color palette

  - PALETTE_9_VIRIDISViridis-based 9-color palette

An object of class `list` of length 7.

## Examples

``` r
if (FALSE) { # \dontrun{
# Use discrete palette
ggplot(data, aes(x, y, color = group)) +
  geom_point() +
  scale_color_manual(values = COLOR_DIS_CON$dis)

# Use continuous palette
ggplot(data, aes(x, y, color = value)) +
  geom_point() +
  scale_color_gradientn(colors = COLOR_DIS_CON$con)

# Use white background palette
ggplot(data, aes(x, y, color = category)) +
  geom_point() +
  scale_color_manual(values = COLOR_LIST$PALETTE_WHITE_BG)
} # }
```
