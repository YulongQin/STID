# Create spatial plots for STID data

This function generates spatial visualizations with discrete or
continuous coloring. Can be used directly with an STID object or with a
custom data frame.

## Usage

``` r
Plot_Spatial(
  STID_obj = NULL,
  meta_key = NULL,
  plot_data = NULL,
  x_colnm = NULL,
  y_colnm = NULL,
  group_by = NULL,
  facet_grpnm = NULL,
  datatype = "discrete",
  col = COLOR_DIS_CON,
  pt_size = 1.1,
  vmin = NULL,
  vmax = "p99",
  title = NULL,
  subtitle = NULL,
  black_bg = FALSE
)
```

## Arguments

- STID_obj:

  An STID object (optional if plot_data provided)

- meta_key:

  Character, metadata key to use (default: NULL)

- plot_data:

  Data frame containing plotting data (optional if STID_obj provided)

- x_colnm:

  Character, column name for x coordinates

- y_colnm:

  Character, column name for y coordinates

- group_by:

  Character, column name for grouping/coloring variable

- facet_grpnm:

  Character, column name for faceting (default: NULL)

- datatype:

  Character, data type - "discrete" or "continuous" (default:
  "discrete")

- col:

  Color palette for visualization (default: COLOR_DIS_CON)

- pt_size:

  Numeric, point size (default: 1.1)

- vmin:

  Numeric or character, minimum value for color scale (default: NULL)

- vmax:

  Numeric or character, maximum value for color scale (default: "p99")

- title:

  Character, plot title (default: NULL)

- subtitle:

  Character, plot subtitle (default: NULL)

- black_bg:

  Logical, whether to use black background (default: FALSE)

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
# Create spatial plot from STID object
p <- Plot_Spatial(
  STID_obj = STID_object,
  group_by = "Pathogen_Score",
  datatype = "continuous",
  black_bg = TRUE
)

# Create spatial plot from custom data frame
p <- Plot_Spatial(
  plot_data = my_data,
  x_colnm = "x_coord",
  y_colnm = "y_coord",
  group_by = "cell_type",
  datatype = "discrete"
)
} # }
```
