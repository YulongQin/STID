# Calculate Niche Organization and Spatial Entropy (OSE)

Performs comprehensive spatial niche analysis including region
segmentation, entropy calculation, and visualization of spatial
organization patterns. This function processes spatial transcriptomics
data to identify organizational niches and calculate adjusted spatial
entropy metrics.

## Usage

``` r
CalNicheOSE(x)
```

## Arguments

- x:

  Placeholder parameter for future expansion

## Value

NULL (invisible), generates multiple output files and visualizations

## Details

The function performs the following analyses:

1.  **Data Preprocessing**: Loads and normalizes spatial transcriptomics
    data

2.  **Spatial Region Segmentation**: Uses iterative k-means-like
    algorithm to partition tissue into spatial regions

3.  **Spatial Entropy Calculation**: Computes adjusted entropy based on
    cell type diversity in local neighborhoods

4.  **Cell Type Proportion Analysis**: Calculates cell type composition
    for each region

5.  **Boundary Detection**: Identifies boundaries between different
    regions

6.  **Visualization**: Generates multiple plots including region maps,
    cell type distributions, pie charts, and entropy heatmaps

## Parameters

The function includes extensive parameter settings based on analysis
mode:

- `mode`: Analysis resolution - "bin50", "bin200", or "pub"

- `binsize`: Spatial bin size (50, 200, or 1 based on mode)

- `window`: Window size for initial region segmentation

- `m`: Weight parameter balancing spatial distance and expression
  distance

- `nPC`: Number of principal components to use

- `min_spot_num`: Minimum spots per region for inclusion

## Outputs

Generates files in outputdata/ and photo/ directories:

- `*_center_df.txt`: Coordinates of region centers

- `*_data.txt`: Cell-level assignments to regions

- `*_entropy_results.txt`: Raw entropy calculations

- `*_entropy_results_fil.txt`: Filtered entropy results

- `*_celltype_percent.txt`: Cell type proportions per region

- `*_celltype_percent_fil.txt`: Filtered cell type proportions

- `*_segment.txt`: Region boundary coordinates

- Multiple visualization PDFs and PNGs

## Examples

``` r
if (FALSE) { # \dontrun{
# Run niche analysis with default settings
CalNicheOSE()

# Modify parameters within the function before running
# Set mode <- "bin200" for lower resolution analysis
# Set only_plot <- TRUE to regenerate plots from existing data
} # }
```
