# Gene and Geneset File Structure

A nested list structure containing file paths for gene and geneset
definitions for human and mouse species. This dataset provides organized
access to various gene sets and gene collections used in immune response
and metabolic pathway analyses.

## Usage

``` r
Gene_Geneset
```

## Format

A nested list with 2 top-level elements:

- Human:

  List containing gene and geneset collections for human

- Mouse:

  List containing gene and geneset collections for mouse

## Source

Internal package data compiled from various public databases:

- GO: Gene Ontology Consortium

- KEGG: Kyoto Encyclopedia of Genes and Genomes

- Reactome: Reactome Pathway Database

- MSigDB: Molecular Signatures Database

## Details

The structure follows this hierarchy:

- Species (Human/Mouse)

- Type (Gene/Geneset)

- Category (GO, KEGG, Reactome, MSigDB, etc.)

- Specific gene set files

## Human Gene Collections

- Human_Infection_Immunity_gene:

  Genes involved in infection and immunity

- Human_Macrophage_gene:

  Macrophage-associated genes

- Human_PRR_gene:

  Pattern recognition receptor genes

## Human Geneset Collections

- GO_BP_Detect_bacterial:

  GO biological process terms for bacterial detection

- GO_BP_Detect_viral:

  GO biological process terms for viral detection

- MSigDB_Hallmark:

  MSigDB hallmark gene sets

- PCD_geneset:

  Programmed cell death gene sets

- KEGG_Detect_bacterial:

  KEGG pathways for bacterial detection

- KEGG_Detect_parasitic:

  KEGG pathways for parasitic detection

- KEGG_Detect_viral:

  KEGG pathways for viral detection

- KEGG_Metabolism:

  KEGG metabolism pathways

- Reactome_Metabolism:

  Reactome metabolism pathways

## Mouse Gene Collections

- Mouse_Infection_Immunity_gene:

  Mouse genes involved in infection and immunity

- Mouse_Macrophage_gene:

  Mouse macrophage-associated genes

- Mouse_PRR_gene:

  Mouse pattern recognition receptor genes

## Mouse Geneset Collections

- GO_BP_Detect_bacterial:

  Mouse GO biological process terms for bacterial detection

- GO_BP_Detect_viral:

  Mouse GO biological process terms for viral detection

- KEGG_Detect_bacterial:

  Mouse KEGG pathways for bacterial detection

- KEGG_Detect_parasitic:

  Mouse KEGG pathways for parasitic detection

- KEGG_Detect_viral:

  Mouse KEGG pathways for viral detection

- KEGG_Metabolism:

  Mouse KEGG metabolism pathways

- HotSpot_geneset:

  Mouse HotSpot gene sets

- MSigDB_Hallmark:

  Mouse MSigDB hallmark gene sets

- PCD_geneset:

  Mouse programmed cell death gene sets

- Reactome_Metabolism:

  Mouse Reactome metabolism pathways

## Examples

``` r
# Access human gene sets
Gene_Geneset$Human$Gene$Human_Macrophage_gene
#> # A tibble: 15 × 3
#>    M1      M2     Fusion 
#>    <chr>   <chr>  <chr>  
#>  1 FCGR3A  CD163  ADAM9  
#>  2 CD86    MRC1   CD44   
#>  3 CD14    ARG1   CD81   
#>  4 NOS2    CHI3L1 DCSTAMP
#>  5 IL6     RETN   OCSTAMP
#>  6 TNF     IL10   STAT1  
#>  7 IL1B    TGFB1  TREM2  
#>  8 IL12B   PPARG  TYROBP 
#>  9 CXCL9   KLF4   NA     
#> 10 CXCL10  IL4R   NA     
#> 11 CD80    CCL17  NA     
#> 12 HLA-DRA CCL22  NA     
#> 13 IRF5    NA     NA     
#> 14 SOCS3   NA     NA     
#> 15 FPR2    NA     NA     

# List all mouse KEGG pathways
names(Gene_Geneset$Mouse$Geneset$KEGG)
#> [1] "Mouse_KEGG_Detect_bacterial_geneset" "Mouse_KEGG_Detect_parasitic_geneset"
#> [3] "Mouse_KEGG_Detect_viral_geneset"     "Mouse_KEGG_Metabolism_geneset"      

# Get all human geneset categories
names(Gene_Geneset$Human$Geneset)
#> [1] "GO"                            "Human_MSigDB_Hallmark_geneset"
#> [3] "Human_PCD_geneset"             "KEGG"                         
#> [5] "Reactome"                     
```
