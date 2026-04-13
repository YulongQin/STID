#' Gene and Geneset File Structure
#'
#' A nested list structure containing file paths for gene and geneset definitions
#' for human and mouse species. This dataset provides organized access to various
#' gene sets and gene collections used in immune response and metabolic pathway
#' analyses.
#'
#' @format A nested list with 2 top-level elements:
#' \describe{
#'   \item{Human}{List containing gene and geneset collections for human}
#'   \item{Mouse}{List containing gene and geneset collections for mouse}
#' }
#'
#' @details
#' The structure follows this hierarchy:
#' \itemize{
#'   \item Species (Human/Mouse)
#'   \item Type (Gene/Geneset)
#'   \item Category (GO, KEGG, Reactome, MSigDB, etc.)
#'   \item Specific gene set files
#' }
#'
#' @section Human Gene Collections:
#' \describe{
#'   \item{Human_Infection_Immunity_gene}{Genes involved in infection and immunity}
#'   \item{Human_Macrophage_gene}{Macrophage-associated genes}
#'   \item{Human_PRR_gene}{Pattern recognition receptor genes}
#' }
#'
#' @section Human Geneset Collections:
#' \describe{
#'   \item{GO_BP_Detect_bacterial}{GO biological process terms for bacterial detection}
#'   \item{GO_BP_Detect_viral}{GO biological process terms for viral detection}
#'   \item{MSigDB_Hallmark}{MSigDB hallmark gene sets}
#'   \item{PCD_geneset}{Programmed cell death gene sets}
#'   \item{KEGG_Detect_bacterial}{KEGG pathways for bacterial detection}
#'   \item{KEGG_Detect_parasitic}{KEGG pathways for parasitic detection}
#'   \item{KEGG_Detect_viral}{KEGG pathways for viral detection}
#'   \item{KEGG_Metabolism}{KEGG metabolism pathways}
#'   \item{Reactome_Metabolism}{Reactome metabolism pathways}
#' }
#'
#' @section Mouse Gene Collections:
#' \describe{
#'   \item{Mouse_Infection_Immunity_gene}{Mouse genes involved in infection and immunity}
#'   \item{Mouse_Macrophage_gene}{Mouse macrophage-associated genes}
#'   \item{Mouse_PRR_gene}{Mouse pattern recognition receptor genes}
#' }
#'
#' @section Mouse Geneset Collections:
#' \describe{
#'   \item{GO_BP_Detect_bacterial}{Mouse GO biological process terms for bacterial detection}
#'   \item{GO_BP_Detect_viral}{Mouse GO biological process terms for viral detection}
#'   \item{KEGG_Detect_bacterial}{Mouse KEGG pathways for bacterial detection}
#'   \item{KEGG_Detect_parasitic}{Mouse KEGG pathways for parasitic detection}
#'   \item{KEGG_Detect_viral}{Mouse KEGG pathways for viral detection}
#'   \item{KEGG_Metabolism}{Mouse KEGG metabolism pathways}
#'   \item{HotSpot_geneset}{Mouse HotSpot gene sets}
#'   \item{MSigDB_Hallmark}{Mouse MSigDB hallmark gene sets}
#'   \item{PCD_geneset}{Mouse programmed cell death gene sets}
#'   \item{Reactome_Metabolism}{Mouse Reactome metabolism pathways}
#' }
#'
#' @source Internal package data compiled from various public databases:
#' \itemize{
#'   \item GO: Gene Ontology Consortium
#'   \item KEGG: Kyoto Encyclopedia of Genes and Genomes
#'   \item Reactome: Reactome Pathway Database
#'   \item MSigDB: Molecular Signatures Database
#' }
#'
#' @examples
#' # Access human gene sets
#' Gene_Geneset$Human$Gene$Human_Macrophage_gene
#'
#' # List all mouse KEGG pathways
#' names(Gene_Geneset$Mouse$Geneset$KEGG)
#'
#' # Get all human geneset categories
#' names(Gene_Geneset$Human$Geneset)
#'
"Gene_Geneset"



