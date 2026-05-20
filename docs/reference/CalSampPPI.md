# Predict Niche Protein-Protein Interactions

Predicts pathogen-host protein-protein interactions using BLAST homology
search and STRING database integration. Constructs and visualizes PPI
networks.

## Usage

``` r
CalSampPPI(
  STID_obj = NULL,
  p.fasta_path = NULL,
  h.fasta_path = NULL,
  p.symbol2protein_path = NULL,
  h.symbol2protein_path = NULL,
  p.features = NULL,
  h.features = NULL,
  BLAST_tool = "blastp",
  dbtype = "prot",
  BLAST_args = "-evalue 1e-10 -max_target_seqs 5 -max_hsps 50 -num_threads 4",
  BLAST_fil_args = "pident > 40, evalue < 1e-10, bitscore > 100, length > 100",
  version_STRING = "12.0",
  score_thre_only = 100,
  score_thre_all = 700,
  return_data = TRUE,
  grp_nm = NULL,
  dir_nm = "M3_CalSampPPI"
)
```

## Arguments

- STID_obj:

  An STID object containing niche analysis results

- p.fasta_path:

  Character, path to pathogen protein FASTA file

- h.fasta_path:

  Character, path to host protein FASTA file

- p.symbol2protein_path:

  Character, pathogen gene symbol to protein ID mapping file

- h.symbol2protein_path:

  Character, host gene symbol to protein ID mapping file

- p.features:

  Character vector, pathogen genes to include

- h.features:

  Character vector, host genes to include

- BLAST_tool:

  Character, BLAST tool to use (default: "blastp")

- dbtype:

  Character, database type (default: "prot")

- BLAST_args:

  Character, BLAST search arguments (default: "-evalue 1e-10
  -max_target_seqs 5 -max_hsps 50 -num_threads 4")

- BLAST_fil_args:

  Character, BLAST result filtering criteria

- version_STRING:

  Character, STRING database version (default: "12.0")

- score_thre_only:

  Numeric, score threshold for pathogen-only network (default: 100)

- score_thre_all:

  Numeric, score threshold for full network (default: 700)

- return_data:

  Logical, whether to return results list (default: TRUE)

- grp_nm:

  Character, group name for output organization (default: NULL)

- dir_nm:

  Character, directory name for output (default: "M3_CalSampPPI")

## Value

If return_data = TRUE, returns a list of PPI results

## Examples

``` r
if (FALSE) { # \dontrun{
# Predict pathogen-host PPIs
results <- CalSampPPI(
  STID_obj = STID_obj,
  p.fasta_path = "pathogen_proteins.fasta",
  h.fasta_path = "host_proteins.fasta",
  p.symbol2protein_path = "pathogen_mapping.txt",
  h.symbol2protein_path = "host_mapping.txt"
)
} # }
```
