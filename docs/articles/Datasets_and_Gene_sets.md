# Datasets and gene sets

## Overview

This vignette briefly summarizes the public spatial transcriptomic
datasets of infectious diseases and curated gene resources used by
`STID`.

``` r
library(tidyverse)
library(Seurat)
library(STID)
```

## Public datasets

Table S1 includes eight public mouse spatial transcriptomic datasets
covering viral, bacterial, and parasitic infection models.

| Pathogen type | Pathogen | Abbrev. | Disease | Tissue | Reference | Technology | Repository/accession | n |
|----|----|----|----|----|----|----|----|----|
| Virus | Japanese encephalitis virus | JEV | Japanese encephalitis | Brain | Ou et al. (2026) | Stereo-seq | STT0000076 | 12 |
| Virus | Venezuelan equine encephalitis virus | VEEV | Venezuelan equine encephalitis | Brain | Rangel et al. (2024) | 10x Visium | GSE275201 | 2 |
| Virus | Reovirus type 1 Lang | Reo-T1L | Viral myocarditis | Heart | Mantri et al. (2022) | 10x Visium | GSE189636 | 4 |
| Bacterium | Mycobacterium tuberculosis | MTB | Tuberculosis | Lung | Zhao et al. (2025) | Stereo-seq V2 | CRA018250 | 8 |
| Bacterium | Klebsiella pneumoniae | Kp | Klebsiella pneumoniae pneumonia | Lung | Xu et al. (2022) | 10x Visium | GSE190225 | 4 |
| Parasite | Echinococcus multilocularis | Em | Alveolar echinococcosis | Liver | Ou et al. (2025) | Stereo-seq | STT0000072 | 12 |
| Parasite | Plasmodium berghei ANKA | N/A | Malaria | Liver | Hildebrandt et al. (2024) | 10x Visium | GSE268068 | 4 |
| Parasite | Trypanosoma brucei brucei | Tbb | African trypanosomiasis | Brain | Quintana et al. (2022) | 10x Visium | GSE200642 | 3 |

Raw datasets are available from the repositories reported in the
original publications and Table S1. Processed datasets and curated gene
and gene-set resources are available from
[Figshare](https://doi.org/10.6084/m9.figshare.31839988).

## Curated genes and gene sets

Table S2 summarizes 12 curated resources. Counts are shown as
`Human/Mouse`.

| Type | Name | Main use | Count |
|----|----|----|----|
| Gene | Infection_Immunity | Infection and immune response genes | 391/345 |
| Gene | PRR | Pattern recognition receptor genes | 137/128 |
| Gene | Macrophage | Macrophage polarization/function genes | 35/35 |
| Gene set | GO_BP_Detect_viral | GO biological process terms related to viral infection | 9/12 |
| Gene set | GO_BP_Detect_bacterial | GO biological process terms related to bacterial infection | 10/7 |
| Gene set | KEGG_Detect_viral | KEGG viral infection pathways | 21/21 |
| Gene set | KEGG_Detect_bacterial | KEGG bacterial infection pathways | 11/7 |
| Gene set | KEGG_Detect_parasitic | KEGG parasitic infection pathways | 6/6 |
| Gene set | KEGG_Metabolism | KEGG metabolism pathways | 85/85 |
| Gene set | Reactome_Metabolism | Reactome metabolism pathways | 301/253 |
| Gene set | MSigDB_Hallmark | MSigDB hallmark gene sets | 50/50 |
| Gene set | PCD | Programmed cell death gene sets | 18/18 |

``` r
STID::Gene_Geneset
#> $Human
#> $Human$Gene
#> $Human$Gene$Human_Infection_Immunity_gene
#> # A tibble: 110 × 17
#>    IFN       IRF   IFI   CXCL  CXCR  CCL   CCR   TNF   CSF   IL    NFKB  Attacin
#>    <chr>     <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr>  
#>  1 IFNA1     IRF1… IFIT1 CXCL1 CXCR1 CCL1  CCR1  TNF   CSF1  IL1R2 NFKB… CAMP   
#>  2 IFNAR1    IRF1  IFIT… CXCL… CXCR… CCL2  CCR2  TNFA… CSF1R IL1R… NFKB… DEFA1  
#>  3 IFNA2     IRF2… IFIT… CXCL2 CXCR2 CCL3… CCRL2 TNFR… CSF2  IL1R1 NFKB… DEFA1B 
#>  4 IFNAR2    IRF2  IFIT… CXCL3 CXCR3 CCL3  CCR3  TNFR… CSF2… IL1F… NFKB… DEFA3  
#>  5 IFNAR2-I… IRF2… IFIH1 CXCL5 CXCR4 CCL3… CCR4  TNFA… CSF2… IL1A  NFKB… DEFA4  
#>  6 IFNA4     IRF2… IFIT2 CXCL6 CXCR5 CCL4  CCR5… TNFA… CSF2… IL1R… NFKB1 DEFA5  
#>  7 IFNA5     IRF2… IFIT… CXCL8 CXCR6 CCL4… CCR5  TNFS… CSF2… ILRU… NFKB… DEFA6  
#>  8 IFNA6     IRF3  IFIT… CXCL9 NA    CCL5  CCR6  TNFR… CSF3R ILDR1 NFKB2 DEFA7P 
#>  9 IFNA7     IRF4  IFIT… CXCL… NA    CCL7  CCR7  TNFA… CSF3  IL1RN REL   DEFA8P 
#> 10 IFNA8     IRF5  IFIT… CXCL… NA    CCL8  CCR8  TNFR… NA    IL1R… RELA  DEFA9P 
#> # ℹ 100 more rows
#> # ℹ 5 more variables: ATE <chr>, ROS_NO <chr>, Clq <chr>, Acute_Phase <chr>,
#> #   Antigen <chr>
#> 
#> $Human$Gene$Human_Macrophage_gene
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
#> 
#> $Human$Gene$Human_PRR_gene
#> # A tibble: 38 × 8
#>    TLRs     TLRs_downstream CASP    NLRs   CLRs    SRs    CRs    FcRs   
#>    <chr>    <chr>           <chr>   <chr>  <chr>   <chr>  <chr>  <chr>  
#>  1 TLR1     MYD88           CASP1   NAIP   MRC1    MSR1   CR1L   FCGRT  
#>  2 TLR2     TICAM1          CASP1P2 NAIPP1 CLECL1P SCARB1 CR2    FCGR1BP
#>  3 TLR3     TICAM2          CASP1P1 NAIPP2 CLECL1P CD36   ITGAM  FCGR1BP
#>  4 TLR4     TICAM2-AS1      CASP2   NAIPP3 CLEC1B  SCARB1 ITGAX  FCGR1CP
#>  5 CD14     TIRAP           CASP3P1 NAIPP4 CLEC1A  MARCO  FCGR1A FCGR1A 
#>  6 TLR5     TRAM1L1         CASP3   NODAL  CLEC2L  CD163  NA     FCGR2A 
#>  7 TLR6     TRAM1           CASP4LP NOD1   CLEC2D  NA     NA     FCGR2C 
#>  8 TLR7     TRAM2           CASP4LP NOD2   CLEC2B  NA     NA     FCGR2B 
#>  9 TLR8-AS1 TRAM2-AS1       CASP4   NLRC3  CLEC2A  NA     NA     FCGR3A 
#> 10 TLR8     TRAF1           CASP5   NLRC4  CLEC3B  NA     NA     FCGR3B 
#> # ℹ 28 more rows
#> 
#> 
#> $Human$Geneset
#> $Human$Geneset$GO
#> $Human$Geneset$GO$Human_GO_BP_Detect_bacterial_geneset
#> # A tibble: 376 × 10
#>    GOBP_ANTIBACTERIAL_HUMORAL_RE…¹ GOBP_ANTIBACTERIAL_I…² GOBP_ANTIBACTERIAL_P…³
#>    <chr>                           <chr>                  <chr>                 
#>  1 ADM                             CXCL6                  ELANE                 
#>  2 ANG                             DEFB136                EVPL                  
#>  3 APP                             IL33                   KLK3                  
#>  4 B2M                             MARCHF2                KLK5                  
#>  5 BPI                             MIRLET7B               KLK7                  
#>  6 BPIFA1                          MPEG1                  LGALS4                
#>  7 CALCA                           NFKB1                  MMP7                  
#>  8 CAMP                            NOD2                   PGC                   
#>  9 CTSG                            RASGRP4                SPINK5                
#> 10 DEFA1                           SLC15A2                NA                    
#> # ℹ 366 more rows
#> # ℹ abbreviated names: ¹​GOBP_ANTIBACTERIAL_HUMORAL_RESPONSE,
#> #   ²​GOBP_ANTIBACTERIAL_INNATE_IMMUNE_RESPONSE,
#> #   ³​GOBP_ANTIBACTERIAL_PEPTIDE_PRODUCTION
#> # ℹ 7 more variables:
#> #   GOBP_CELLULAR_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN <chr>,
#> #   GOBP_DETECTION_OF_MOLECULE_OF_BACTERIAL_ORIGIN <chr>, …
#> 
#> $Human$Geneset$GO$Human_GO_BP_Detect_viral_geneset
#> # A tibble: 437 × 9
#>    GOBP_CELLULAR_RESPONSE_TO_VIRUS GOBP_DEFENSE_RESPONS…¹ GOBP_NEGATIVE_REGULA…²
#>    <chr>                           <chr>                  <chr>                 
#>  1 ADAR                            ABCC9                  ATG12                 
#>  2 ARF1                            ABCF3                  ATG5                  
#>  3 ATF2                            ACOD1                  C1QBP                 
#>  4 BAX                             ADAR                   FGL2                  
#>  5 CALR                            ADARB1                 FOXP3                 
#>  6 CCL19                           AGBL4                  ILRUN                 
#>  7 CCL5                            AGBL5                  IRGM                  
#>  8 CHUK                            AICDA                  ITCH                  
#>  9 CXCL10                          AIM2                   MICB                  
#> 10 DDX3X                           AIMP1                  MIR26B                
#> # ℹ 427 more rows
#> # ℹ abbreviated names: ¹​GOBP_DEFENSE_RESPONSE_TO_VIRUS,
#> #   ²​GOBP_NEGATIVE_REGULATION_OF_DEFENSE_RESPONSE_TO_VIRUS
#> # ℹ 6 more variables:
#> #   GOBP_POSITIVE_REGULATION_OF_DEFENSE_RESPONSE_TO_VIRUS_BY_HOST <chr>,
#> #   GOBP_REGULATION_BY_VIRUS_OF_VIRAL_PROTEIN_LEVELS_IN_HOST_CELL <chr>,
#> #   GOBP_REGULATION_OF_DEFENSE_RESPONSE_TO_VIRUS <chr>, …
#> 
#> 
#> $Human$Geneset$Human_MSigDB_Hallmark_geneset
#> # A tibble: 200 × 50
#>    HALLMARK_ADIPOGENESIS HALLMARK_ALLOGRAFT_REJECTION HALLMARK_ANDROGEN_RESPONSE
#>    <chr>                 <chr>                        <chr>                     
#>  1 ABCA1                 AARS1                        ABCC4                     
#>  2 ABCB8                 ABCE1                        ABHD2                     
#>  3 ACAA2                 ABI1                         ACSL3                     
#>  4 ACADL                 ACHE                         ACTN1                     
#>  5 ACADM                 ACVR2A                       ADAMTS1                   
#>  6 ACADS                 AKT1                         ADRM1                     
#>  7 ACLY                  APBB1                        AKAP12                    
#>  8 ACO2                  B2M                          AKT1                      
#>  9 ACOX1                 BCAT1                        ALDH1A3                   
#> 10 ADCY6                 BCL10                        ANKH                      
#> # ℹ 190 more rows
#> # ℹ 47 more variables: HALLMARK_ANGIOGENESIS <chr>,
#> #   HALLMARK_APICAL_JUNCTION <chr>, HALLMARK_APICAL_SURFACE <chr>,
#> #   HALLMARK_APOPTOSIS <chr>, HALLMARK_BILE_ACID_METABOLISM <chr>,
#> #   HALLMARK_CHOLESTEROL_HOMEOSTASIS <chr>, HALLMARK_COAGULATION <chr>,
#> #   HALLMARK_COMPLEMENT <chr>, HALLMARK_DNA_REPAIR <chr>,
#> #   HALLMARK_E2F_TARGETS <chr>, …
#> 
#> $Human$Geneset$Human_PCD_geneset
#> # A tibble: 580 × 18
#>    Apoptosis Pyroptosis Ferroptosis Autophagy Necroptosis Cuproptosis
#>    <chr>     <chr>      <chr>       <chr>     <chr>       <chr>      
#>  1 AATF      BAK1       ABCC1       ABL1      GLUD1       NFE2L2     
#>  2 ABL1      BAX        ACACA       ABL2      GLUD2       NLRP3      
#>  3 ACAA2     CASP1      ACO1        ACER2     ALOX15      ATP7B      
#>  4 ACKR3     CASP3      ACSF2       ADRA1A    FTH1        ATP7A      
#>  5 ACVR1     CASP4      ACSL1       ADRB2     PYG         SLC31A1    
#>  6 ACVR1B    CASP5      ACSL3       AKT1      CAPN1       FDX1       
#>  7 ADORA1    CASP6      ACSL4       AMBRA1    CASP1       LIAS       
#>  8 AEN       CASP8      ACSL5       ATF6      GLNA        LIPT1      
#>  9 AGT       CASP9      ACSL6       ATG101    BAX         LIPT2      
#> 10 AGTR2     CHMP2A     AIFM2       ATG13     BCL2        DLD        
#> # ℹ 570 more rows
#> # ℹ 12 more variables: Parthanatos <chr>, Entoticcelldeath <chr>,
#> #   Netoticcelldeath <chr>, `Lysosome-dependentcelldeath` <chr>,
#> #   Alkaliptosis <chr>, Oxeiptosis <chr>, NETosis <chr>,
#> #   Immunogenic_cell_death <chr>, Anoikis <chr>, Paraptosis <chr>,
#> #   Methuosis <chr>, Entosis <chr>
#> 
#> $Human$Geneset$KEGG
#> $Human$Geneset$KEGG$Human_KEGG_Detect_bacterial_geneset
#> # A tibble: 253 × 11
#>    Bacterial invasion of epithel…¹ Vibrio cholerae infe…² Epithelial cell sign…³
#>    <chr>                           <chr>                  <chr>                 
#>  1 ARPC5                           ATP6V1FP2              ATP6V1FP2             
#>  2 ARPC4                           TCIRG1                 ADAM10                
#>  3 ARPC3                           CFTR                   TCIRG1                
#>  4 ARPC1B                          ADCY3                  NOD1                  
#>  5 ACTR3                           KDELR1                 CHUK                  
#>  6 ACTR2                           SEC61B                 ATP6V1G3              
#>  7 ARPC2                           KDELR2                 MAPK14                
#>  8 WASF2                           KDELR3                 CSK                   
#>  9 MAD2L2                          ADCY9                  IGSF5                 
#> 10 ARPC1A                          ATP6V1G3               ATP6V0E2              
#> # ℹ 243 more rows
#> # ℹ abbreviated names: ¹​`Bacterial invasion of epithelial cells`,
#> #   ²​`Vibrio cholerae infection`,
#> #   ³​`Epithelial cell signaling in Helicobacter pylori infection`
#> # ℹ 8 more variables: `Pathogenic Escherichia coli infection` <chr>,
#> #   Shigellosis <chr>, `Salmonella infection` <chr>, Pertussis <chr>,
#> #   Legionellosis <chr>, `Yersinia infection` <chr>, …
#> 
#> $Human$Geneset$KEGG$Human_KEGG_Detect_parasitic_geneset
#> # A tibble: 112 × 6
#>    Leishmaniasis `Chagas disease` African trypanosomiasi…¹ Malaria Toxoplasmosis
#>    <chr>         <chr>            <chr>                    <chr>   <chr>        
#>  1 IGH           AKT3             IGH                      KLRC4-… AKT3         
#>  2 TAB1          TLR6             IDO2                     COMP    PPIF         
#>  3 FCGR3B        ADCY1            F2RL1                    CR1     LAMC3        
#>  4 CR1           P3R3URF-PIK3R3   PLCB1                    CR1L    TAB1         
#>  5 CR1L          CHUK             GNAQ                     CSF3    CHUK         
#>  6 MAPK14        MAPK14           HBA1                     KLRK1   CCR5         
#>  7 CYBA          TICAM1           HBA2                     ACKR1   MAPK14       
#>  8 CYBB          ACE              HBB                      GYPA    PIK3R6       
#>  9 EEF1A1        AKT1             HPR                      GYPB    AKT1         
#> 10 EEF1A2        AKT2             APOA1                    GYPC    AKT2         
#> # ℹ 102 more rows
#> # ℹ abbreviated name: ¹​`African trypanosomiasis`
#> # ℹ 1 more variable: Amoebiasis <chr>
#> 
#> $Human$Geneset$KEGG$Human_KEGG_Detect_viral_geneset
#> # A tibble: 333 × 21
#>    `Viral life cycle - HIV-1` Virion - Human immunodefi…¹ Virion - Flavivirus …²
#>    <chr>                      <chr>                       <chr>                 
#>  1 PDCD6IP                    CLEC4M                      CLEC4M                
#>  2 APOBEC3A_B                 CCR5                        HAVCR1                
#>  3 CDK9                       CD209                       CD209                 
#>  4 SERINC3                    CXCR4                       PHB1                  
#>  5 CPSF6                      CD4                         MXRA8                 
#>  6 PSIP1                      NA                          TYRO3                 
#>  7 CCR5                       NA                          NA                    
#>  8 CHMP4B                     NA                          NA                    
#>  9 CREBBP                     NA                          NA                    
#> 10 APOBEC3D                   NA                          NA                    
#> # ℹ 323 more rows
#> # ℹ abbreviated names: ¹​`Virion - Human immunodeficiency virus`,
#> #   ²​`Virion - Flavivirus and Alphavirus`
#> # ℹ 18 more variables:
#> #   `Virion - Ebolavirus, Lyssavirus and Morbillivirus` <chr>,
#> #   `Virion - Herpesvirus` <chr>, `Virion - Adenovirus` <chr>,
#> #   `Virion - Rotavirus` <chr>, `Virion - Hepatitis viruses` <chr>, …
#> 
#> $Human$Geneset$KEGG$Human_KEGG_Metabolism_geneset
#> # A tibble: 136 × 85
#>    `Glycolysis / Gluconeogenesis` Citrate cycle (TCA cy…¹ Pentose phosphate pa…²
#>    <chr>                          <chr>                   <chr>                 
#>  1 AKR1A1                         CS                      GLYCTK                
#>  2 ADH1A                          DLAT                    FBP1                  
#>  3 ADH1B                          DLD                     PRPS1L1               
#>  4 ADH1C                          DLST                    ALDOA                 
#>  5 ADH4                           FH                      ALDOB                 
#>  6 ADH5                           IDH1                    RPIA                  
#>  7 ADH6                           IDH2                    ALDOC                 
#>  8 GALM                           IDH3A                   SHPK                  
#>  9 ADH7                           IDH3B                   G6PD                  
#> 10 LDHAL6A                        IDH3G                   PGLS                  
#> # ℹ 126 more rows
#> # ℹ abbreviated names: ¹​`Citrate cycle (TCA cycle)`,
#> #   ²​`Pentose phosphate pathway`
#> # ℹ 82 more variables: `Pentose and glucuronate interconversions` <chr>,
#> #   `Fructose and mannose metabolism` <chr>, `Galactose metabolism` <chr>,
#> #   `Ascorbate and aldarate metabolism` <chr>, `Fatty acid biosynthesis` <chr>,
#> #   `Fatty acid elongation` <chr>, `Fatty acid degradation` <chr>, …
#> 
#> 
#> $Human$Geneset$Reactome
#> $Human$Geneset$Reactome$Human_Reactome_Metabolism_geneset
#> # A tibble: 374 × 301
#>    Mitochondrial iron-sulfur clu…¹ The citric acid (TCA…² Reversible hydration…³
#>    <chr>                           <chr>                  <chr>                 
#>  1 FDX2                            LRPPRC                 CA5B                  
#>  2 ISCA2                           TRAP1                  CA14                  
#>  3 HSCB                            ATP5PD                 CA13                  
#>  4 FDX1                            ATP5MG                 CA1                   
#>  5 FDXR                            ME3                    CA2                   
#>  6 ISCU                            UQCR11                 CA3                   
#>  7 FXN                             COX20                  CA4                   
#>  8 GLRX5                           NDUFA11                CA5A                  
#>  9 SLC25A37                        COX4I1                 CA6                   
#> 10 LYRM4                           COX5B                  CA7                   
#> # ℹ 364 more rows
#> # ℹ abbreviated names: ¹​`Mitochondrial iron-sulfur cluster biogenesis`,
#> #   ²​`The citric acid (TCA) cycle and respiratory electron transport`,
#> #   ³​`Reversible hydration of carbon dioxide`
#> # ℹ 298 more variables: `Inositol phosphate metabolism` <chr>,
#> #   `Metabolism of nucleotides` <chr>,
#> #   `Integration of energy metabolism` <chr>, …
#> 
#> 
#> 
#> 
#> $Mouse
#> $Mouse$Gene
#> $Mouse$Gene$Mouse_Infection_Immunity_gene
#> # A tibble: 85 × 17
#>    IFN      IRF    IFI   CXCL  CXCR  CCL   CCR   TNF   CSF   Il    NFKB  Attacin
#>    <chr>    <chr>  <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr>  
#>  1 Ifna-ps1 Irf1   Ifit… Cxcl1 Cxcr1 Ccl1  Ccr1  Lta   Csf1  Il1rn Nfkb1 Camp   
#>  2 Ifnab    Irf2   Ifit… Cxcl2 Cxcr2 Ccl2  Ccr1… Tnf   Csf1r Il1a  Nfkb2 Defa1  
#>  3 Ifne     Irf2b… Ifit1 Cxcl3 Cxcr3 Ccl3  Ccrl2 Ltb   Csf1… Il1b  Nfkb… Defa2  
#>  4 Ifng     Irf2b… Ifit… Cxcl5 Cxcr4 Ccl4  Ccr2  Tnfs… Csf2  Il1b… Nfkb… Defa3  
#>  5 Ifnk     Irf2b… Ifih1 Cxcl9 Cxcr5 Ccl5  Ccr3  Cd40… Csf2… Il1f… Nfkb… Defa5  
#>  6 Ifnz     Irf3   Ifit… Cxcl… Cxcr6 Ccl6  Ccr4  Faslg Csf2… Il1r1 Nfkb… Defb1  
#>  7 Ifna1    Irf4   Ifit2 Cxcl… NA    Ccl7  Ccr5  Cd70  Csf2… Il1r2 Nfkb… Defb2  
#>  8 Ifnar1   Irf5   Ifit3 Cxcl… NA    Ccl8  Ccr6  Tnfs… Csf3r Il1r… Nfkb… Defb3  
#>  9 Ifnb1    Irf6   Ifit… Cxcl… NA    Ccl9  Ccr7  Tnfs… Csf3  Il1r… Rel   Defb4  
#> 10 Ifngas1  Irf7   Ifit… Cxcl… NA    Ccl11 Ccr8  Tnfs… NA    Il1r… Rela  Defb5  
#> # ℹ 75 more rows
#> # ℹ 5 more variables: ATE <chr>, ROS_NO <chr>, Clq <chr>, Acute_Phase <chr>,
#> #   Antigen <chr>
#> 
#> $Mouse$Gene$Mouse_Macrophage_gene
#> # A tibble: 15 × 3
#>    M1     M2     Fusion 
#>    <chr>  <chr>  <chr>  
#>  1 Fcgr3  Cd163  Adam9  
#>  2 Cd86   Mrc1   Cd44   
#>  3 Cd14   Arg1   Cd81   
#>  4 Nos2   Chil3  Dcstamp
#>  5 Il6    Retnla Ocstamp
#>  6 Tnf    Il10   Stat1  
#>  7 Il1b   Tgfb1  Trem2  
#>  8 Il12b  Pparg  Tyrobp 
#>  9 Cxcl9  Klf4   NA     
#> 10 Cxcl10 Il4ra  NA     
#> 11 Cd80   Ccl17  NA     
#> 12 H2-Aa  Ccl22  NA     
#> 13 Irf5   NA     NA     
#> 14 Socs3  NA     NA     
#> 15 Fpr2   NA     NA     
#> 
#> $Mouse$Gene$Mouse_PRR_gene
#> # A tibble: 37 × 8
#>    TLRs  TLRs_downstream CASP     NLRs      CLRs   SRs    CRs   FcRs  
#>    <chr> <chr>           <chr>    <chr>     <chr>  <chr>  <chr> <chr> 
#>  1 Tlr1  Myd88           Casp1    Citta     Mrc1   Msr1   Cr1l  Fcgrt 
#>  2 Tlr2  Ticam1          Casp2    Naip1     Clec1a Scaf1  Cr2   Fcgr1 
#>  3 Tlr3  Ticam2          Casp3    Naip2     Clec1b Cd36   Itgam Fcgr2b
#>  4 Tlr4  Tirap           Casp4    Naip3     Clec2m Scarb1 Itgax Fcgr3 
#>  5 Cd14  Tram1           Casp6    Naip3-ps1 Clec2d Marco  Fcgr1 Fcgr4 
#>  6 Tlr5  Tram2           Casp7    Naip5     Clec2e Cd163  NA    NA    
#>  7 Tlr6  Tram1l1         Casp8    Naip6     Clec2f NA     NA    NA    
#>  8 Tlr7  Traf1           Casp8ap2 Nod1      Clec2g NA     NA    NA    
#>  9 Tlr8  Trafd1          Casp9    Nod2      Clec2h NA     NA    NA    
#> 10 Tlr9  Traf2           Casp12   Nodal     Clec2i NA     NA    NA    
#> # ℹ 27 more rows
#> 
#> 
#> $Mouse$Geneset
#> $Mouse$Geneset$GO
#> $Mouse$Geneset$GO$Mouse_GO_BP_Detect_bacterial_geneset
#> # A tibble: 484 × 7
#>    GOBP_ANTIBACTERIAL_HUMORAL_RE…¹ GOBP_ANTIBACTERIAL_I…² GOBP_ANTIBACTERIAL_P…³
#>    <chr>                           <chr>                  <chr>                 
#>  1 AY761185                        Defb42                 Elane                 
#>  2 Adm                             Il33                   Evpl                  
#>  3 Ang                             Marchf2                Ivl                   
#>  4 Ang2                            Mpeg1                  Klk5                  
#>  5 Ang4                            Muc19                  Klk7                  
#>  6 Ang5                            Nfkb1                  Lgals4                
#>  7 Ang6                            Nod2                   Mmp7                  
#>  8 App                             Rasgrp4                Nod2                  
#>  9 B2m                             Slc15a2                Pgc                   
#> 10 Bpi                             Tfeb                   Ppl                   
#> # ℹ 474 more rows
#> # ℹ abbreviated names: ¹​GOBP_ANTIBACTERIAL_HUMORAL_RESPONSE,
#> #   ²​GOBP_ANTIBACTERIAL_INNATE_IMMUNE_RESPONSE,
#> #   ³​GOBP_ANTIBACTERIAL_PEPTIDE_PRODUCTION
#> # ℹ 4 more variables: GOBP_DETECTION_OF_MOLECULE_OF_BACTERIAL_ORIGIN <chr>,
#> #   GOBP_REGULATION_OF_ANTIBACTERIAL_PEPTIDE_PRODUCTION <chr>,
#> #   GOBP_RESPONSE_TO_BACTERIAL_LIPOPROTEIN <chr>, …
#> 
#> $Mouse$Geneset$GO$Mouse_GO_BP_Detect_viral_geneset
#> # A tibble: 367 × 12
#>    GOBP_CELLULAR_RESPONSE_TO_VIRUS GOBP_DEFENSE_RESPONS…¹ GOBP_DETECTION_OF_VI…²
#>    <chr>                           <chr>                  <chr>                 
#>  1 Adar                            Abcc9                  Ncr1                  
#>  2 Arf1                            Abcf3                  Rigi                  
#>  3 Atf2                            Acod1                  Serinc3               
#>  4 Bax                             Agbl4                  Serinc5               
#>  5 Ccl19                           Agbl5                  Zcchc3                
#>  6 Ccl19-ps1                       Aicda                  NA                    
#>  7 Ccl19-ps3                       Aim2                   NA                    
#>  8 Ccl19-ps4                       Akap1                  NA                    
#>  9 Ccl19-ps5                       Apobec1                NA                    
#> 10 Ccl19-ps6                       Apobec3                NA                    
#> # ℹ 357 more rows
#> # ℹ abbreviated names: ¹​GOBP_DEFENSE_RESPONSE_TO_VIRUS,
#> #   ²​GOBP_DETECTION_OF_VIRUS
#> # ℹ 9 more variables: GOBP_INTRACELLULAR_TRANSPORT_OF_VIRUS <chr>,
#> #   GOBP_NEGATIVE_REGULATION_OF_DEFENSE_RESPONSE_TO_VIRUS <chr>,
#> #   GOBP_POSITIVE_REGULATION_OF_DEFENSE_RESPONSE_TO_VIRUS_BY_HOST <chr>,
#> #   GOBP_RECEPTOR_MEDIATED_ENDOCYTOSIS_OF_VIRUS_BY_HOST_CELL <chr>, …
#> 
#> 
#> $Mouse$Geneset$KEGG
#> $Mouse$Geneset$KEGG$Mouse_KEGG_Detect_bacterial_geneset
#> # A tibble: 252 × 7
#>    Bacterial invasion of epithe…¹ `Salmonella infection` Pertussis Legionellosis
#>    <chr>                          <chr>                  <chr>     <chr>        
#>  1 Dnm3                           Dynlt1f                Ticam1    Apaf1        
#>  2 Actb                           Dynlt1c                Nod1      Arf1         
#>  3 Actg1                          Ahnak2                 Tirap     Arf2         
#>  4 Rhoa                           Dynlt1a                Rhoa      Bnip3        
#>  5 Arpc1b                         Brk1                   Serping1  C3           
#>  6 Ctnna1                         Exoc5                  C1qa      Casp1        
#>  7 Ctnna2                         Nckap1l                C1qb      Casp3        
#>  8 Ctnnb1                         Nod1                   C1qc      Casp7        
#>  9 Cav1                           Raf1                   C2        Casp8        
#> 10 Cav2                           Dync2h1                C3        Casp9        
#> # ℹ 242 more rows
#> # ℹ abbreviated name: ¹​`Bacterial invasion of epithelial cells`
#> # ℹ 3 more variables: `Yersinia infection` <chr>,
#> #   `Staphylococcus aureus infection` <chr>, Tuberculosis <chr>
#> 
#> $Mouse$Geneset$KEGG$Mouse_KEGG_Detect_parasitic_geneset
#> # A tibble: 109 × 6
#>    Leishmaniasis `Chagas disease` African trypanosomiasi…¹ Malaria Toxoplasmosis
#>    <chr>         <chr>            <chr>                    <chr>   <chr>        
#>  1 H2-Ea         Ticam1           Hbb-bs                   Klrb1   H2-Ea        
#>  2 C3            Calr4            Hbb-bt                   Hbb-bs  Pik3r6       
#>  3 Cr1l          Ace              Hba-a2                   Hbb-bt  Ppif         
#>  4 Cyba          Akt1             Apoa1                    Hba-a2  Akt1         
#>  5 Cybb          Akt2             F2rl1                    Cd36    Akt2         
#>  6 Eef1a1        Bdkrb2           Fas                      Cd81    Alox5        
#>  7 Eef1a2        C1qa             Fasl                     Comp    Birc3        
#>  8 Elk1          C1qb             Gnaq                     Cr1l    Birc2        
#>  9 Fcgr1         C1qc             Hba-a1                   Csf3    Xiap         
#> 10 Fcgr3         C3               Hbb-b1                   Ackr1   Bad          
#> # ℹ 99 more rows
#> # ℹ abbreviated name: ¹​`African trypanosomiasis`
#> # ℹ 1 more variable: Amoebiasis <chr>
#> 
#> $Mouse$Geneset$KEGG$Mouse_KEGG_Detect_viral_geneset
#> # A tibble: 359 × 21
#>    `Viral life cycle - HIV-1` Virion - Human immunodefi…¹ Virion - Flavivirus …²
#>    <chr>                      <chr>                       <chr>                 
#>  1 Supt4b                     Cd4                         Cd209c                
#>  2 Psip1                      Cxcr4                       Cd209d                
#>  3 Xpo1                       Ccr5                        Cd209e                
#>  4 Cdk9                       Cd209c                      Cd209a                
#>  5 Vps4a                      Cd209d                      Havcr1                
#>  6 Bicd1                      Cd209e                      Timd2                 
#>  7 Ccnt1                      Cd209a                      Phb1                  
#>  8 Cd4                        Cd209f                      Tyro3                 
#>  9 Cxcr4                      Cd209b                      1700071K01Rik         
#> 10 Ccr5                       Cd209g                      Cd209f                
#> # ℹ 349 more rows
#> # ℹ abbreviated names: ¹​`Virion - Human immunodeficiency virus`,
#> #   ²​`Virion - Flavivirus and Alphavirus`
#> # ℹ 18 more variables:
#> #   `Virion - Ebolavirus, Lyssavirus and Morbillivirus` <chr>,
#> #   `Virion - Herpesvirus` <chr>, `Virion - Adenovirus` <chr>,
#> #   `Virion - Rotavirus` <chr>, `Virion - Hepatitis viruses` <chr>, …
#> 
#> $Mouse$Geneset$KEGG$Mouse_KEGG_Metabolism_geneset
#> # A tibble: 141 × 85
#>    `Glycolysis / Gluconeogenesis` Citrate cycle (TCA cy…¹ Pentose phosphate pa…²
#>    <chr>                          <chr>                   <chr>                 
#>  1 Gck                            Acly                    H6pd                  
#>  2 Ldhal6b                        Aco1                    Pgd                   
#>  3 Aldh7a1                        Aco2                    Prps2                 
#>  4 Adh1                           Cs                      Aldoa                 
#>  5 Adh7                           Dld                     Aldoc                 
#>  6 Adh5                           Fh1                     Fbp2                  
#>  7 Gm29667                        Idh1                    Fbp1                  
#>  8 Aldh2                          Idh3g                   G6pd2                 
#>  9 Aldh3a1                        Idh3b                   G6pdx                 
#> 10 Aldh3a2                        Mdh2                    Gpi1                  
#> # ℹ 131 more rows
#> # ℹ abbreviated names: ¹​`Citrate cycle (TCA cycle)`,
#> #   ²​`Pentose phosphate pathway`
#> # ℹ 82 more variables: `Pentose and glucuronate interconversions` <chr>,
#> #   `Fructose and mannose metabolism` <chr>, `Galactose metabolism` <chr>,
#> #   `Ascorbate and aldarate metabolism` <chr>, `Fatty acid biosynthesis` <chr>,
#> #   `Fatty acid elongation` <chr>, `Fatty acid degradation` <chr>, …
#> 
#> 
#> $Mouse$Geneset$Mouse_HotSpot_geneset
#> # A tibble: 0 × 2
#> # ℹ 2 variables: Inflammatory <chr>, Immunity <chr>
#> 
#> $Mouse$Geneset$Mouse_MSigDB_Hallmark_geneset
#> # A tibble: 200 × 50
#>    HALLMARK_ADIPOGENESIS HALLMARK_ALLOGRAFT_REJECTION HALLMARK_ANDROGEN_RESPONSE
#>    <chr>                 <chr>                        <chr>                     
#>  1 Abca1                 Aars1                        Abcc4                     
#>  2 Abcb8                 Abce1                        Abhd2                     
#>  3 Acaa2                 Abi1                         Acsl3                     
#>  4 Acadl                 Ache                         Actn1                     
#>  5 Acadm                 Acvr2a                       Adamts1                   
#>  6 Acads                 Akt1                         Adrm1                     
#>  7 Acly                  Apbb1                        Akap12                    
#>  8 Aco2                  B2m                          Akt1                      
#>  9 Acox1                 Bcat1                        Aldh1a3                   
#> 10 Adcy6                 Bcl10                        Ank                       
#> # ℹ 190 more rows
#> # ℹ 47 more variables: HALLMARK_ANGIOGENESIS <chr>,
#> #   HALLMARK_APICAL_JUNCTION <chr>, HALLMARK_APICAL_SURFACE <chr>,
#> #   HALLMARK_APOPTOSIS <chr>, HALLMARK_BILE_ACID_METABOLISM <chr>,
#> #   HALLMARK_CHOLESTEROL_HOMEOSTASIS <chr>, HALLMARK_COAGULATION <chr>,
#> #   HALLMARK_COMPLEMENT <chr>, HALLMARK_DNA_REPAIR <chr>,
#> #   HALLMARK_E2F_TARGETS <chr>, …
#> 
#> $Mouse$Geneset$Mouse_PCD_geneset
#> # A tibble: 577 × 18
#>    Apoptosis Pyroptosis Ferroptosis Autophagy Necroptosis Cuproptosis
#>    <chr>     <chr>      <chr>       <chr>     <chr>       <chr>      
#>  1 Cd74      Gsdmc      Sqle        Prkaa1    Dnm1l       Dlst       
#>  2 Zfp385a   Gsdmc2     Fdft1       Snx6      Xiap        Atp7a      
#>  3 Slc25a5   Gsdmc3     Abcc1       Srebf2    Ifngr2      Dbt        
#>  4 Sod1      Gsdmc4     Lpcat3      Gsk3b     Ifnar1      Atp7b      
#>  5 Gsk3b     Gm10053    Ftl1-ps1    Hsp90aa1  Ifnar2      Pdha1      
#>  6 Fxn       Chmp7      Sat1        Trim27    Tlr3        Dlat       
#>  7 Ppp1r13b  Chmp3      Gclm        Atp6v1a   Il33        Fdx1       
#>  8 Appl1     Scaf11     Nox1        Pik3c3    Ftl1-ps1    Slc31a1    
#>  9 Snai2     Tirap      Slc39a14    Gata4     Ifng        Nlrp3      
#> 10 Mir17     Gsdma3     Zeb1        Tspo      Aifm1       Gcsh       
#> # ℹ 567 more rows
#> # ℹ 12 more variables: Parthanatos <chr>, Entoticcelldeath <chr>,
#> #   Netoticcelldeath <chr>, Lysosome.dependentcelldeath <chr>,
#> #   Alkaliptosis <chr>, Oxeiptosis <chr>, NETosis <chr>,
#> #   Immunogenic_cell_death <chr>, Anoikis <chr>, Paraptosis <chr>,
#> #   Methuosis <chr>, Entosis <chr>
#> 
#> $Mouse$Geneset$Reactome
#> $Mouse$Geneset$Reactome$Mouse_Reactome_Metabolism_geneset
#> # A tibble: 260 × 253
#>    Mitochondrial iron-sulfur clu…¹ The citric acid (TCA…² Reversible hydration…³
#>    <chr>                           <chr>                  <chr>                 
#>  1 Hscb                            Ndufb4b                Car1                  
#>  2 Fdx1                            Sco2                   Car2                  
#>  3 Fdxr                            Ndufb11                Car3                  
#>  4 Fxn                             Ldhal6b                Car4                  
#>  5 Nfs1                            Me2                    Car5a                 
#>  6 Lyrm4                           Me3                    Car6                  
#>  7 Iscu                            Glo1                   Car7                  
#>  8 Fdx2                            Cox6b1                 Car9                  
#>  9 Isca1                           Etfb                   Car14                 
#> 10 Glrx5                           Etfa                   Car5b                 
#> # ℹ 250 more rows
#> # ℹ abbreviated names: ¹​`Mitochondrial iron-sulfur cluster biogenesis`,
#> #   ²​`The citric acid (TCA) cycle and respiratory electron transport`,
#> #   ³​`Reversible hydration of carbon dioxide`
#> # ℹ 250 more variables: `Inositol phosphate metabolism` <chr>,
#> #   `Metabolism of nucleotides` <chr>,
#> #   `Integration of energy metabolism` <chr>, …
```

## Session information

``` r
sessionInfo()
#> R version 4.2.0 (2022-04-22 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 22000)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
#> [2] LC_CTYPE=Chinese (Simplified)_China.utf8   
#> [3] LC_MONETARY=Chinese (Simplified)_China.utf8
#> [4] LC_NUMERIC=C                               
#> [5] LC_TIME=Chinese (Simplified)_China.utf8    
#> 
#> attached base packages:
#> [1] tcltk     stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] DynDoc_1.74.0       Biobase_2.58.0      BiocGenerics_0.44.0
#> [4] widgetTools_1.74.0 
#> 
#> loaded via a namespace (and not attached):
#>   [1] ggvenn_0.1.19               ica_1.0-3                  
#>   [3] ggdensity_1.0.0             foreach_1.5.2              
#>   [5] lmtest_0.9-40               crayon_1.5.3               
#>   [7] MASS_7.3-60                 nlme_3.1-162               
#>   [9] GOSemSim_2.24.0             rlang_1.1.7                
#>  [11] XVector_0.38.0              HDO.db_0.99.1              
#>  [13] ROCR_1.0-11                 irlba_2.3.5.1              
#>  [15] proto_1.0.0                 BiocParallel_1.30.4        
#>  [17] rjson_0.2.21                bit64_4.0.5                
#>  [19] glue_1.8.0                  sctransform_0.4.1          
#>  [21] parallel_4.2.0              spatstat.sparse_3.0-3      
#>  [23] AnnotationDbi_1.60.2        SeuratDisk_0.0.0.9021      
#>  [25] dotCall64_1.1-1             DOSE_3.24.2                
#>  [27] spatstat.geom_3.2-9         tidyselect_1.2.1           
#>  [29] SummarizedExperiment_1.28.0 distances_0.1.10           
#>  [31] UCell_2.0.1                 SeuratObject_5.3.0         
#>  [33] fitdistrplus_1.1-11         XML_3.99-0.16.1            
#>  [35] tidyr_1.3.1                 zoo_1.8-12                 
#>  [37] chron_2.3-61                xtable_1.8-4               
#>  [39] RcppHNSW_0.6.0              magrittr_2.0.3             
#>  [41] evaluate_1.0.1              ggplot2_4.0.2              
#>  [43] cli_3.6.5                   zlibbioc_1.44.0            
#>  [45] dbscan_1.1-11               rstudioapi_0.15.0          
#>  [47] miniUI_0.1.1.1              furrr_0.3.1                
#>  [49] sp_2.1-3                    bslib_0.8.0                
#>  [51] bmp_0.3                     fastmatch_1.1-3            
#>  [53] treeio_1.22.0               fastDummies_1.7.3          
#>  [55] shiny_1.13.0                xfun_0.49                  
#>  [57] gson_0.1.0                  cluster_2.1.4              
#>  [59] caTools_1.18.2              tidygraph_1.2.3            
#>  [61] tkWidgets_1.74.0            KEGGREST_1.36.0            
#>  [63] tibble_3.2.1                ggrepel_0.9.5              
#>  [65] ape_5.8                     listenv_0.9.1              
#>  [67] Biostrings_2.66.0           png_0.1-8                  
#>  [69] future_1.33.1               withr_3.0.2                
#>  [71] bitops_1.0-7                ggforce_0.4.1              
#>  [73] plyr_1.8.9                  GSEABase_1.60.0            
#>  [75] pillar_1.9.0                gplots_3.1.3.1             
#>  [77] cachem_1.1.0                otel_0.2.0                 
#>  [79] fs_1.6.3                    hdf5r_1.3.9                
#>  [81] clusterProfiler_4.6.2       hash_2.2.6.3               
#>  [83] paletteer_1.5.0             DelayedMatrixStats_1.20.0  
#>  [85] vctrs_0.6.5                 generics_0.1.3             
#>  [87] gsubfn_0.7                  tools_4.2.0                
#>  [89] tweenr_2.0.2                fgsea_1.22.0               
#>  [91] proxy_0.4-27                DelayedArray_0.24.0        
#>  [93] fastmap_1.2.0               compiler_4.2.0             
#>  [95] abind_1.4-8                 httpuv_1.6.12              
#>  [97] plotly_4.10.4               GenomeInfoDbData_1.2.9     
#>  [99] gridExtra_2.3               lattice_0.21-8             
#> [101] ggnewscale_0.4.9            deldir_2.0-4               
#> [103] utf8_1.2.4                  later_1.3.2                
#> [105] dplyr_1.1.4                 STRINGdb_2.18.0            
#> [107] jsonlite_1.8.8              concaveman_1.2.0           
#> [109] scales_1.4.0                graph_1.76.0               
#> [111] tidytree_0.4.4              pbapply_1.7-2              
#> [113] sparseMatrixStats_1.10.0    lazyeval_0.2.2             
#> [115] promises_1.5.0              doParallel_1.0.17          
#> [117] R.utils_2.12.3              goftest_1.2-3              
#> [119] spatstat.utils_3.0-4        reticulate_1.45.0          
#> [121] rmarkdown_2.29              pkgdown_2.2.0              
#> [123] cowplot_1.1.3               textshaping_0.3.6          
#> [125] Rtsne_0.17                  dichromat_2.0-0.1          
#> [127] downloader_0.4              uwot_0.2.2                 
#> [129] igraph_2.0.3                survival_3.5-5             
#> [131] yaml_2.3.10                 plotrix_3.8-2              
#> [133] systemfonts_1.0.4           htmltools_0.5.8.1          
#> [135] memoise_2.0.1               Seurat_5.4.0               
#> [137] graphlayouts_1.1.1          IRanges_2.32.0             
#> [139] viridisLite_0.4.2           digest_0.6.35              
#> [141] assertthat_0.2.1            mime_0.12                  
#> [143] rappdirs_0.3.3              tiff_0.1-12                
#> [145] spam_2.10-0                 RSQLite_2.3.6              
#> [147] sqldf_0.4-11                yulab.utils_0.2.4          
#> [149] future.apply_1.11.3         data.table_1.15.4          
#> [151] blob_1.2.4                  S4Vectors_0.36.2           
#> [153] R.oo_1.26.0                 ragg_1.3.0                 
#> [155] splines_4.2.0               rematch2_2.1.2             
#> [157] RCurl_1.98-1.14             colorspace_2.1-0           
#> [159] mnormt_2.1.1                GenomicRanges_1.50.2       
#> [161] aplot_0.1.10                sass_0.4.9                 
#> [163] Rcpp_1.0.11                 RANN_2.6.1                 
#> [165] ggh4x_0.2.8                 enrichplot_1.18.3          
#> [167] fansi_1.0.6                 parallelly_1.37.1          
#> [169] R6_2.6.1                    grid_4.2.0                 
#> [171] ggridges_0.5.6              lifecycle_1.0.5            
#> [173] zip_2.3.0                   rBLAST_1.3.1               
#> [175] ggsignif_0.6.4              jquerylib_0.1.4            
#> [177] Matrix_1.6-5                qvalue_2.30.0              
#> [179] desc_1.4.3                  RcppAnnoy_0.0.22           
#> [181] RColorBrewer_1.1-3          iterators_1.0.14           
#> [183] spatstat.explore_3.2-7      stringr_1.5.1              
#> [185] S7_0.2.0                    htmlwidgets_1.6.4          
#> [187] Mfuzz_2.56.0                polyclip_1.10-6            
#> [189] purrr_1.0.2                 shadowtext_0.1.2           
#> [191] gridGraphics_0.5-1          globals_0.16.3             
#> [193] readbitmap_0.1.5            patchwork_1.3.0            
#> [195] spatstat.random_3.2-3       progressr_0.14.0           
#> [197] codetools_0.2-19            matrixStats_1.1.0          
#> [199] GO.db_3.16.0                FNN_1.1.4.1                
#> [201] gtools_3.9.5                psych_2.4.3                
#> [203] SingleCellExperiment_1.20.1 RSpectra_0.16-1            
#> [205] R.methodsS3_1.8.2           GenomeInfoDb_1.34.9        
#> [207] imager_0.45.8               gtable_0.3.6               
#> [209] DBI_1.2.2                   stats4_4.2.0               
#> [211] ggfun_0.1.1                 tensor_1.5                 
#> [213] httr_1.4.7                  KernSmooth_2.23-22         
#> [215] stringi_1.8.3               anndata_0.7.5.6            
#> [217] reshape2_1.4.4              farver_2.1.2               
#> [219] annotate_1.76.0             viridis_0.6.4              
#> [221] ggtree_3.6.2                AUCell_1.28.0              
#> [223] ggplotify_0.1.1             scattermore_1.2            
#> [225] bit_4.0.5                   shinycssloaders_1.1.0      
#> [227] scatterpie_0.2.1            jpeg_0.1-10                
#> [229] MatrixGenerics_1.10.0       spatstat.data_3.0-4        
#> [231] ggraph_2.1.0                pkgconfig_2.0.3            
#> [233] ggprism_1.0.5               STID_0.0.0.9000            
#> [235] knitr_1.49
```
