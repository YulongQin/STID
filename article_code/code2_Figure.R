

{
  source("./code0_library.R")
  source("./STID/R/Utilities.R")
  source("./STID/R/STID-class.R")
  source("./STID/R/Preprocessing.R")
  source("./STID/R/SpotDetect.R")
  source("./STID/R/NicheDetect.R")
  source("./STID/R/SingleSampAnalysis.R")
  source("./STID/R/MultiSampAnalysis.R")
  source("./STID/R/RPyCall.R")
  source("./STID/R/zzz.R")
}
# 先加载R包和代码后，再导入rds，否则会报错

#### 一、Figure2 ####
###  1.A-D：CE ####
# 运行得到CE的STID_obj

### 1.1 CorrectBackgroud
options("parallel_workers" = 1)
STID_obj_before <- STID_obj
STID_obj_after <- CorrectBackgroud(STID_obj = STID_obj_before, bg_samp_id = c("DPI_0_1"), bg_features = pathogen_genes, 
                                   PosThres_prob = 0.95, 
                                   assay_id = "Spatial", layer_id = "counts",
                                   grp_nm = "Final_CE_0.95", dir_nm = "M1_CorrectBackgroud")


### 1.2 SpotDetect_Gene
high_exp_genes <- GetTopGenes(STID_obj_before, top_n = 10, pattern = "EmuJ-",
                              grp_by_samp = F, grp_by_celltype = F,
                              assay_id = "Spatial", layer_id = "counts")
STID_obj_before <- GetGeneStat(STID_obj = STID_obj_before, features = high_exp_genes,prefix = "top10_gene",func = "sum") %>% 
  AddMetaColumn(STID_obj = STID_obj_before,
                add_data = ., # data.frame
                meta_key = "raw", # string
                igrnore_rownm = FALSE)
STID_obj_before <- GetGeneStat(STID_obj = STID_obj_before, features = pathogen_genes,prefix = "all_gene",func = "sum") %>%
  AddMetaColumn(STID_obj = STID_obj_before,
                add_data = ., # data.frame
                meta_key = "raw",
                igrnore_rownm = FALSE)
STID_obj_after <- GetGeneStat(STID_obj = STID_obj_after, features = high_exp_genes,prefix = "top10_gene",func = "sum") %>%
  AddMetaColumn(STID_obj = STID_obj_after,
                add_data = ., # data.frame
                meta_key = "raw",
                igrnore_rownm = FALSE)
STID_obj_after <- GetGeneStat(STID_obj = STID_obj_after, features = pathogen_genes,prefix = "all_gene",func = "sum") %>%
  AddMetaColumn(STID_obj = STID_obj_after,
                add_data = ., # data.frame
                meta_key = "raw",
                igrnore_rownm = FALSE)
meta_data <- STID_obj_before@meta.data

# >仅在Figure1中有效
STID_obj_before <- SpotDetect_Gene(STID_obj_before,
                                   features = high_exp_genes,
                                   feature_colnm = grep("top10_gene",colnames(STID_obj_before@meta.data),value = T),
                                   PosThres_prob = 0, PosThres_count = 1,
                                   # col = c("#3D3576","#92D74D"),
                                   # col = c("#8DBDDC","#FD9F8F"),
                                   col = c("grey20","#80FFFF"),
                                   black_bg = T,pt_size = 1,
                                   # blur_method = "isoblur",
                                   blur_method = NULL,
                                   blur_n = 1,blur_sigma = 0.5, 
                                   plot_method = "single",
                                   grp_nm = "CE_correct_before_top10_gene_black")
tmp <- GetSapThreshold(STID_obj_before,
                       meta_key = "M1_SpotDetect_Gene_CE_correct_before_top10_gene_black")
STID_obj_before <- SpotDetect_Gene(STID_obj_before,
                                   features = high_exp_genes,
                                   feature_colnm = grep("all_gene",colnames(STID_obj_before@meta.data),value = T),
                                   PosThres_prob = 0, PosThres_count = 1,
                                   # col = c("#3D3576","#92D74D"),
                                   # col = c("#8DBDDC","#FD9F8F"),
                                   col = c("grey20","#80FFFF"),
                                   black_bg = T,pt_size = 1,
                                   # blur_method = "isoblur",
                                   blur_method = NULL,
                                   blur_n = 1,blur_sigma = 0.5, 
                                   plot_method = "single",
                                   grp_nm = "CE_correct_before_all_gene_black")
STID_obj_after <- SpotDetect_Gene(STID_obj_after,
                                  features = high_exp_genes,
                                  feature_colnm = grep("top10_gene",colnames(STID_obj_after@meta.data),value = T),
                                  PosThres_prob = 0, PosThres_count = 1,
                                  # col = c("#3D3576","#92D74D"),
                                  # col = c("#8DBDDC","#FD9F8F"),
                                  col = c("grey20","#80FFFF"),
                                  black_bg = T,pt_size = 1,
                                  # blur_method = "isoblur",
                                  blur_method = NULL,
                                  blur_n = 1,blur_sigma = 0.5, 
                                  plot_method = "single",
                                  grp_nm = "CE_correct_after_top10_gene_black")
STID_obj_after <- SpotDetect_Gene(STID_obj_after,
                                  features = high_exp_genes,
                                  feature_colnm = grep("all_gene",colnames(STID_obj_after@meta.data),value = T),
                                  PosThres_prob = 0, PosThres_count = 1,
                                  # col = c("#3D3576","#92D74D"),
                                  # col = c("#8DBDDC","#FD9F8F"),
                                  col = c("grey20","#80FFFF"),
                                  black_bg = T,pt_size = 1,
                                  # blur_method = "isoblur",
                                  blur_method = NULL,
                                  blur_n = 1,blur_sigma = 0.5, 
                                  plot_method = "single",
                                  grp_nm = "CE_correct_after_all_gene_black")


### 2.E-K: CE ####
## 2.1 CE ####
# gene
STID_obj_after <- SpotDetect_Gene(STID_obj_after,
                                  features = high_exp_genes,
                                  feature_colnm = grep("all_gene",colnames(STID_obj_after@meta.data),value = T),
                                  PosThres_prob = 0, PosThres_count = 3,
                                  # col = c("#3D3576","#92D74D"),
                                  # col = c("#8DBDDC","#FD9F8F"),
                                  col = COLOR_DIS_CON,
                                  black_bg = F,pt_size = 1,
                                  # blur_method = "isoblur",
                                  blur_method = NULL,
                                  blur_n = 1,blur_sigma = 0.5, 
                                  plot_method = "single",
                                  grp_nm = "CE_correct_after_all_gene_white")
STID_obj_after <- SpotDetect_Gene(STID_obj_after,
                                  features = host_genes,
                                  feature_colnm = grep("all_gene",colnames(STID_obj_after@meta.data),value = T),
                                  PosThres_prob = 0, PosThres_count = 3, # D4_2为3，D15_2为10
                                  # col = c("grey95","red"),
                                  col = COLOR_DIS_CON,
                                  # col = c("grey95","#6A1B9A"),
                                  # col = c("grey95","#6A1B9A"),
                                  # col = c("grey95","#512DA8"),
                                  black_bg = F,pt_size = 1,
                                  blur_method = NULL, blur_n = 1,blur_sigma = 0.5, 
                                  plot_method = "single",
                                  grp_nm = "CE_correct_after_host_gene_white")

# geneset
geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/Mouse_PCD_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
geneset_list <- geneset_list[5]
names(geneset_list)
STID_obj_after <- SpotDetect_Geneset(STID_obj_after,
                                     geneset_list = geneset_list,
                                     score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                                     PosThres_prob = 0.25, PosThres_score = 0, # ！！！！调低PosThres_prob？
                                     pt_size = 1,
                                     col = COLOR_DIS_CON,
                                     black_bg = F, blur_method = NULL,
                                     plot_method = "single",
                                     grp_nm = "CE_correct_after_PCD_white")

geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/KEGG/Mouse_KEGG_Detect_parasitic_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
geneset_list <- geneset_list[4]
names(geneset_list)
STID_obj_after <- SpotDetect_Geneset(STID_obj_after,
                                     geneset_list = geneset_list,
                                     score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                                     PosThres_prob = 0.25, PosThres_score = 0,
                                     pt_size = 1,
                                     col = COLOR_DIS_CON,
                                     black_bg = F, blur_method = NULL,
                                     plot_method = "single",
                                     grp_nm = "CE_correct_after_KEGG_Parasite_white")
saveRDS(STID_obj_after,file = "./rds/STID_obj_CE_SpotDetect.rds")


## 2.2 TB ####
# 需要先做背景校正Supp Figure1: A-D：TB 

# gene
STID_obj_after <- SpotDetect_Gene(STID_obj_after,
                                  features = high_exp_genes,
                                  feature_colnm = grep("all_gene",colnames(STID_obj_after@meta.data),value = T),
                                  PosThres_prob = 0, PosThres_count = 0,
                                  # col = c("#3D3576","#92D74D"),
                                  # col = c("#8DBDDC","#FD9F8F"),
                                  col = COLOR_DIS_CON,
                                  black_bg = F,pt_size = 0.5,
                                  # blur_method = "isoblur",
                                  blur_method = NULL,
                                  blur_n = 1,blur_sigma = 0.5, 
                                  plot_method = "single",
                                  grp_nm = "TB_correct_after_all_gene_white")
STID_obj_after <- SpotDetect_Gene(STID_obj_after,
                                  features = host_genes,
                                  feature_colnm = grep("all_gene",colnames(STID_obj_after@meta.data),value = T),
                                  PosThres_prob = 0, PosThres_count = 0,
                                  # col = c("grey95","red"),
                                  col = COLOR_DIS_CON,
                                  # col = c("grey95","#6A1B9A"),
                                  # col = c("grey95","#6A1B9A"),
                                  # col = c("grey95","#512DA8"),
                                  black_bg = F,pt_size = 0.5,
                                  blur_method = NULL, blur_n = 1,blur_sigma = 0.5, 
                                  plot_method = "single",
                                  grp_nm = "TB_correct_after_host_gene_white")

# geneset
geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/Mouse_PCD_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
geneset_df <- geneset_df[10]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj_after <- SpotDetect_Geneset(STID_obj_after,
                                     geneset_list = geneset_list,
                                     score_method = "AddModuleScore", n_iter = 5, nbin = 24,seed = 100,
                                     PosThres_prob = 0.7, PosThres_score = 0,
                                     pt_size = 0.5,
                                     col = COLOR_DIS_CON,
                                     black_bg = F, blur_method = NULL,
                                     plot_method = "single",
                                     grp_nm = "TB_correct_after_PCD_white")

geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/KEGG/Mouse_KEGG_Detect_bacterial_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
geneset_df <- geneset_df[7]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj_after <- SpotDetect_Geneset(STID_obj_after,
                                     geneset_list = geneset_list,
                                     score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                                     PosThres_prob = 0.7, PosThres_score = 0,
                                     pt_size = 0.5,
                                     col = COLOR_DIS_CON,
                                     black_bg = F, blur_method = NULL,
                                     plot_method = "single",
                                     grp_nm = "TB_correct_after_KEGG_bacterial_white")
STID_obj_after


## 2.3 JEV：D5_1/D3_1 ####
# gene
pathogen_genes <- c("NS5","C","NS3","NS1","E","Prm","NS4aAlt","NS4bAlt","NS2a","NS2b") # 只有包含至少两个 JEV 基因或三个病毒读数的 bin 才被视为 JEV 阳性
STID_obj <- GetGeneStat(STID_obj = STID_obj, features = pathogen_genes,prefix = "all_gene",func = "sum") %>% 
  AddMetaColumn(STID_obj = STID_obj,
                add_data = ., # data.frame
                meta_key = "raw", # string
                igrnore_rownm = FALSE)
pathogen_genes <- c("NS5") # 节约时间
STID_obj <- SpotDetect_Gene(STID_obj,
                            features = pathogen_genes,
                            feature_colnm = grep("all_gene",colnames(STID_obj@meta.data),value = T),
                            PosThres_prob = 0, PosThres_count = 1,
                            # col = c("#3D3576","#92D74D"),
                            # col = c("#8DBDDC","#FD9F8F"),
                            col = COLOR_DIS_CON,
                            black_bg = F,pt_size = 0.25,
                            # blur_method = "isoblur",
                            blur_method = NULL,
                            blur_n = 1,blur_sigma = 0.5, 
                            plot_method = "single",
                            grp_nm = "JEV_correct_before_all_gene_white")
STID_obj <- SpotDetect_Gene(STID_obj,
                            # features = host_genes,
                            # features = c("Ccl5","Irf7"),
                            features = host_genes,
                            feature_colnm = grep("all_gene",colnames(STID_obj@meta.data),value = T),
                            PosThres_prob = 0, PosThres_count = 4,
                            # col = c("grey95","red"),
                            col = COLOR_DIS_CON,
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#512DA8"),
                            black_bg = F,pt_size = 0.25,
                            blur_method = NULL, blur_n = 1,blur_sigma = 0.5, 
                            plot_method = "single",
                            grp_nm = "JEV_correct_before_host_gene_white")

# geneset
geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/Mouse_PCD_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
geneset_df <- geneset_df[c(2,5)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,seed = 10,
                               PosThres_prob = 0.75, PosThres_score = 0, # 最终确定是0.75
                               pt_size = 0.25,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "JEV_correct_before_PCD_white")

#>
geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/GO/Mouse_GO_BP_Detect_viral_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
colnames(geneset_df) <- gsub("GOBP_","",colnames(geneset_df))
geneset_df <- geneset_df[c(2,11)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                               PosThres_prob = 0.75, PosThres_score = 0, # 最终改为0.75
                               pt_size = 0.25,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "JEV_correct_before_GO_viral_white")
# saveRDS(STID_obj,file = "./rds/STID_obj_JEV_SpotDetect_D5_1.rds")
saveRDS(STID_obj,file = "./rds/STID_obj_JEV_SpotDetect_D3_1.rds")


## 2.4 Tbb ####
# gene
STID_obj <- GetGeneStat(STID_obj = STID_obj, features = host_genes,prefix = "all_gene",func = "sum") %>% 
  AddMetaColumn(STID_obj = STID_obj,
                add_data = ., # data.frame
                meta_key = "raw", # string
                igrnore_rownm = FALSE)
STID_obj <- SpotDetect_Gene(STID_obj,
                            features = pathogen_genes,
                            feature_colnm = grep("all_gene",colnames(STID_obj@meta.data),value = T),
                            PosThres_prob = 0, PosThres_count = 2,
                            # col = c("#3D3576","#92D74D"),
                            # col = c("#8DBDDC","#FD9F8F"),
                            col = COLOR_DIS_CON,
                            black_bg = F,pt_size = 2.2,vmax = "p99",
                            # blur_method = "isoblur",
                            blur_method = NULL,
                            blur_n = 1.5,blur_sigma = 0.5, 
                            plot_method = "single",
                            grp_nm = "Tbb_correct_before_all_gene_white")
STID_obj <- SpotDetect_Gene(STID_obj,
                            features = host_genes,
                            feature_colnm = grep("all_gene",colnames(STID_obj@meta.data),value = T),
                            PosThres_prob = 0, PosThres_count = 2,
                            # col = c("grey95","red"),
                            col = COLOR_DIS_CON,
                            black_bg = F,pt_size = 2.2,
                            blur_method = NULL, blur_n = 1,blur_sigma = 0.5, 
                            plot_method = "single",
                            grp_nm = "Tbb_correct_before_host_gene_white")


# geneset，后续选择连续值，不受设置的阈值的影响
geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/Mouse_PCD_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
geneset_list <- geneset_list[5]
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,seed = 10,
                               PosThres_prob = 0, PosThres_score = 7.5,
                               pt_size = 2.2,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "Tbb_correct_before_PCD_white")


geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/KEGG/Mouse_KEGG_Detect_parasitic_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
colnames(geneset_df) <- gsub("GOBP_","",colnames(geneset_df))
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
geneset_list <- geneset_list[2]
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                               PosThres_prob = 0, PosThres_score = 1.75, # 最终固定为1.75
                               pt_size = 2.2,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "Tbb_correct_before_KEGG_parasite_white")


#### 二、Supp Figure1 ####
### 1.A-D：TB ####
### 1.1 CorrectBackgroud
STID_obj_before <- STID_obj
STID_obj_after <- CorrectBackgroud(STID_obj = STID_obj_before, 
                                   bg_samp_id = c("WT_1"), 
                                   # bg_samp_id = c("PI1d_1"), 
                                   bg_features = pathogen_genes, 
                                   PosThres_prob = 0.9, adjust_UMI = F,
                                   assay_id = "Spatial", layer_id = "counts",
                                   grp_nm = "Final_TB_0.95_noadj_20260425", dir_nm = "M1_CorrectBackgroud")

### 1.2 SpotDetect_Gene
# gene
high_exp_genes <- GetTopGenes(STID_obj_before, top_n = 2, pattern = "^Rv",
                              grp_by_samp = F, grp_by_celltype = F,
                              assay_id = "Spatial", layer_id = "counts")
STID_obj_before <- GetGeneStat(STID_obj = STID_obj_before, features = high_exp_genes,prefix = "top10_gene",func = "sum") %>% 
  AddMetaColumn(STID_obj = STID_obj_before,
                add_data = ., # data.frame
                meta_key = "raw", # string
                igrnore_rownm = FALSE)
STID_obj_before <- GetGeneStat(STID_obj = STID_obj_before, features = pathogen_genes,prefix = "all_gene",func = "sum") %>% 
  AddMetaColumn(STID_obj = STID_obj_before,
                add_data = ., # data.frame
                meta_key = "raw",
                igrnore_rownm = FALSE)
STID_obj_after <- GetGeneStat(STID_obj = STID_obj_after, features = high_exp_genes,prefix = "top10_gene",func = "sum") %>%
  AddMetaColumn(STID_obj = STID_obj_after,
                add_data = ., # data.frame
                meta_key = "raw",
                igrnore_rownm = FALSE)
STID_obj_after <- GetGeneStat(STID_obj = STID_obj_after, features = pathogen_genes,prefix = "all_gene",func = "sum") %>%
  AddMetaColumn(STID_obj = STID_obj_after,
                add_data = ., # data.frame
                meta_key = "raw",
                igrnore_rownm = FALSE)
meta_data <- STID_obj_before@meta.data
STID_obj_before <- SpotDetect_Gene(STID_obj_before,
                                   features = high_exp_genes,
                                   feature_colnm = grep("all_gene",colnames(STID_obj_before@meta.data),value = T),
                                   PosThres_prob = 0, PosThres_count = 0,
                                   # col = c("#3D3576","#92D74D"),
                                   # col = c("#8DBDDC","#FD9F8F"),
                                   col = c("grey20","#80FFFF"),
                                   black_bg = T,pt_size = 0.75, # ctrl是0.75
                                   # blur_method = "isoblur",
                                   blur_method = NULL,
                                   blur_n = 1,blur_sigma = 0.5, 
                                   plot_method = "single",
                                   grp_nm = "TB_correct_before_all_gene_black_20260425")
STID_obj_after <- SpotDetect_Gene(STID_obj_after,
                                  features = high_exp_genes,
                                  feature_colnm = grep("all_gene",colnames(STID_obj_after@meta.data),value = T),
                                  PosThres_prob = 0, PosThres_count = 0,
                                  # col = c("#3D3576","#92D74D"),
                                  # col = c("#8DBDDC","#FD9F8F"),
                                  col = c("grey20","#80FFFF"),
                                  black_bg = T,pt_size = 0.75, # ctrl是0.75
                                  # blur_method = "isoblur",
                                  blur_method = NULL,
                                  blur_n = 1,blur_sigma = 0.5, 
                                  plot_method = "single",
                                  grp_nm = "TB_correct_after_all_gene_black_20260425")
STID_obj_after


### 2.Supp Figure1: D ####

## 2.1 Lang ####
# gene
STID_obj <- GetGeneStat(STID_obj = STID_obj, features = pathogen_genes,prefix = "all_gene",func = "sum") %>% 
  AddMetaColumn(STID_obj = STID_obj,
                add_data = ., # data.frame
                meta_key = "raw", # string
                igrnore_rownm = FALSE)
if(0){
  SpotDetect_Gene(STID_obj,
                  features = pathogen_genes,
                  feature_colnm = grep("all_gene",colnames(STID_obj@meta.data),value = T), # 没有任何基因
                  PosThres_prob = 0, PosThres_count = 0,
                  # col = c("#3D3576","#92D74D"),
                  # col = c("#8DBDDC","#FD9F8F"),
                  col = COLOR_DIS_CON,
                  black_bg = F,pt_size = 2.4,
                  # blur_method = "isoblur",
                  blur_method = NULL,
                  blur_n = 1,blur_sigma = 0.5, 
                  plot_method = "single",
                  grp_nm = "Lang_correct_before_all_gene_white")
}
options("parallel_workers" = 1)
STID_obj <- SpotDetect_Gene(STID_obj,
                            features = host_genes,
                            feature_colnm = grep("all_gene",colnames(STID_obj@meta.data),value = T),
                            PosThres_prob = 0, PosThres_count = 4,
                            # col = c("grey95","red"),
                            col = COLOR_DIS_CON,
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#512DA8"),
                            black_bg = F,pt_size = 2.4,
                            blur_method = NULL, blur_n = 1,blur_sigma = 0.5, 
                            plot_method = "single",
                            grp_nm = "Lang_correct_before_host_gene_white")
meta_data <- GetMetaData(STID_obj,meta_key = "M1_SpotDetect_Gene_Lang_correct_before_host_gene_white")[[1]]

# geneset
geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/Mouse_PCD_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
geneset_df <- geneset_df[c(1,5)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,seed = 10,
                               PosThres_prob = 0.8, PosThres_score = 0, # 0.8
                               pt_size = 2,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "Lang_correct_before_PCD_white")

geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/GO/Mouse_GO_BP_Detect_viral_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
colnames(geneset_df) <- gsub("GOBP_","",colnames(geneset_df))
# geneset_df <- geneset_df[c(2,11)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
geneset_list <- geneset_list[11]
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                               PosThres_prob = 0.8, PosThres_score = 0, # 0.8
                               pt_size = 2,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "Lang_correct_before_GO_viral_white")
meta_data <- GetMetaData(STID_obj,meta_key = "M1_SpotDetect_Geneset_Lang_correct_before_GO_viral_white")[[1]]

if(0){
  # 平滑
  geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/GO/Mouse_GO_BP_Detect_viral_geneset.txt",
                           sep = "\t",row.names = NULL,header = T,na.strings = "")
  colnames(geneset_df) <- gsub("GOBP_","",colnames(geneset_df))
  geneset_df <- geneset_df[c(2,11)]
  geneset_list <- lapply(geneset_df, function(x) na.omit(x))
  names(geneset_list)
  SpotDetect_Geneset(STID_obj,
                     geneset_list = geneset_list,
                     score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                     PosThres_prob = 0.8, PosThres_score = 0,
                     pt_size = 2,
                     col = COLOR_DIS_CON,
                     black_bg = F, blur_method = "isoblur",sigma = 0.2,
                     plot_method = "single",
                     grp_nm = "Lang_correct_before_GO_viral_isoblur")
}


## 2.2 VVE ####
# gene
STID_obj <- SpotDetect_Gene(STID_obj,
                            features = host_genes,
                            feature_colnm = NULL,
                            PosThres_prob = 0, PosThres_count = 4,
                            # col = c("grey95","red"),
                            col = COLOR_DIS_CON,
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#512DA8"),
                            black_bg = F,pt_size = 2.2,
                            blur_method = NULL, blur_n = 1,blur_sigma = 0.5, 
                            plot_method = "single",
                            grp_nm = "VVE_correct_before_host_gene_white")

# geneset
geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/Mouse_PCD_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
# geneset_df <- geneset_df[c(2,5)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,seed = 10,
                               PosThres_prob = 0.8, PosThres_score = 0,
                               pt_size = 2,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "VVE_correct_before_PCD_white")

geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/GO/Mouse_GO_BP_Detect_viral_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
colnames(geneset_df) <- gsub("GOBP_","",colnames(geneset_df))
# geneset_df <- geneset_df[c(2,11)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                               PosThres_prob = 0.8, PosThres_score = 0,
                               pt_size = 2,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "VVE_correct_before_GO_viral_white")

## 2.3 Kp ####
# gene
STID_obj <- SpotDetect_Gene(STID_obj,
                            features = host_genes,
                            feature_colnm = NULL,
                            PosThres_prob = 0, PosThres_count = 4,
                            # col = c("grey95","red"),
                            # col = viridis(5) ,
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#512DA8"),
                            vmax = "p95",
                            black_bg = F,pt_size = 2.2,
                            blur_method = NULL, blur_n = 1,blur_sigma = 0.5, 
                            plot_method = "single",
                            grp_nm = "Kp_correct_before_host_gene_white")

# geneset
geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/Mouse_PCD_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
geneset_df <- geneset_df[c(1,5)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,seed = 10,
                               PosThres_prob = 0.8, PosThres_score = 0,
                               pt_size = 2, # 不要用2.2，就是2.0最合适
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "Kp_correct_before_PCD_white")


geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/GO/Mouse_GO_BP_Detect_bacterial_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
colnames(geneset_df) <- gsub("GOBP_","",colnames(geneset_df))
# geneset_df <- geneset_df[c(2,11)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
geneset_list <- geneset_list[1]
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                               PosThres_prob = 0.85, PosThres_score = 0,
                               pt_size = 2, # 不要用2.2，就是2.0最合适
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "Kp_correct_before_GO_bacterial_white")


## 2.4 Malaria ####
# gene
STID_obj <- GetGeneStat(STID_obj = STID_obj, pattern = NULL, features = pathogen_genes,
                        prefix = "all_gene",func = "sum") %>% 
  AddMetaColumn(STID_obj = STID_obj,
                add_data = ., # data.frame
                meta_key = "raw", # string
                igrnore_rownm = FALSE)
STID_obj <- GetGeneStat(STID_obj = STID_obj, pattern = NULL, features = host_genes,
                        prefix = "ISG",func = "sum") %>% 
  AddMetaColumn(STID_obj = STID_obj,
                add_data = ., # data.frame
                meta_key = "raw", # string
                igrnore_rownm = FALSE)
STID_obj <- SpotDetect_Gene(STID_obj,
                            features = host_genes,
                            feature_colnm = c(grep("all_gene",colnames(STID_obj@meta.data),value = T),
                                              grep("ISG",colnames(STID_obj@meta.data),value = T)),
                            PosThres_prob = 0, PosThres_count = 30,
                            # col = c("grey95","red"),
                            col = COLOR_DIS_CON,
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#512DA8"),
                            black_bg = F,pt_size = 2,
                            blur_method = NULL, blur_n = 1,blur_sigma = 0.5, 
                            plot_method = "single",
                            grp_nm = "Ma_correct_before_host_gene_white")


# geneset，后续选择连续值，不受设置的阈值的影响
geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/Mouse_PCD_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
# geneset_df <- geneset_df[c(2,5)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,seed = 10,
                               PosThres_prob = 0.75, PosThres_score = 0,
                               pt_size = 2,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "Ma_correct_before_PCD_white")


geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/KEGG/Mouse_KEGG_Detect_parasitic_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
colnames(geneset_df) <- gsub("GOBP_","",colnames(geneset_df))
# geneset_df <- geneset_df[c(2,11)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                               PosThres_prob = 0.75, PosThres_score = 0,
                               pt_size = 2,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "Ma_correct_before_KEGG_parasite_white")


#### 三、Figure3 ####
### 1.BCDEFGHIJK ####
## 1.1 TB：Microbe ####
STID_obj_after %>% print()
meta_data <- GetMetaData(
  STID_obj_after,
  meta_key = "M1_SpotDetect_Gene_TB_correct_after_all_gene_white"
)[[1]]
STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj_after, meta_key = "M1_SpotDetect_Gene_TB_correct_after_all_gene_white", 
                                   loop_id = "PI8W_1",
                                   spatial_scale_method = "region", region_detect_method = "convex", update_spots = F,
                                   ROI_size = NULL, k_kNNdist = 4,
                                   pos_colnm = "Label_all_gene_nFeature(sum)", 
                                   description = NULL,grp_nm = "TB_DBSCAN", dir_nm = "M2_NicheDetect_STS")
STID_obj_detect %>% print()
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_TB_DBSCAN",
                           add_coord = F)[[1]]

# > plot
SEVEN_DARK <- c("#AFBF41FF" ,"#50C49FFF"  ,"#FCB11C" ,"#F81B02FF" ,"#FC7715FF","#3B95C4FF" ,"#B560D4FF") # 7色，深色1，常用
SEVEN_LIGHT <- c("#BBBFA1", "#B2C4BE", "#FCDB9A","#F88A7E", "#FCC093",  "#9DB7C4", "#D1CAD4") # 7色，淡色3，折线图
pdf("./photo/Z_other/20260408_Figure3_result/TB/TB_microbe_edge_raw.pdf",width = 30,height = 10)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "sample", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/TB/TB_microbe_region_raw.pdf",width = 30,height = 10)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "All_ROI_label2",
             facet_grpnm = "sample", datatype = "discrete",
             col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()


# > Plot_DistLine_Exp
STID_obj_detect %>% print()
tmp <- GetMetaData(STID_obj_detect,
                   meta_key = list(c("M1_SpotDetect_Geneset_TB_correct_after_PCD_white")))[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/TB/TB_microbe_DistLine_Exp_host.pdf",width = 4.8,height = 6)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = c("Rvnr01","Igkc"), feature_colnm = "all_gene_nFeature(sum)", 
                  loop_id = "PI8W_1", col = c("#F81B02FF"  ,"#3B95C4FF","#F81B02FF" ) ,
                  meta_key = list(c("M1_SpotDetect_Gene_TB_correct_after_all_gene_white",
                                    "M2_NicheDetect_STS_TB_DBSCAN")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/TB/TB_microbe_DistLine_Exp_geneset.pdf",width = 4.8,height = 3)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Tuberculosis"),
                  loop_id = "PI8W_1", col = "#3B95C4FF",exp_scale = F,
                  meta_key = list(c("M1_SpotDetect_Geneset_TB_correct_after_KEGG_bacterial_white",
                                    "M2_NicheDetect_STS_TB_DBSCAN")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/TB/TB_microbe_DistLine_Exp_PCD.pdf",width = 5.8, height = 3)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Lysosome.dependentcelldeath"), 
                  loop_id = "PI8W_1", col = "#3B95C4FF",
                  meta_key = list(c("M1_SpotDetect_Geneset_TB_correct_after_PCD_white",
                                    "M2_NicheDetect_STS_TB_DBSCAN")))
dev.off()

## 1.2 TB：host ####
STID_obj_detect %>% print()
meta_data <- GetMetaData(
  STID_obj_detect,
  meta_key = "M1_SpotDetect_Geneset_TB_correct_after_KEGG_bacterial_white"
)[[1]]
colnames(meta_data)
STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj_detect, meta_key = "M1_SpotDetect_Geneset_TB_correct_after_KEGG_bacterial_white", 
                                   loop_id = "PI8W_1",
                                   spatial_scale_method = "region", region_detect_method = "convex", update_spots = F,
                                   density_thres = 0.4,
                                   ROI_size = 200,
                                   pos_colnm = "Label_Tuberculosis", 
                                   description = NULL,grp_nm = "TB_host_DBSCAN", dir_nm = "M2_NicheDetect_STS")
STID_obj_detect %>% print()
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_TB_host_DBSCAN",
                           add_coord = F)[[1]]


# > plot
# SEVEN_DARK <- c("#B560D4FF" ,"#F81B02FF","#FCB11C" ,"#FC7715FF") # 7色，深色1，常用
# SEVEN_LIGHT <- c("#D1CAD4","#F88A7E","#FCDB9A",  "#FCC093") # 7色，淡色3，折线图
SEVEN_DARK <- c("#B560D4FF" ,"#FC7715FF" ,"#FCB11C" ,"#AFBF41FF" ,"#50C49FFF" ,"#3B95C4FF" ,"#F81B02FF")
SEVEN_LIGHT <- c("#D1CAD4", "#FCC093", "#FCDB9A", "#BBBFA1", "#B2C4BE", "#9DB7C4", "#F88A7E")
pdf("./photo/Z_other/20260408_Figure3_result/TB/TB_host_edge_raw.pdf",width = 30,height = 10)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "sample", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/TB/TB_host_region_raw.pdf",width = 30,height = 10)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "All_ROI_label2",
             facet_grpnm = "sample", datatype = "discrete",
             col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

# > Plot_DistLine_Exp
STID_obj_detect %>% print()
tmp <- GetMetaData(STID_obj_detect,
                   meta_key = list(c("M1_SpotDetect_Gene_TB_correct_after_all_gene_white")))[[1]]
colnames(tmp)
pdf("./photo/Z_other/20260408_Figure3_result/TB/TB_host_DistLine_Exp_host.pdf",width = 4.8,height = 6)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = c("Rvnr01","Igkc"), feature_colnm = "all_gene_nFeature(sum)", 
                  loop_id = "PI8W_1", col = c("#F81B02FF"  ,"#3B95C4FF","#F81B02FF" ) ,
                  meta_key = list(c("M1_SpotDetect_Gene_TB_correct_after_all_gene_white",
                                    "M2_NicheDetect_STS_TB_host_DBSCAN")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/TB/TB_host_DistLine_Exp_geneset.pdf",width = 4.8,height = 3)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Tuberculosis"),
                  loop_id = "PI8W_1", col = "#3B95C4FF",exp_scale = F,
                  meta_key = list(c("M1_SpotDetect_Geneset_TB_correct_after_KEGG_bacterial_white",
                                    "M2_NicheDetect_STS_TB_host_DBSCAN")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/TB/TB_host_DistLine_Exp_PCD.pdf",width = 5.8, height = 3)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Lysosome.dependentcelldeath"), 
                  loop_id = "PI8W_1", col = "#3B95C4FF",
                  meta_key = list(c("M1_SpotDetect_Geneset_TB_correct_after_PCD_white",
                                    "M2_NicheDetect_STS_TB_host_DBSCAN")))
dev.off()

#> CompareNiche
STID_obj_detect
pdf("./photo/Z_other/20260408_Figure3_result/TB/TB_CompareNiche_microbe_host.pdf",width = 5,height = 8)
CompareNiche(STID_obj = STID_obj_detect, 
             meta_key1 = "M2_NicheDetect_STS_TB_DBSCAN",
             meta_key2 = "M2_NicheDetect_STS_TB_host_DBSCAN")
dev.off()
saveRDS(STID_obj_detect,"./outputdata/STID_obj_detect_TB.rds")


## 1.3 JEV：microbe ####
# 先运行## 2.3 JEV 
STID_obj
meta_data <- GetMetaData(
  STID_obj,
  meta_key = "M1_SpotDetect_Gene_JEV_correct_before_all_gene_white"
)[[1]]


#>TSO_ROI
# 还没开发好
if(0){
  STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj, meta_key = "M1_SpotDetect_Gene_JEV_correct_before_all_gene_white", 
                                     spatial_scale_method = "spot", update_spots = T,
                                     ROI_size = NULL,
                                     pos_colnm = "Label_all_gene_nFeature(sum)", 
                                     description = NULL,grp_nm = "STS_JEV_microbe_spot", dir_nm = "M2_NicheDetect_STS")
}


#>DBSCAN
STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj, meta_key = "M1_SpotDetect_Gene_JEV_correct_before_all_gene_white", 
                                   spatial_scale_method = "region", region_detect_method = "convex", update_spots = F,
                                   ROI_size = NULL,, density_thres = 1, # 原始是1，后来改为0.9
                                   pos_colnm = "Label_all_gene_nFeature(sum)", 
                                   description = NULL,grp_nm = "STS_JEV_microbe_region", dir_nm = "M2_NicheDetect_STS")
STID_obj_detect
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_STS_JEV_microbe_region",
                           add_coord = F)[[1]]

#> plot
# SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C"  ,"#50C49FFF" ,"#3B95C4FF" ,"#B560D4FF")
# SEVEN_LIGHT <- c("#F88A7E", "#FCC093", "#FCDB9A", "#BBBFA1", "#9DB7C4", "#D1CAD4")
SEVEN_DARK <- c("#50C49FFF" ,"#FC7715FF" ,"#FCB11C"  ,"#F81B02FF" ,"#3B95C4FF" ,"#B560D4FF")
SEVEN_LIGHT <- c("#BBBFA1",  "#FCC093", "#FCDB9A","#F88A7E", "#9DB7C4", "#D1CAD4")

pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_microbe_edge_raw.pdf",width = 15,height = 15)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "new_samp",
             datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_microbe_region_raw.pdf",width = 15,height = 15)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "All_ROI_label2",
             facet_grpnm = "new_samp",
             datatype = "discrete",
             col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

#> Plot_DistLine_Exp
STID_obj_detect %>% print()
tmp <- GetMetaData(STID_obj_detect,
                   meta_key = list(c("M2_NicheDetect_STS_STS_JEV_microbe_region")))[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_microbe_DistLine_Exp_host.pdf",width = 4.8,height = 6)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = c("NS5","Ccl2"), feature_colnm = "all_gene_nFeature(sum)", 
                  loop_id = "D5_1", col = c("#F81B02FF"  ,"#3B95C4FF","#F81B02FF" ) ,
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Gene_JEV_correct_before_all_gene_white",
                                                        "M2_NicheDetect_STS_STS_JEV_microbe_region")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_microbe_DistLine_Exp_geneset.pdf",width = 5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("RESPONSE_TO_VIRUS"),
                  loop_id = "D5_1", col = "#3B95C4FF",
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_GO_viral_white",
                                                        "M2_NicheDetect_STS_STS_JEV_microbe_region")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_microbe_DistLine_Exp_PCD.pdf",width = 4.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Necroptosis"), 
                  loop_id = "D5_1", col = "#3B95C4FF",
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_PCD_white",
                                                        "M2_NicheDetect_STS_STS_JEV_microbe_region")))
dev.off()


# > 根据病原最低点，扩充距离
STID_obj_detect %>% print()
tmp <- GetMetaData(STID_obj_detect,
                   meta_key = "M2_NicheDetect_STS_STS_JEV_microbe_region")[[1]]
STID_obj_expand <- NicheExpand(STID_obj_detect, meta_key = "M2_NicheDetect_STS_STS_JEV_microbe_region",
                               pos_colnm = "ROI_label",center_colnm = "ROI_center",expand_dist =30) 
STID_obj_expand %>% print()
meta_key <- "M2_NicheExpand_20260408_231317"
expand_meta <- GetMetaData(STID_obj_expand, meta_key = meta_key,add_coord = T)[[1]]

#> plot
# SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C"   ,"#3B95C4FF" ,"#B560D4FF")
# SEVEN_LIGHT <- c("#F88A7E", "#FCC093", "#FCDB9A", "#9DB7C4", "#D1CAD4")
SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C"  ,"#50C49FFF" ,"#3B95C4FF" ,"#B560D4FF")
SEVEN_LIGHT <- c("#F88A7E", "#FCC093", "#FCDB9A", "#BBBFA1", "#9DB7C4", "#D1CAD4")
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_microbe_edge_expand.pdf",width = 15,height = 15)
Plot_Spatial(plot_data = expand_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "new_samp", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_microbe_region_expand.pdf",width = 15,height = 15)
Plot_Spatial(plot_data = expand_meta,x_colnm = "x",y_colnm = "y",group_by = "All_ROI_label2",
             facet_grpnm = "new_samp", datatype = "discrete",
             col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

#> Plot_DistLine_Exp
STID_obj_expand %>% print()
tmp <- GetMetaData(STID_obj_expand,
                   meta_key = list(c(meta_key)))[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_microbe_DistLine_Exp_host_expand.pdf",width = 4.8,height = 6)
Plot_DistLine_Exp(STID_obj = STID_obj_expand, features = c("NS5","Ccl2"), feature_colnm = "all_gene_nFeature(sum)", 
                  loop_id = "D5_1", col = c("#F81B02FF"  ,"#3B95C4FF","#F81B02FF" ) ,
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Gene_JEV_correct_before_all_gene_white",
                                                        meta_key)))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_microbe_DistLine_Exp_geneset_expand.pdf",width = 5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_expand, features = NULL, feature_colnm = c("RESPONSE_TO_VIRUS"),
                  loop_id = "D5_1", col = "#3B95C4FF",
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_GO_viral_white",
                                                        meta_key)))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_microbe_DistLine_Exp_PCD_expand.pdf",width = 4.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_expand, features = NULL, feature_colnm = c("Necroptosis"), 
                  loop_id = "D5_1", col = "#3B95C4FF",
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_PCD_white",
                                                        meta_key)))
dev.off()


## 1.4 JEV：host ####
STID_obj_detect
meta_data <- GetMetaData(
  STID_obj_detect,
  meta_key = "M1_SpotDetect_Geneset_JEV_correct_before_GO_viral_white"
)[[1]]
STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj_detect, meta_key = "M1_SpotDetect_Geneset_JEV_correct_before_GO_viral_white", 
                                   spatial_scale_method = "region", region_detect_method = "convex", update_spots = T,
                                   pos_colnm = "Label_RESPONSE_TO_VIRUS",
                                   ROI_size = NULL,
                                   density_thres = 0.3,
                                   description = NULL,grp_nm = "STS_JEV_host_region", dir_nm = "M2_NicheDetect_STS")
STID_obj_detect
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_STS_JEV_host_region",add_coord = F)[[1]]

#> plot
SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C"  ,"#B560D4FF")
SEVEN_LIGHT <- c("#F88A7E", "#FCC093", "#FCDB9A", "#D1CAD4")
# SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C"  ,"#50C49FFF" ,"#3B95C4FF" ,"#B560D4FF")
# SEVEN_LIGHT <- c("#F88A7E", "#FCC093", "#FCDB9A", "#BBBFA1", "#9DB7C4", "#D1CAD4")
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_host_edge_raw.pdf",width = 15,height = 15)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "new_samp", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_host_region_raw.pdf",width = 15,height = 15)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "All_ROI_label2",
             facet_grpnm = "new_samp", datatype = "discrete",
             col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

#> Plot_DistLine_Exp
STID_obj_detect %>% print()
tmp <- GetMetaData(STID_obj_detect,
                   meta_key = list(c("M2_NicheDetect_STS_STS_JEV_host_region")))[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_host_DistLine_Exp_host.pdf",width = 4.8,height = 6)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = c("NS5","Ccl2"), feature_colnm = "all_gene_nFeature(sum)", 
                  loop_id = "D5_1", col = c("#F81B02FF"  ,"#3B95C4FF","#F81B02FF" ) ,
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Gene_JEV_correct_before_all_gene_white",
                                                        "M2_NicheDetect_STS_STS_JEV_host_region")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_host_DistLine_Exp_geneset.pdf",width = 5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("RESPONSE_TO_VIRUS"),
                  loop_id = "D5_1", col = "#3B95C4FF",
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_GO_viral_white",
                                                        "M2_NicheDetect_STS_STS_JEV_host_region")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_host_DistLine_Exp_PCD.pdf",width = 4.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Necroptosis"), 
                  loop_id = "D5_1", col = "#3B95C4FF",
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_PCD_white",
                                                        "M2_NicheDetect_STS_STS_JEV_host_region")))
dev.off()


### 1.5 JEV：CompareNiche ####
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_CompareNiche_microbe_host.pdf",width = 5,height = 8)
CompareNiche(STID_obj = STID_obj_detect, 
             meta_key1 = "M2_NicheDetect_STS_STS_JEV_microbe_region",
             meta_key2 = "M2_NicheDetect_STS_STS_JEV_host_region",
             bins = 15)
dev.off()

#> DBSCAN concave
STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj_detect, meta_key = "M1_SpotDetect_Gene_JEV_correct_before_all_gene_white", 
                                   spatial_scale_method = "region", region_detect_method = "concave", update_spots = F,
                                   ROI_size = NULL, density_thres = 1, concavity = 2.5,
                                   pos_colnm = "Label_all_gene_nFeature(sum)", 
                                   description = NULL,grp_nm = "STS_JEV_microbe_concave", dir_nm = "M2_NicheDetect_STS")

#> 
expand_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_STS_JEV_microbe_concave",add_coord = T)[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_microbe_concave_edge_raw.pdf",width = 15,height = 15)
Plot_Spatial(plot_data = expand_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "grp", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

STID_obj_detect %>% print()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_CompareNiche_microbe_region_concave.pdf",width = 5,height = 8)
CompareNiche(STID_obj = STID_obj_detect, 
             meta_key1 = "M2_NicheDetect_STS_STS_JEV_microbe_region",
             meta_key2 = "M2_NicheDetect_STS_STS_JEV_microbe_concave")
dev.off()
saveRDS(STID_obj_detect,file = "./rds/STID_obj_detect_JEV.rds")


## 1.6 CE: D4_2, D91_1 ####
if(0){
  table(STID_obj_expand@meta.data$anno)
  celltypes <- STID_obj_expand@meta.data$anno %>% levels()
  col_lasso <- COLOR_LIST$PALETTE_WHITE_BG[1:10]
  names(col_lasso) <- celltypes
  col_lasso2 <- col_lasso
  col_lasso["Hepatocytes"] <- "grey95"
}else{
  col_lasso <- c(
    "HsPCs" = "#E41A1C",
    "Hepatocytes" = "grey95",
    "Infla Heps" = "#4DAF4A",
    "Fibroblasts" = "#984EA3",
    "Cho/Spp1+ cells" = "#FFFF33",
    "Spp1+ MoMFs" = "#FF7F00",
    "MoKCs" = "#377EB8",
    "Neutrophils" = "#F781BF",
    "B/plasma cells" = "#A65628",
    "Others" = "#8DA0CB"
  )
  col_lasso2 <- c(
    "HsPCs" = "#E41A1C",
    "Hepatocytes" = "#66C2A5",
    "Infla Heps" = "#4DAF4A",
    "Fibroblasts" = "#984EA3",
    "Cho/Spp1+ cells" = "#FFFF33",
    "Spp1+ MoMFs" = "#FF7F00",
    "MoKCs" = "#377EB8",
    "Neutrophils" = "#F781BF",
    "B/plasma cells" = "#A65628",
    "Others" = "#8DA0CB"
  )
}

STID_obj_lasso <- NicheDetect_Lasso(STID_obj_SS_CE, meta_key = "coord",group_by = "anno",
                                    col = col_lasso,grp_nm = "test")
STID_obj_lasso %>% print()
meta_key <- "M2_NicheDetect_Lasso_test"
lasso_meta <- GetMetaData(STID_obj_lasso, meta_key = meta_key, # 随时更换
                          add_coord = F)[[1]]
saveRDS(STID_obj_lasso,file = "./rds/STID_obj_CE_lasso.rds")


#> 细胞注释展示图
pdf("./photo/Z_other/20260408_Figure3_result/CE_Lasso_Spatial.pdf",width = 20,height = 5)
Plot_Spatial(STID_obj = STID_obj_lasso,
             facet_grpnm = "batch", datatype = "discrete",
             col = list(dis = col_lasso,con = NULL),
             pt_size = 1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

# > plot
SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF","#FCB11C") # 7色，深色1，常用
SEVEN_LIGHT <- c("#F88A7E", "#FCC093","#FCDB9A") # 7色，淡色3，折线图
pdf("./photo/Z_other/20260408_Figure3_result/CE_Lasso_edge_raw.pdf",width = 20,height = 5)
Plot_Spatial(plot_data = lasso_meta,x_colnm = "x", y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "batch", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/CE_Lasso_region_raw.pdf",width = 20,height = 5)
Plot_Spatial(plot_data = lasso_meta,x_colnm = "x", y_colnm = "y",group_by = "All_ROI_label2",
             facet_grpnm = "batch", datatype = "discrete",
             col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()


# > 宿主反应和病原同时展示图
host_meta <- GetMetaData(STID_obj_lasso,meta_key = list(c("M1_SpotDetect_Gene_CE_correct_after_host_gene_white")))[[1]]
host_meta <- host_meta %>% 
  mutate(grp = if_else(Label_Il1b == "pos",
                       if_else(`Label_all_gene_nFeature(sum)` == "pos","co","Host"),
                       if_else(`Label_all_gene_nFeature(sum)` == "pos","Pathogen","none"))) %>% 
  mutate(grp = factor(grp,levels = c("none","Pathogen","Host","co")))
table(host_meta$grp)
pdf("./photo/Z_other/20260408_Figure3_result/CE_host_pathogen_Spatial.pdf",width = 20,height = 5)
Plot_Spatial(plot_data = host_meta,x_colnm = "x", y_colnm = "y",group_by ="grp",
             facet_grpnm = "batch", datatype = "discrete",
             col = list(dis = c("grey95","#FCB11C","#F81B02FF","#5BBCD6"),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

#>
tmp1 <- GetMetaData(STID_obj_lasso,meta_key = list(c("M1_SpotDetect_Gene_CE_correct_after_host_gene_white",meta_key)))[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/CE_DistLine_Exp_host.pdf",width = 8,height = 9)
Plot_DistLine_Exp(STID_obj = STID_obj_lasso, features = c("EmuJ-002209100","Spp1","Il1b"), feature_colnm = "all_gene_nFeature(sum)",
                  # loop_id = "DPI_4_2", 
                  col = c("#F81B02FF","#3B95C4FF","#3B95C4FF","#F81B02FF") ,
                  facet_grpnm = "batch",meta_key = list(c("M1_SpotDetect_Gene_CE_correct_after_host_gene_white",meta_key)))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/CE_DistLine_Exp_PCD.pdf",width = 7.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_lasso, features = NULL, feature_colnm = c("Necroptosis"), 
                  # loop_id = "DPI_4_2", 
                  col = c("#3B95C4FF"),
                  facet_grpnm = "batch",meta_key = list(c("M1_SpotDetect_Geneset_CE_correct_after_PCD_white",meta_key)))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/CE_DistLine_Exp_geneset.pdf",width = 7.2,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_lasso, features = NULL, feature_colnm = c("Malaria"), 
                  # loop_id = "DPI_4_2", 
                  col = c("#3B95C4FF"),
                  facet_grpnm = "batch",meta_key = list(c("M1_SpotDetect_Geneset_CE_correct_after_KEGG_Parasite_white",meta_key)))
dev.off()

#>>> 根据宿主反应区域，扩充microbe距离
tmp <- GetMetaData(STID_obj_lasso,meta_key = meta_key)[[1]]
STID_obj_expand <- NicheExpand(STID_obj_lasso, meta_key = meta_key,
                               pos_colnm = "ROI_label",center_colnm = "ROI_center",expand_dist = 8,
                               grp_nm = "CE") # 8 比较合适
STID_obj_expand %>% print()
meta_key <- "M2_NicheExpand_CE"
expand_meta <- GetMetaData(STID_obj_expand, meta_key = meta_key,add_coord = T)[[1]]

SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF","#FCB11C") # 7色，深色1，常用
SEVEN_LIGHT <- c("#F88A7E", "#FCC093","#FCDB9A") # 7色，淡色3，折线图
pdf("./photo/Z_other/20260408_Figure3_result/CE_Lasso_edge_expand.pdf",width = 20,height = 5)
Plot_Spatial(plot_data = expand_meta,x_colnm = "x", y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "batch", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/CE_Lasso_region_expand.pdf",width = 20,height = 5)
Plot_Spatial(plot_data = expand_meta,x_colnm = "x", y_colnm = "y",group_by = "All_ROI_label2",
             facet_grpnm = "batch", datatype = "discrete",
             col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

#>
STID_obj_expand %>% print()
tmp <- GetMetaData(STID_obj_expand,meta_key = meta_key)[[1]]
pdf()
Plot_DistLine_Exp(STID_obj = STID_obj_expand, features = c("EmuJ-002209100","Spp1","Il1b"), feature_colnm = "all_gene_nFeature(sum)",
                  # loop_id = "DPI_4_2",
                  col = c("#F81B02FF","#3B95C4FF","#3B95C4FF","#F81B02FF") ,
                  facet_grpnm = "batch",meta_key = list(c("M1_SpotDetect_Gene_CE_correct_after_all_gene_white",meta_key)))
Plot_DistLine_Exp(STID_obj = STID_obj_expand, features = NULL, feature_colnm = c("Necroptosis"), 
                  # loop_id = "DPI_4_2", 
                  col = c("#3B95C4FF"),
                  facet_grpnm = "batch",meta_key = list(c("M1_SpotDetect_Geneset_CE_correct_after_PCD_white",meta_key)))
Plot_DistLine_Exp(STID_obj = STID_obj_expand, features = NULL, feature_colnm = c("Malaria"), 
                  # loop_id = "DPI_4_2", 
                  col = c("#3B95C4FF"),
                  facet_grpnm = "batch",meta_key = list(c("M1_SpotDetect_Geneset_CE_correct_after_KEGG_Parasite_white",meta_key)))

#>
STID_obj_expand %>% print()
tmp <- GetMetaData(STID_obj_expand,meta_key = meta_key)[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/CE_DistLine_Exp_host_expand.pdf",width = 8,height = 9)
Plot_DistLine_Exp(STID_obj = STID_obj_expand, features = c("EmuJ-002209100","Spp1","Il1b"), feature_colnm = "all_gene_nFeature(sum)",
                  # loop_id = "DPI_4_2", 
                  col = c("#F81B02FF","#3B95C4FF","#3B95C4FF","#F81B02FF") ,
                  facet_grpnm = "batch",meta_key = list(c("M1_SpotDetect_Gene_CE_correct_after_host_gene_white",meta_key)))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/CE_DistLine_Exp_PCD_expand.pdf",width = 7.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_expand, features = NULL, feature_colnm = c("Necroptosis"), 
                  # loop_id = "DPI_4_2", 
                  col = c("#3B95C4FF"),
                  facet_grpnm = "batch",meta_key = list(c("M1_SpotDetect_Geneset_CE_correct_after_PCD_white",meta_key)))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/CE_DistLine_Exp_geneset_expand.pdf",width = 7.2,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_expand, features = NULL, feature_colnm = c("Malaria"), 
                  # loop_id = "DPI_4_2", 
                  col = c("#3B95C4FF"),
                  facet_grpnm = "batch",meta_key = list(c("M1_SpotDetect_Geneset_CE_correct_after_KEGG_Parasite_white",meta_key)))
dev.off()
saveRDS(STID_obj_expand,file = "./rds/STID_obj_CE_lasso_expand.rds")



### 2.JEV D3_1 microbe ####
# 先gene和geneset运行2.3 JEV：D5_1/D3_1


#> NicheDetect_Spot
STID_obj %>% print()
STID_obj <- NicheDetect_Spot(STID_obj = STID_obj,pos_colnm = "Label_all_gene_nFeature(sum)", 
                             meta_key = "M1_SpotDetect_Gene_JEV_correct_before_all_gene_white", 
                             description = NULL,
                             grp_nm = "D3_1")

# Plot_DistLine_Exp
STID_obj %>% print()
meta_key <- "M2_NicheDetect_Spot_D3_1"
tmp <- GetMetaData(STID_obj,meta_key = list(c(meta_key)))[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_D3_microbe_DistLine_Exp_host_raw.pdf",width = 4.8,height = 6)
Plot_DistLine_Exp(STID_obj = STID_obj, features = c("NS5","Ccl2"), feature_colnm = "all_gene_nFeature(sum)", 
                  loop_id = "D3_1", col = c("#F81B02FF"  ,"#3B95C4FF","#F81B02FF" ) ,
                  distance_scale = F,exp_scale = F,
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Gene_JEV_correct_before_all_gene_white",
                                                        meta_key)))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_D3_microbe_DistLine_Exp_geneset_raw.pdf",width = 5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj, features = NULL, feature_colnm = c("RESPONSE_TO_VIRUS"),
                  loop_id = "D3_1", col = "#3B95C4FF",
                  distance_scale = F,exp_scale = F,
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_GO_viral_white",
                                                        meta_key)))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/JEV/JEV_D3_microbe_DistLine_Exp_PCD_raw.pdf",width = 4.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj, features = NULL, feature_colnm = c("Necroptosis"), 
                  loop_id = "D3_1", col = "#3B95C4FF",
                  distance_scale = F,exp_scale = F,
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_PCD_white",
                                                        meta_key)))
dev.off()
# 目前不需要放拓展的



#### 四、supp Figure2 ####
### 1.Supp Figure2-FH ####
## 1.1 Tbb：microbe ####
STID_obj
meta_data <- GetMetaData(
  STID_obj,
  meta_key = "M1_SpotDetect_Gene_Tbb_correct_before_all_gene_white"
)[[1]]
colnames(meta_data)

STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj, meta_key = "M1_SpotDetect_Gene_Tbb_correct_before_all_gene_white", 
                                   spatial_scale_method = "region", region_detect_method = "convex", update_spots = F,
                                   pos_colnm = "Label_Tb927.6.4280", 
                                   description = NULL,grp_nm = "DBSCAN_Tbb_microbe", dir_nm = "M2_NicheDetect_STS")
STID_obj_detect
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_DBSCAN_Tbb_microbe",
                           add_coord = F)[[1]]

#> plot
SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF" )
SEVEN_LIGHT <- c("#F88A7E", "#FCC093")
pdf("./photo/Z_other/20260408_Figure3_result/Tbb/Tbb_microbe_edge_raw.pdf",width = 14,height = 5)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "group", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 2,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Tbb/Tbb_microbe_region_raw.pdf",width = 14,height = 5)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "All_ROI_label2",
             facet_grpnm = "group", datatype = "discrete",
             col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
             pt_size = 2,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

#> Plot_DistLine_Exp
STID_obj_detect %>% print()
tmp <- GetMetaData(STID_obj_detect,
                   meta_key = list(c("M1_SpotDetect_Geneset_Tbb_correct_before_PCD_white","M2_NicheDetect_STS_DBSCAN_Tbb_microbe")))[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/Tbb/Tbb_microbe_DistLine_Exp_host.pdf",width = 4.8,height = 6)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = c("Tb927.6.4280","Chil3"), feature_colnm = NULL, 
                  loop_id = "Inf_d45pi", col = c("#F81B02FF"  ,"#FCB11C" ) ,
                  facet_grpnm = "group",meta_key = "M2_NicheDetect_STS_DBSCAN_Tbb_microbe")
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Tbb/Tbb_microbe_DistLine_Exp_PCD.pdf",width = 4.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Necroptosis"), 
                  loop_id = "Inf_d45pi", col = "#3B95C4FF",
                  facet_grpnm = "group",meta_key = list(c("M1_SpotDetect_Geneset_Tbb_correct_before_PCD_white",
                                                          "M2_NicheDetect_STS_DBSCAN_Tbb_microbe")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Tbb/Tbb_microbe_DistLine_Exp_geneset.pdf",width = 4.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Chagas.disease"), 
                  loop_id = "Inf_d45pi", col = "#3B95C4FF",
                  facet_grpnm = "group",meta_key = list(c("M1_SpotDetect_Geneset_Tbb_correct_before_KEGG_parasite_white",
                                                          "M2_NicheDetect_STS_DBSCAN_Tbb_microbe")))
dev.off()


## 1.2 Tbb：host ####
STID_obj_detect %>% print()
meta_data <- GetMetaData(
  STID_obj_detect,
  meta_key = "M1_SpotDetect_Geneset_Tbb_correct_before_KEGG_parasite_white"
)[[1]]
colnames(meta_data)
STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj_detect, meta_key = "M1_SpotDetect_Geneset_Tbb_correct_before_KEGG_parasite_white", 
                                   spatial_scale_method = "region", region_detect_method = "convex", update_spots = F,
                                   pos_colnm = "Label_Chagas.disease",
                                   description = NULL,grp_nm = "DBSCAN_Tbb_host", dir_nm = "M2_NicheDetect_STS")
STID_obj_detect
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_DBSCAN_Tbb_host",
                           add_coord = F)[[1]]

#> plot
SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF" )
SEVEN_LIGHT <- c("#F88A7E", "#FCC093")
pdf("./photo/Z_other/20260408_Figure3_result/Tbb/Tbb_host_edge_raw.pdf",width = 14,height = 5)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "group", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 2,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Tbb/Tbb_host_region_raw.pdf",width = 14,height = 5)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "All_ROI_label2",
             facet_grpnm = "group", datatype = "discrete",
             col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
             pt_size = 2,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

#> Plot_DistLine_Exp
STID_obj_detect %>% print()
tmp <- GetMetaData(STID_obj_detect,
                   meta_key = list(c("M1_SpotDetect_Geneset_Tbb_correct_before_PCD_white","M2_NicheDetect_STS_DBSCAN_Tbb_host")))[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/Tbb/Tbb_host_DistLine_Exp_host.pdf",width = 4.8,height = 6)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = c("Tb927.6.4280","Chil3"), feature_colnm = NULL, 
                  loop_id = "Inf_d45pi", col = c("#F81B02FF"  ,"#FCB11C" ) ,
                  facet_grpnm = "group",meta_key = "M2_NicheDetect_STS_DBSCAN_Tbb_host")
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Tbb/Tbb_host_DistLine_Exp_PCD.pdf",width = 4.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Necroptosis"), 
                  loop_id = "Inf_d45pi", col = "#3B95C4FF",
                  facet_grpnm = "group",meta_key = list(c("M1_SpotDetect_Geneset_Tbb_correct_before_PCD_white",
                                                          "M2_NicheDetect_STS_DBSCAN_Tbb_host")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Tbb/Tbb_host_DistLine_Exp_geneset.pdf",width = 4.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Chagas.disease"), 
                  loop_id = "Inf_d45pi", col = "#3B95C4FF",
                  facet_grpnm = "group",meta_key = list(c("M1_SpotDetect_Geneset_Tbb_correct_before_KEGG_parasite_white",
                                                          "M2_NicheDetect_STS_DBSCAN_Tbb_host")))
dev.off()


#> CompareNiche
STID_obj_detect
pdf("./photo/Z_other/20260408_Figure3_result/Tbb/Tbb_CompareNiche_microbe_host.pdf",width = 5,height = 8)
CompareNiche(STID_obj = STID_obj_detect, 
             meta_key1 = "M2_NicheDetect_STS_DBSCAN_Tbb_microbe",
             meta_key2 = "M2_NicheDetect_STS_DBSCAN_Tbb_host")
dev.off()
saveRDS(STID_obj_detect,file = "./rds/STID_obj_detect_Tbb.rds")


## 1.3 kp：host ####
STID_obj %>% print()
meta_data <- GetMetaData(
  STID_obj,
  meta_key = "M1_SpotDetect_Geneset_Kp_correct_before_GO_bacterial_white"
)[[1]]
colnames(meta_data)
STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj, meta_key = "M1_SpotDetect_Geneset_Kp_correct_before_GO_bacterial_white",
                                   spatial_scale_method = "region", region_detect_method = "convex", update_spots = F,
                                   pos_colnm = "Label_ANTIBACTERIAL_HUMORAL_RESPONSE",density_thres = 1,
                                   description = NULL,grp_nm = "DBSCAN_kp_host", dir_nm = "M2_NicheDetect_STS")
# STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj, meta_key = "M1_SpotDetect_Gene_Kp_correct_before_host_gene_white",
#                                   spatial_scale_method = "region", region_detect_method = "convex", update_spots = F,
#                                   pos_colnm = "Label_Ccl20", ROI_size = 20,density_thres = 0.9,
#                                   description = NULL,grp_nm = "DBSCAN_kp_host", dir_nm = "M2_NicheDetect_STS")
STID_obj_detect
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_DBSCAN_kp_host",
                           add_coord = F)[[1]]
# SEVEN_DARK = c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C" ,"#AFBF41FF" ,"#50C49FFF" ,"#3B95C4FF" ,"#B560D4FF") 
# SEVEN_LIGHT = c("#F88A7E", "#FCC093", "#FCDB9A", "#BBBFA1", "#B2C4BE", "#9DB7C4", "#D1CAD4")
SEVEN_DARK = c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C"  ,"#50C49FFF" ,"#3B95C4FF" ,"#B560D4FF") 
SEVEN_LIGHT = c("#F88A7E", "#FCC093", "#FCDB9A", "#B2C4BE", "#9DB7C4", "#D1CAD4")
detect_meta <- detect_meta %>% 
  filter(orig.ident == "A3")
pdf("./photo/Z_other/20260408_Figure3_result/Kp/Kp_host_edge_raw.pdf",width = 4.2,height = 4)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y", group_by = "ROI_region",
             facet_grpnm = "orig.ident", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1.5,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Kp/Kp_host_region_raw.pdf",width = 4.2,height = 4)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "All_ROI_label2",
             facet_grpnm = "orig.ident", datatype = "discrete",
             col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
             pt_size = 1.5,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

#> Plot_DistLine_Exp
STID_obj_detect %>% print()
pathogen_genes
host_genes
tmp <- GetMetaData(STID_obj_detect,
                   meta_key = list(c("M2_NicheDetect_STS_DBSCAN_kp_host")))[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/Kp/Kp_host_DistLine_Exp_host.pdf",width = 4.5, height = 5)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = host_genes, feature_colnm = NULL, 
                  loop_id = "A3", col = c("#3B95C4FF"  ,"#3B95C4FF" ) ,
                  facet_grpnm = "orig.ident",meta_key = "M2_NicheDetect_STS_DBSCAN_kp_host")
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Kp/Kp_host_DistLine_Exp_PCD.pdf",width = 6,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("ANTIBACTERIAL_HUMORAL_RESPONSE"), 
                  loop_id = "A3", col = c("#3B95C4FF","#3B95C4FF"),
                  facet_grpnm = "orig.ident",meta_key = list(c("M1_SpotDetect_Geneset_Kp_correct_before_GO_bacterial_white",
                                                               "M2_NicheDetect_STS_DBSCAN_kp_host")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Kp/Kp_host_DistLine_Exp_geneset.pdf",width = 4.5,height = 5)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Apoptosis","Necroptosis"), 
                  loop_id = "A3", col = c("#3B95C4FF","#3B95C4FF"),
                  facet_grpnm = "orig.ident",meta_key = list(c("M1_SpotDetect_Geneset_Kp_correct_before_PCD_white",
                                                               "M2_NicheDetect_STS_DBSCAN_kp_host")))
dev.off()
saveRDS(STID_obj_detect,file = "./rds/STID_obj_detect_kp_host.rds")

## 1.4 Lang：host ####
STID_obj %>% print()
meta_data <- GetMetaData(
  STID_obj,
  meta_key = "M1_SpotDetect_Geneset_Lang_correct_before_GO_viral_white"
)[[1]]
colnames(meta_data)
STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj, meta_key = "M1_SpotDetect_Geneset_Lang_correct_before_GO_viral_white", 
                                   # loop_id = "T1L_D7",
                                   spatial_scale_method = "region", region_detect_method = "convex", update_spots = F,
                                   density_thres = 0.5,ROI_size = 12,
                                   pos_colnm = "Label_RESPONSE_TO_VIRUS",
                                   description = NULL,grp_nm = "DBSCAN_Lang_host", dir_nm = "M2_NicheDetect_STS")
STID_obj_detect
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_DBSCAN_Lang_host",
                           add_coord = F)[[1]]
# SEVEN_DARK = c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C" ,"#AFBF41FF" ,"#50C49FFF" ,"#3B95C4FF" ,"#B560D4FF","grey80")
# SEVEN_LIGHT = c("#F88A7E", "#FCC093", "#FCDB9A", "#BBBFA1", "#B2C4BE", "#9DB7C4", "#D1CAD4","grey95")
SEVEN_DARK = c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C" ,"#AFBF41FF"  ,"#3B95C4FF" ,"#B560D4FF")
SEVEN_LIGHT = c("#F88A7E", "#FCC093", "#FCDB9A", "#BBBFA1", "#9DB7C4", "#D1CAD4")
detect_meta <- detect_meta %>% 
  filter(sample == "T1L_D7")
pdf("./photo/Z_other/20260408_Figure3_result/Lang/Lang_host_edge_raw.pdf",width = 3.5,height = 4)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "sample", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1.5,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Lang/Lang_host_region_raw.pdf",width = 3.5,height = 4)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "All_ROI_label2",
             facet_grpnm = "sample", datatype = "discrete",
             col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
             pt_size = 1.5,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()


#> Plot_DistLine_Exp
STID_obj_detect %>% print()
pathogen_genes
host_genes
tmp <- GetMetaData(STID_obj_detect,
                   meta_key = list(c("M2_NicheDetect_STS_DBSCAN_Lang_host")))[[1]]
pdf("./photo/Z_other/20260408_Figure3_result/Lang/Lang_host_DistLine_Exp_host.pdf",width = 4.2,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = "Ly6a", feature_colnm = NULL, 
                  loop_id = "T1L_D7", col = c("#3B95C4FF" ) ,
                  facet_grpnm = "sample",
                  meta_key = "M2_NicheDetect_STS_DBSCAN_Lang_host")
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Lang/Lang_host_DistLine_Exp_PCD.pdf",width = 5.2,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("RESPONSE_TO_VIRUS"), 
                  loop_id = "T1L_D7", col = c("#3B95C4FF","#3B95C4FF"),
                  facet_grpnm = "sample",meta_key = list(c("M2_NicheDetect_STS_DBSCAN_Lang_host",
                                                           "M1_SpotDetect_Geneset_Lang_correct_before_GO_viral_white")))
dev.off()
pdf("./photo/Z_other/20260408_Figure3_result/Lang/Lang_host_DistLine_Exp_geneset.pdf",width = 4.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Necroptosis"), 
                  loop_id = "T1L_D7", col = c("#3B95C4FF","#3B95C4FF"),
                  facet_grpnm = "sample",meta_key = list(c("M1_SpotDetect_Geneset_Lang_correct_before_PCD_white",
                                                           "M2_NicheDetect_STS_DBSCAN_Lang_host")))
dev.off()
saveRDS(STID_obj_detect,file = "./rds/STID_obj_detect_Lang.rds")


## 1.5 Lang & CE Spot ####
STID_obj_detect %>% print
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M1_SpotDetect_Geneset_Lang_correct_before_GO_viral_white",
                           add_coord = F)[[1]]
colnames(detect_meta)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "Label_RESPONSE_TO_VIRUS",
             facet_grpnm = "sample", datatype = "discrete",
             col = list(dis = c("grey95","#EB1E2C"),con = NULL),
             pt_size = 1.5,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj_detect, meta_key = "M1_SpotDetect_Geneset_Lang_correct_before_GO_viral_white", 
                                   # loop_id = "T1L_D7",
                                   spatial_scale_method = "spot", region_detect_method = "convex", update_spots = F,
                                   density_thres = 0.5,ROI_size = 10,
                                   pos_colnm = "Label_RESPONSE_TO_VIRUS",
                                   description = NULL,grp_nm = "spot_Lang_host", dir_nm = "M2_NicheDetect_STS")
STID_obj_detect
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_spot_Lang_host",
                           add_coord = F)[[1]]
SEVEN_DARK = c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C" ,"#AFBF41FF"  ,"#3B95C4FF" ,"#B560D4FF")
SEVEN_LIGHT = c("#F88A7E", "#FCC093", "#FCDB9A", "#BBBFA1", "#9DB7C4", "#D1CAD4")
detect_meta <- detect_meta %>% 
  filter(sample == "T1L_D7")
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "sample", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1.5,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)


#### 五、Figure4: CE ####
if(0){
  table(STID_obj_expand@meta.data$anno)
  celltypes <- STID_obj_expand@meta.data$anno %>% levels()
  col_lasso <- COLOR_LIST$PALETTE_WHITE_BG[1:10]
  names(col_lasso) <- celltypes
  col_lasso2 <- col_lasso
  col_lasso["Hepatocytes"] <- "grey95"
}else{
  col_lasso <- c(
    "HsPCs" = "#E41A1C",
    "Hepatocytes" = "grey95",
    "Infla Heps" = "#4DAF4A",
    "Fibroblasts" = "#984EA3",
    "Cho/Spp1+ cells" = "#FFFF33",
    "Spp1+ MoMFs" = "#FF7F00",
    "MoKCs" = "#377EB8",
    "Neutrophils" = "#F781BF",
    "B/plasma cells" = "#A65628",
    "Others" = "#8DA0CB"
  )
  col_lasso2 <- c(
    "HsPCs" = "#E41A1C",
    "Hepatocytes" = "#66C2A5",
    "Infla Heps" = "#4DAF4A",
    "Fibroblasts" = "#984EA3",
    "Cho/Spp1+ cells" = "#FFFF33",
    "Spp1+ MoMFs" = "#FF7F00",
    "MoKCs" = "#377EB8",
    "Neutrophils" = "#F781BF",
    "B/plasma cells" = "#A65628",
    "Others" = "#8DA0CB"
  )
}


### 1.BC: CreateSSNiche ####
suppressMessages({
  library(vroom)
})

## 1.1 组织和细胞边界图
geneset_df <- vroom(file = "./inputdata/Gene_Geneset2/Mouse/Geneset/Stratify/Mouse_Stratify_liver_geneset.txt",delim = "\t")
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj_expand <- SpotDetect_Geneset(STID_obj_expand,
                                      geneset_list = geneset_list,
                                      score_method = "MeanExp", n_iter = 5, nbin = 24,seed = 10,
                                      PosThres_prob = 0, PosThres_score = 0,
                                      pt_size = 0.25,
                                      col = COLOR_DIS_CON,
                                      black_bg = F,
                                      blur_method = NULL,
                                      # blur_method = "isoblur",
                                      plot_method = "single",
                                      grp_nm = "CE_detect_geneset_stratify")
STID_obj_expand %>% print()

#> 组织边界图
meta_key <- "M2_NicheExpand_CE"
meta_data <- GetMetaData(STID_obj_expand,meta_key = list(c("M1_SpotDetect_Geneset_CE_detect_geneset_stratify",
                                                           meta_key)))[[1]]
meta_data <- meta_data %>% 
  rownames_to_column(var = "cell_id") %>%
  group_by(batch) %>%
  mutate(PV_CV = Periportal - Pericentral) %>% 
  mutate(label = cut_number(PV_CV, n = 9,labels = F)) %>% 
  mutate(label9 = paste0("Layer_",as.numeric(label))) %>% 
  mutate(label3 = ifelse(label9 %in% c("Layer_1","Layer_2","Layer_3"),"Pericentral",
                         ifelse(label9 %in% c("Layer_4","Layer_5","Layer_6"),"Midzone","Periportal"))) %>% 
  mutate(label3 = factor(label3,levels = c("Pericentral","Midzone","Periportal"))) %>% 
  mutate(merge_label = if_else(ROI_region == "edge","edge",label9,missing = label9)) %>% 
  column_to_rownames(var = "cell_id")
STID_obj_expand <- AddMetaData(STID_obj = STID_obj_expand,meta_key = "hep_layer9",add_data = meta_data)
table(meta_data[c("batch","label3")])
table(meta_data[c("batch","label9")])
table(meta_data[c("batch","label")])
table(meta_data[c("batch","merge_label")])

Plot_Spatial(plot_data = meta_data,x_colnm = "x", y_colnm = "y",group_by = "label3",
             facet_grpnm = "batch", datatype = "discrete",
             col = list(dis = c("#F81B02FF"  ,"#FCDB9A","#9DB7C4" ),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
Plot_Spatial(plot_data = meta_data,x_colnm = "x", y_colnm = "y",group_by = "label9",
             facet_grpnm = "batch", datatype = "discrete",
             col = list(dis = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF", "#5DC863FF", "#AADC32FF", "#FDE725FF") %>% rev(),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_liver_layer_edge_merge.pdf",width = 20,height = 5)
Plot_Spatial(plot_data = meta_data,x_colnm = "x", y_colnm = "y",group_by = "merge_label",
             facet_grpnm = "batch", datatype = "discrete",
             col = list(dis = c("red",c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF", "#5DC863FF", "#AADC32FF", "#FDE725FF") %>% rev()),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()


#> 细胞边界图
celltype_data <- GetMetaData(STID_obj_expand, meta_key = meta_key)[[1]]
celltype_data <- celltype_data %>% 
  mutate(anno_merge = if_else(ROI_region == "edge","edge",anno,missing = anno)) %>% 
  mutate(anno_merge = factor(anno_merge,levels = c("edge",levels(anno))))
table(celltype_data[c("batch","anno_merge")])

pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_celltype_spatial.pdf",width = 20,height = 5)
Plot_Spatial(plot_data = celltype_data,x_colnm = "x", y_colnm = "y",group_by = "anno", # !!!
             facet_grpnm = "batch", datatype = "discrete",
             col = list(dis = c(edge = "red",col_lasso),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_celltype_edge_merge.pdf",width = 20,height = 5)
Plot_Spatial(plot_data = celltype_data,x_colnm = "x", y_colnm = "y",group_by = "anno_merge",
             facet_grpnm = "batch", datatype = "discrete",
             col = list(dis = c(edge = "red",col_lasso),con =NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()


## 1.2 组织和细胞条形图：CreateSingleSampNiche
#> STID_obj_SS
STID_obj_expand %>% print()
meta_data <- GetMetaData(STID_obj_expand,meta_key = list(c(meta_key,"hep_layer9")))[[1]]
colnames(meta_data)
STID_obj_SS <- CreateSingleSampNiche(
  STID_obj = STID_obj_expand,
  # loop_id = "DPI_4_2",
  niche_key = "Niche",
  meta_key = list(c(meta_key,"hep_layer9")),
  ROI_type = "ROI",
  pos_colnm = "ROI_label",
  center_colnm = "ROI_center",
  edge_colnm = "ROI_edge",
  all_label_colnm = "All_ROI_label",
  all_dist_colnm = "All_Dist2ROIcenter",
  other_colnm = "label9",
  description = NULL
)
STID_obj_SS %>% print()
STID_obj_SS <- AddSSNicheCells(
  STID_obj = STID_obj_SS,
  # loop_id = "DPI_4_2", # vector
  meta_key = "hep_layer9", # string
  select_colnm = "label3", # vector
  niche_key = "Niche"
)
SS_cells <- GetSSNicheCells(STID_obj_SS,
                            loop_id = "DPI_4_2",
                            niche_key = "Niche")[[1]]  
saveRDS(STID_obj_SS,file = "./rds/STID_obj_SS_CE.rds")


#> 分层共定位结果：结构
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_liver_layer_niche_composition.pdf",width = 3,height = 3)
CalSampComp(STID_obj = STID_obj_SS, niche_key = "Niche",group_by = "label9",
            loop_id = "LoopAllSamp",
            # col = COLOR_LIST[["PALETTE_9_CLASSIC"]],
            col =c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF", "#5DC863FF", "#AADC32FF", "#FDE725FF" ) %>% rev(),
            return_data = F
)
dev.off()
CalSampComp(STID_obj = STID_obj_SS, niche_key = "Niche",group_by = "label3",
            loop_id = "LoopAllSamp",
            col =  c("#F81B02FF"  ,"#FCDB9A","#9DB7C4" ),
            return_data = F
)

#> 分层共定位结果：细胞
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_celltype_niche_composition.pdf",width = 3,height = 3)
CalSampComp(STID_obj = STID_obj_SS, niche_key = "Niche",group_by = "anno",
            loop_id = "LoopAllSamp",
            col = col_lasso,
            return_data = F)
dev.off()


### 2.DE：免疫富集指数 ####
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_niche_agg_index.pdf",width = 6,height = 5)
tmp <- CalSampCAI(STID_obj = STID_obj_SS,
                  niche_key = "Niche",
                  group_by = "anno",
                  loop_id = "LoopAllSamp",
                  k_neighbors = 8,
                  min_agg_size = 10,
                  dist_thres = 1,
                  col = col_lasso)
dev.off()

### 3.F: 细胞分布图 ####
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_celltype_DistLine_Ratio.pdf",width = 4,height = 8)
Plot_DistLine_Ratio(STID_obj = STID_obj_SS,
                    celltypes = c("Hepatocytes","Spp1+ MoMFs","Neutrophils","Fibroblasts","Cho/Spp1+ cells"),
                    group_by = "anno",
                    loop_id = "DPI_4_2",
                    facet_grpnm = "batch",meta_key = meta_key,coord_interval_ratio = 8)
dev.off()


### 3.+ 共定位图 ####
tmp <- CalSampCoLoc(STID_obj = STID_obj_SS,
                    niche_key = "Niche",
                    group_by = "anno",
                    # group_use = c(" Hepatocytes","Infla Heps","Fibroblasts","Spp1+ MoMFs","MoKCs","Neutrophils","B/plasma cells"),
                    loop_id = "DPI_4_2",
                    method = "mistyR",
                    mistyR_params = c(15,10),
                    return_data = TRUE)

# 基因共定位，只有mistyR可用
if(0){
  tmp <- CalSampCoLoc(STID_obj = STID_obj_SS,
                      niche_key = "Niche",
                      group_by = NULL,
                      features = c("Ccl3","Ccl4","Ccr1","Cxcl3","Cxcr2","Il1b","Spp1","Icam1","Itgb2"), # 基因共定位
                      feature_colnm = NULL,
                      loop_id = "DPI_4_2",
                      return_data = TRUE) 
}
configure_conda("squidpy")
tmp <- CalSampCoLoc(STID_obj = STID_obj_SS, 
                    # niche_key = "Niche",
                    group_by = "anno",
                    # group_use = c(" Hepatocytes","Infla Heps","Fibroblasts","Spp1+ MoMFs","MoKCs","Neutrophils","B/plasma cells"),
                    loop_id = "DPI_4_2",
                    method = "squidpy",
                    return_data = TRUE)

#> 空间图共定位
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_celltype_coloc.pdf",width = 5,height = 5)
Plot_SpatialCoLoc(STID_obj = STID_obj_SS,
                  niche_key = "Niche",
                  group_by = "anno",
                  group_use = c("Neutrophils","Spp1+ MoMFs"),
                  loop_id = "LoopAllSamp")
dev.off()
Plot_SpatialCoLoc(STID_obj = STID_obj_SS,
                  niche_key = "Niche",
                  group_by = "anno",
                  group_use = c("Neutrophils","Infla Heps"),
                  loop_id = "LoopAllSamp")
Plot_SpatialCoLoc(STID_obj = STID_obj_SS,
                  niche_key = "Niche",
                  group_by = "anno",
                  group_use = c("Spp1+ MoMFs","Hepatocytes"),
                  loop_id = "LoopAllSamp")
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_Ccl3_ccr1_coloc.pdf",width = 5,height = 5)
Plot_SpatialCoLoc(STID_obj = STID_obj_SS,
                  niche_key = "Niche",
                  group_by = NULL,
                  group_use = NULL,
                  features = c("Ccl3","Ccr1"),
                  loop_id = "LoopAllSamp")
dev.off()
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_Cxcl3_cxcr2_coloc.pdf",width = 5,height = 5)
Plot_SpatialCoLoc(STID_obj = STID_obj_SS,
                  niche_key = "Niche",
                  group_by = NULL,
                  group_use = NULL,
                  features = c("Cxcl3","Cxcr2"),
                  loop_id = "LoopAllSamp")
dev.off()
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_Ccl3_ccr5_coloc.pdf",width = 5,height = 5)
Plot_SpatialCoLoc(STID_obj = STID_obj_SS,
                  niche_key = "Niche",
                  group_by = NULL,
                  group_use = NULL,
                  features = c("Ccl3","Ccr5"),
                  loop_id = "LoopAllSamp")
dev.off()


#> 寄生虫基因与宿主基因共定位
# 没有任何基因是和寄生虫比较高度共定位的
if(0){
  loop_genes <- tmp_data %>% 
    arrange(-Correlation) %>% 
    filter(Correlation > 0) %>% 
    pull(Variable2)
  for (i_gene in loop_genes) {
    Plot_SpatialCoLoc(STID_obj = STID_obj_SS,
                      niche_key = "Niche",
                      group_by = NULL,
                      group_use = NULL,
                      features = i_gene,
                      feature_colnm = "all_gene_nCount(sum)",
                      exp_thres = 10,
                      loop_id = "DPI_4_2")
  }
  
  #>
  STID_obj_SS %>% print()
  test <- STID_obj_SS
  meta_data <- GetMetaData(STID_obj_SS, meta_key = list(c("M1_SpotDetect_Gene_CE_correct_after_host_gene_white",
                                                          meta_key)))[[1]]
  test@meta.data <- meta_data
  test <- subset(test, subset = !is.na(ROI_label))
  colnames(meta_data)
  
  degs_niche <- FindMarkers(object = test,
                            ident.1 = "neg",
                            ident.2 = "pos",
                            group.by = "Label_all_gene_nFeature(sum)",
                            only.pos = FALSE)
  degs_niche_list <- .Process_DEGs(degs_niche, grp_id = "Overall", logfc_thres = 0.05, 
                                   padj_thres = 1, topGeneN = 10)
  degs_niche <- degs_niche_list$degs
  tmp <- degs_niche %>% 
    filter(change == "Up") %>% 
    filter(pct.2<0.1) %>% 
    mutate(diff_pct = pct.1 - pct.2) %>% 
    arrange(-diff_pct)
  loop_genes <- tmp$gene %>% head(40)
  for (i_gene in loop_genes) {
    Plot_SpatialCoLoc(STID_obj = STID_obj_SS,
                      niche_key = "Niche",
                      group_by = NULL,
                      group_use = NULL,
                      features = i_gene,
                      feature_colnm = "all_gene_nCount(sum)",
                      exp_thres = 2,
                      loop_id = "DPI_4_2")
  }
}


### 4.GHIJ: CalSampDEGs ####
# 去掉非mRNA
parse_gtf <- fread(file = "./inputdata/gencode.vM35.annotation.gtf_parsed.txt",sep = "\t",header = T)
mRNA_gene <- parse_gtf$gene_name[parse_gtf$gene_type == "protein_coding"]
mRNA_gene <- unique(mRNA_gene)
head(mRNA_gene)
length(mRNA_gene)
table(mRNA_gene %in% rownames(STID_obj_SS))
pathogen_genes <- grep("EmuJ-",rownames(stRNA),value = T)

nrow(STID_obj_SS) # 基因数
STID_obj_SS <- subset(STID_obj_SS,features = c(mRNA_gene,pathogen_genes))
nrow(STID_obj_SS) 

STID_obj_SS
SS_DEGs <- CalSampDEGs(STID_obj = STID_obj_SS,
                       niche_key = "Niche",
                       assay_id = "Spatial",
                       layer_id = "counts",
                       loop_id = "LoopAllSamp",
                       padj_thres = 0.05,
                       logfc_thres = 1, # log2(1.5)不好
                       adjust_method = "BH",
                       col = col_lasso2,
                       group_by = "anno",
                       group_value = c("Infla Heps","Fibroblasts","Neutrophils","MoKCs"), # "Spp1+ MoMFs"不显著
                       remove_genes = c(grep("^Gm",rownames(STID_obj_SS),value = T))
)
colnames(SS_DEGs$DPI_4_2$data$Overall_DEGs)
Overall_DEGs <- SS_DEGs$DPI_4_2$data$Overall_DEGs
DEGs_gene <- SS_DEGs$DPI_4_2$data$Overall_DEGs %>% 
  filter(change %in% c("Up","Down")) %>%
  # filter(avg_log2FC > 0.25 | avg_log2FC < -0.25) %>%
  filter(p_val_adj<0.05) %>%
  pull(gene) %>% unique()
# DEGs_gene <- SS_DEGs$DPI_4_2$data$Overall_DEGs$gene %>% unique()
length(DEGs_gene)
c("Ccl3","Ccl4","Ccr5") %in% DEGs_gene


#> 富集图
up_gene <- SS_DEGs$DPI_4_2$data$Overall_DEGs$gene[SS_DEGs$DPI_4_2$data$Overall_DEGs$change == "Up"]
down_gene <- SS_DEGs$DPI_4_2$data$Overall_DEGs$gene[SS_DEGs$DPI_4_2$data$Overall_DEGs$change == "Down"]
length(up_gene);length(down_gene)
tmp2 <- GeneEnrichment(STID_obj = STID_obj_SS,
                       DEGs = NULL,
                       up_gene = up_gene,
                       down_gene = down_gene,
                       method = "GO_KEGG", # GO_KEGG, 
                       go_ont = "BP",
                       return_data = T)
if(0){
  # 太耗时了
  tmp2 <- GeneEnrichment(STID_obj = STID_obj_SS,
                         DEGs = Overall_DEGs,
                         up_gene = NULL,
                         down_gene = NULL,
                         method = "GSEA_GO_KEGG", # GSEA
                         go_ont = "BP",
                         return_data = T)
}


#> 展示宿主基因和通路
STID_obj_SS <- SpotDetect_Gene(STID_obj_SS,
                               features = c("Cxcl3","Ccl3","Cxcr2","Ccr1","Ccr5","Il36g","Cstdc6","L3mbtl1","Tsx","Scube2"),
                               feature_colnm = NULL,
                               PosThres_prob = 0, PosThres_count = 3, # D4_2为3，D15_2为10
                               # col = c("grey95","red"),
                               col = COLOR_DIS_CON,
                               # col = c("grey95","#6A1B9A"),
                               # col = c("grey95","#6A1B9A"),
                               # col = c("grey95","#512DA8"),
                               black_bg = F,pt_size = 1,
                               blur_method = NULL, blur_n = 1,blur_sigma = 0.5, 
                               plot_method = "single",
                               grp_nm = "host_Cxcl3")
go_terms <- c(
  "GO:0050900","GO:0060326","GO:0097529","GO:0002274","GO:0050727",
  "GO:0030595","GO:0032103","GO:0019221","GO:0022407","GO:0007159"
)
geneset_df <- GOID2Genelist(GOID = go_terms, STID_obj = STID_obj_SS, host_org = NULL,
                            ont = "BP", keyType = "SYMBOL")
colnames(geneset_df)
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
# geneset_list <- geneset_list[4]
names(geneset_list)
STID_obj_SS <- SpotDetect_Geneset(STID_obj_SS,
                                  geneset_list = geneset_list,
                                  score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                                  PosThres_prob = 0.25, PosThres_score = 0, # ！！！！调低PosThres_prob？
                                  pt_size = 1,
                                  col = COLOR_DIS_CON,
                                  black_bg = F, blur_method = NULL,
                                  plot_method = "single",
                                  grp_nm = "host_inflam_response")
STID_obj_SS %>% print()
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_DEGs_gene_DistLine_Exp.pdf",width = 4,height = 9)
Plot_DistLine_Exp(STID_obj = STID_obj_SS, features = c("Cxcl3","Ccl3","Cxcr2","Ccr1","Ccr5"), feature_colnm = NULL, 
                  loop_id = "DPI_4_2",
                  col = c("#377EB8" ,"#377EB8","#377EB8","#377EB8","#377EB8" ) ,
                  # facet_grpnm = "grp",
                  meta_key = list(c("M1_SpotDetect_Gene_host_Cxcl3",
                                    meta_key)))
dev.off()
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_DEGs_GO_DistLine_Exp.pdf",width = 5.5,height = 5)
Plot_DistLine_Exp(STID_obj = STID_obj_SS, features = NULL, 
                  # feature_colnm = colnames(geneset_df),
                  feature_colnm = c("myeloid leukocyte migration","myeloid leukocyte activation"),
                  loop_id = "DPI_4_2",
                  col = c("#377EB8" ,"#377EB8","#377EB8","#377EB8","#377EB8","#377EB8" ,"#377EB8","#377EB8","#377EB8","#377EB8","#377EB8" ,"#377EB8","#377EB8","#377EB8","#377EB8" ) ,
                  # facet_grpnm = "grp",
                  meta_key = list(c("M1_SpotDetect_Geneset_host_inflam_response",
                                    meta_key)))
dev.off()

### 5.IJ:cellchat ####
#> cellchat
tmp <- CalSampCellComm(STID_obj = STID_obj_SS,
                       niche_key = "Niche",
                       group_by = NULL,
                       assay_id = "Spatial",
                       layer_id = "counts", # ??
                       loop_id = "LoopAllSamp",
                       col = col_lasso2,
                       is_Spatial = T, # 是否进行空间cellchat，默认TRUE
                       spatial.factors = NULL, 
                       interaction.range = 250,
                       return_data = TRUE,
                       grp_nm = NULL,dir_nm = "M3_CalSampCellComm"
)

CellComm_data <- readRDS("./rds/STID_CellComm_data_CE.rds")
pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_cellchat_plot.pdf",width = 5,height = 4.5)
Plot_NicheCellComm(
  STID_obj = STID_obj_SS,
  CellComm_data = CellComm_data,
  samp_mode = "SS",
  loop_id = "LoopAllSamp",
  sources.use = c("Infla Heps","Neutrophils","Spp1+ MoMFs"),
  targets.use = c("Neutrophils","Spp1+ MoMFs","MoKCs","Fibroblasts"),
  # signaling = c("CXCL","CCL","SAA","IL1","COMPLEMENT","TNF","ANNEXIN","ANGPTL"),
  signaling = c("CXCL","CCL","IL1"),
  pairLR.use = NULL,
  col = col_lasso2
)
dev.off()

pdf("./photo/Z_other/20260413_Figure4_result/CE/CE_cellchat_plot_all.pdf",width = 6,height = 6)
Plot_NicheCellComm(
  STID_obj = STID_obj_SS,
  CellComm_data = CellComm_data,
  samp_mode = "SS",
  loop_id = "LoopAllSamp",
  # sources.use = c("Hepatocytes","Infla Heps","Neutrophils"),
  # targets.use = c("Neutrophils","Spp1+ MoMFs"),
  # signaling = c("CXCL","CCL","SAA","IL1","COMPLEMENT","TNF","ANNEXIN","ANGPTL"),
  signaling = c("CXCL","CCL","SAA","IL1","COMPLEMENT","TNF"),
  pairLR.use = NULL,
  col = col_lasso2
)
dev.off()


### 6.NicheNet ####
ligand_target_matrix <- readRDS("./inputdata/Demo/CalSampGRN/NicheNet_V2/ligand_target_matrix_nsga2r_final_mouse.rds")
lr_network <- readRDS("./inputdata/Demo/CalSampGRN/NicheNet_V2/lr_network_mouse_21122021.rds")
weighted_networks <- readRDS("./inputdata/Demo/CalSampGRN/NicheNet_V2/weighted_networks_nsga2r_final_mouse.rds")
gr_network <- readRDS("./inputdata/Demo/CalSampGRN/NicheNet_V2/gr_network_mouse_21122021.rds")
sig_networks <- readRDS("./inputdata/Demo/CalSampGRN/NicheNet_V2/signaling_network_mouse_21122021.rds")
ligand_tf_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final_mouse.rds"))
ref_data <- list(
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks,
  gr_network = gr_network,
  sig_networks = sig_networks,
  ligand_tf_matrix = ligand_tf_matrix
)

#>
Overall_DEGs <- read.table(file = "F:/Scientific_research/ShiXiaoFeng-InfectiousST/01_Project/2025/20250124_STID/Analysis/outputdata/M3_CalSampDEGs/20260414_215238/DPI_4_2_Niche_vs_Bystander_DEGs.txt",
                           sep = "\t",header = T)
DEGs_gene <- Overall_DEGs %>% 
  filter(avg_log2FC > 1) %>%
  filter(p_val_adj<0.05) %>%
  pull(gene) %>% unique()
length(DEGs_gene)
res_NicheGRN <- CalSampGRN(STID_obj = STID_obj_SS,
                           niche_key = "Niche",
                           group_by = NULL,
                           ref_data = ref_data,
                           sender_celltypes = c("Infla Heps","Neutrophils","Spp1+ MoMFs"),
                           # sender_celltypes = c("Neutrophils"),
                           receiver_celltypes = c("Neutrophils","Spp1+ MoMFs","MoKCs","Fibroblasts"), # 不应该限制输出细胞类型，因为DEGs_gene是整个Niche的
                           # receiver_celltypes = c("Spp1+ MoMFs"),
                           target_features = DEGs_gene,
                           expression_pct = 0.1,
                           top_ligand_num = 10,
                           top_target_num = 10,
                           loop_id = "DPI_4_2",
                           return_data = TRUE
)
tmp <- res_NicheGRN$DPI_4_2$data$ligand_activities


### 7.K: cor ####
STID_obj_SS %>% print()
high_exp_genes <- GetTopGenes(STID_obj_SS, top_n = 7, pattern = "EmuJ-",
                              grp_by_samp = F, grp_by_celltype = F,
                              assay_id = "Spatial", layer_id = "counts")
STID_obj_SS <- AddSSNicheCells(
  STID_obj = STID_obj_SS,
  meta_key = "M1_SpotDetect_Gene_CE_correct_after_host_gene_white", # string
  select_colnm = "all_gene_nCount(sum)", # vector
  niche_key = "Niche"
)
tmp <- GetSSNicheCells(STID_obj_SS,loop_id = "LoopAllSamp",niche_key = "Niche")[[1]]

#> 不扩大DEGs，得到的正相关的基因特别少，没法做富集分析
Overall_DEGs <- read.table(file = "F:/Scientific_research/ShiXiaoFeng-InfectiousST/01_Project/2025/20250124_STID/Analysis/outputdata/M3_CalSampDEGs/20260414_215238/DPI_4_2_Niche_vs_Bystander_DEGs.txt",
                           sep = "\t",header = T)
DEGs_gene <- Overall_DEGs %>% 
  filter(avg_log2FC > 1 | avg_log2FC < -1) %>%
  filter(p_val_adj<0.05) %>%
  pull(gene) %>% unique()
DEGs_up <- Overall_DEGs %>% 
  filter(avg_log2FC > 1) %>%
  filter(p_val_adj<0.05) %>%
  pull(gene) %>% unique()
DEGs_down <- Overall_DEGs %>% 
  filter(avg_log2FC < -1) %>%
  filter(p_val_adj<0.05) %>%
  pull(gene) %>% unique()
length(DEGs_gene);length(DEGs_up);length(DEGs_down)

res_GeneCor <- CalSampGeneCor(STID_obj = STID_obj_SS,
                              loop_id = "DPI_4_2",
                              niche_key = "Niche", # 需要使用Niche，否则基因少
                              # meta_key = "M1_SpotDetect_Gene_CE_correct_after_host_gene_white",
                              p.features = high_exp_genes,
                              # h.features = grep("^Gm",rownames(STID_obj_SS),value = T,invert = T),
                              h.features = DEGs_gene,
                              p.feature_colnm = "all_gene_nCount(sum)",
                              method = "spearman",
                              assay_id = "Spatial",
                              layer_id = "counts",
                              return_data = TRUE
)

#> 正相关主要富集的通路
tmp <- res_GeneCor$DPI_4_2$data
sign_data_up <- tmp3$DPI_4_2$data %>% 
  filter(`P-value` < 0.05 & Correlation >0) %>%
  group_by(Variable1) %>%
  arrange(desc(Correlation)) %>%
  slice_head(n = 3)
sign_data_down <- tmp3$DPI_4_2$data %>% 
  filter(`P-value` < 0.05 & Correlation <0) %>% 
  group_by(Variable1) %>%
  arrange(Correlation) %>%
  slice_head(n = 3)

# 这里的一个bug就是这些和病原正相关的基因，可能是在Niche中下调的
# 可以这么解释，这些基因的下降主要是结构破坏导致的基础功能丧失，但是在丧失功能的Niche中，仍然会有一些基因和
# 病原表达正相关，说明这些基因的表达和病原的表达是同步的，可能是病原导致了这些宿主响应的变化
tmp_data <- tmp3$DPI_4_2$data %>% 
  filter(Variable1 == "all_gene_nCount(sum)")
up_gene <- tmp_data %>% arrange(-Correlation) %>% 
  filter(Correlation > 0 & `P-value` < 0.05) %>%
  # filter(Variable2 %in% DEGs_up) %>%
  head(100) %>%
  pull(Variable2)
up_gene[up_gene %in% DEGs_up] # 没有正相关的基因在DEGs的up中
down_gene <- tmp_data %>% arrange(Correlation) %>%
  filter(Correlation < 0 & `P-value` < 0.05) %>%
  # filter(Variable2 %in% DEGs_down) %>%
  head(100) %>%
  pull(Variable2)
length(up_gene);length(down_gene)
tmp2 <- GeneEnrichment(up_gene = up_gene,
                       down_gene = down_gene,
                       host_org = "mouse",
                       enrichment_method = "GO_KEGG", # GO_KEGG, 
                       go_ont = "BP",
                       return_data = FALSE,
                       grp_nm = "CalSampGeneCor")

### 8.L: PPI ####
high_exp_genes <- GetTopGenes(STID_obj_SS, top_n = 20, pattern = "EmuJ-",
                              grp_by_samp = F, grp_by_celltype = F,
                              assay_id = "Spatial", layer_id = "counts")
p.fasta_path <- "./inputdata/Demo/CalSampPPI/EmuJ_protein.fa"
p.symbol2protein_path <- "./inputdata/Demo/CalSampPPI/EmuJ_symbol2protein.txt"
h.fasta_path <- "./inputdata/Demo/CalSampPPI/mouse_protein.fa"
h.symbol2protein_path <- "./inputdata/Demo/CalSampPPI/mouse_symbol2protein.txt"

DEGs_gene <- Overall_DEGs %>% 
  filter(abs(avg_log2FC) > 1) %>% # log2(1.5
  filter(p_val_adj<0.05) %>%
  pull(gene) %>% unique()
length(DEGs_gene)

res_NichePPI <- CalSampPPI(
  STID_obj = STID_obj_SS,
  p.fasta_path = p.fasta_path, h.fasta_path = h.fasta_path,
  p.symbol2protein_path = p.symbol2protein_path, h.symbol2protein_path = h.symbol2protein_path,
  p.features = high_exp_genes,
  h.features = DEGs_gene,
  BLAST_tool = "blastp", dbtype = "prot",
  BLAST_args = "-evalue 1e-10 -max_target_seqs 5 -max_hsps 50 -num_threads 4",
  BLAST_fil_args = "pident > 40, evalue < 1e-10, bitscore>100,length > 100",
  version_STRING = "12.0",
  score_thre_only = 100,score_thre_all = 800
)
PPI_gene <- grep("p\\.",res_NichePPI$data$PPI_network_all$nodes$name,value = T,invert = T)
length(PPI_gene)
GeneEnrichment( STID_obj = STID_obj_SS,
                up_gene = PPI_gene,
                down_gene = NULL,
                enrichment_method = "GO_KEGG", # GO_KEGG, 
                go_ont = "BP",
                # col = viridis(100, option = "H"),
                return_data = FALSE,
                grp_nm = "PPI_gene")


#### 六、supp Figure3: CE + JEV ####
### 1.A-H: CE ####
# 这一部分属于主图的附图，主要代码在主图


### 2.JEV ####
# 这一部分完全是与主图独立的分析
STID_obj_detect %>% print()

#> 组织颜色
col_lasso <- c('OLF' = '#FFDEAD', 'CTX_HPF' = '#8FBC8F', 'HPF' = '#A0522D', 
               'TH' = '#7FFFAA', 'HY' = '#FFC0CB', 'CNU' = '#FF8C00', 
               'MB' = '#000080', 'HB' = '#9932CC', 'CB' = '#87CEEB', 
               'FB' = '#FFFF00', 'MEN' = '#FF0000', 'CHP' = '#006400', 
               'UK' = '#D2B48C', 'MB_HY' = '#7B68EE', 'MB_HY_HB' = '#EF6FD0', 
               'HB_CB' = '#A9A9A9', 'MB_HB' = '#FF1493', 'CTX' = '#80AD80', 
               'MB_CNU' = '#C1F1A4', 'MB_CNU_HB' = '#B6A7E2', 'TH_HY' = '#F7D89A')
col_lasso2 <- c('OLF' = '#FFDEAD', 'CTX_HPF' = '#8FBC8F', 'HPF' = '#A0522D', 
                'TH' = '#7FFFAA', 'HY' = '#FFC0CB', 'CNU' = '#FF8C00', 
                'MB' = '#000080', 'HB' = '#9932CC', 'CB' = '#87CEEB', 
                'FB' = '#FFFF00', 'MEN' = '#FF0000', 'CHP' = '#006400', 
                'UK' = '#D2B48C', 'MB_HY' = '#7B68EE', 'MB_HY_HB' = '#EF6FD0', 
                'HB_CB' = '#A9A9A9', 'MB_HB' = '#FF1493', 'CTX' = '#80AD80', 
                'MB_CNU' = '#C1F1A4', 'MB_CNU_HB' = '#B6A7E2', 'TH_HY' = '#F7D89A')


#> 细胞颜色
col_lasso <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#66C2A5",
               "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
names(col_lasso) <- unique(STID_obj_detect@meta.data$new_cell) %>% sort()
col_lasso2 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#F781BF", "#66C2A5",
                "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
names(col_lasso2) <- c(
  "Adipocytes", "Astrocytes", "Dendritic cells", "Endothelial cells", 
  "Epithelial cells", "Fibroblasts", "Macrophages", "Microglia", 
  "Monocytes", "Neurons", "NK cells", "Oligodendrocytes", "T cells"
)


# >>> STID_obj_SS
STID_obj_detect %>% print()
meta_data <- GetMetaData(STID_obj_detect,meta_key = list(c("M2_NicheDetect_STS_STS_JEV_microbe_region")))[[1]]
colnames(meta_data)
STID_obj_SS <- CreateSingleSampNiche(
  STID_obj = STID_obj_detect,
  # loop_id = "DPI_4_2",
  niche_key = "Niche_microbe",
  meta_key = list(c("M2_NicheDetect_STS_STS_JEV_microbe_region")),
  ROI_type = "ROI",
  pos_colnm = "ROI_label",
  center_colnm = "ROI_center",
  edge_colnm = "ROI_edge",
  all_label_colnm = "All_ROI_label",
  all_dist_colnm = "All_Dist2ROIcenter",
  description = NULL
)
STID_obj_SS <- CreateSingleSampNiche(
  STID_obj = STID_obj_SS,
  # loop_id = "DPI_4_2",
  niche_key = "Niche_host",
  meta_key = list(c("M2_NicheDetect_STS_STS_JEV_host_region")),
  ROI_type = "ROI",
  pos_colnm = "ROI_label",
  center_colnm = "ROI_center",
  edge_colnm = "ROI_edge",
  all_label_colnm = "All_ROI_label",
  all_dist_colnm = "All_Dist2ROIcenter",
  description = NULL
)
STID_obj_SS <- AddSSNicheCells(
  STID_obj = STID_obj_SS,
  meta_key = "raw", # string
  select_colnm = "new_tissue", # vector
  niche_key = "Niche_microbe"
)
STID_obj_SS <- AddSSNicheCells(
  STID_obj = STID_obj_SS,
  meta_key = "raw", # string
  select_colnm = "new_tissue", # vector
  niche_key = "Niche_host"
)
STID_obj_SS %>% print()
SS_cells <- GetSSNicheCells(STID_obj_SS,niche_key = "Niche_microbe")[[1]]  
table(SS_cells$new_tissue)
unique(SS_cells$new_tissue) %in% names(col_lasso)


#> 分层共定位结果：tissue
meta_data <- GetMetaData(STID_obj_SS, meta_key = list(c("raw")))[[1]]
meta_microbe <- GetMetaData(STID_obj_SS, meta_key = list(c("M2_NicheDetect_STS_STS_JEV_microbe_region")))[[1]]
meta_host <- GetMetaData(STID_obj_SS, meta_key = list(c("M2_NicheDetect_STS_STS_JEV_host_region")))[[1]]
pdf("./photo/Z_other/20260413_Figure4_result/JEV/JEV_microbe_host_Spatial_tissue.pdf",width = 15,height = 15)
Plot_Spatial(plot_data = meta_data,x_colnm = "x", y_colnm = "y",group_by = "new_tissue",
             # facet_grpnm = "batch", 
             datatype = "discrete",
             col = list(dis = alpha(col_lasso,0.2),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
Plot_Spatial(plot_data = meta_microbe,x_colnm = "x", y_colnm = "y",group_by = "ROI_edge",
             datatype = "discrete",
             col = list(dis = c("#0000FF00","red"),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
Plot_Spatial(plot_data = meta_host,x_colnm = "x", y_colnm = "y",group_by = "ROI_edge",
             datatype = "discrete",
             col = list(dis = c("#0000FF00","blue"),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()

#>
pdf("./photo/Z_other/20260413_Figure4_result/JEV/JEV_microbe_host_NicheComposition_tissue.pdf",width = 5,height = 5)
CalSampComp(STID_obj = STID_obj_SS, niche_key = "Niche_microbe",group_by = "new_tissue",
            loop_id = "LoopAllSamp",
            col =col_lasso,
            return_data = F
)
CalSampComp(STID_obj = STID_obj_SS, niche_key = "Niche_host",group_by = "new_tissue",
            loop_id = "LoopAllSamp",
            col =col_lasso,
            return_data = F
)
dev.off()

#> 分层共定位结果：cell
pdf("./photo/Z_other/20260413_Figure4_result/JEV/JEV_microbe_host_Spatial_cell.pdf",width = 15,height = 15)
Plot_Spatial(plot_data = meta_data,x_colnm = "x", y_colnm = "y",group_by = "new_cell",
             # facet_grpnm = "batch", 
             datatype = "discrete",
             col = list(dis = alpha(col_lasso,0.2),con = NULL),
             pt_size = 1.1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
pdf("./photo/Z_other/20260413_Figure4_result/JEV/JEV_microbe_host_NicheComposition_cell.pdf",width = 5,height = 5)
CalSampComp(STID_obj = STID_obj_SS, niche_key = "Niche_microbe",group_by = "new_cell",
            loop_id = "LoopAllSamp",
            col =col_lasso,
            return_data = F
)
CalSampComp(STID_obj = STID_obj_SS, niche_key = "Niche_host",group_by = "new_cell",
            loop_id = "LoopAllSamp",
            col =col_lasso,
            return_data = F
)
dev.off()

## 免疫富集指数 
pdf("./photo/Z_other/20260413_Figure4_result/JEV/JEV_microbe_host_NicheAggIndex_cell.pdf",width = 5,height = 5)
tmp <- CalSampCAI(STID_obj = STID_obj_SS,
                  niche_key = "Niche_microbe",
                  group_by = "new_cell",
                  loop_id = "LoopAllSamp",
                  k_neighbors = 8,
                  min_agg_size = 10,
                  dist_thres = 1,
                  col = col_lasso)
tmp <- CalSampCAI(STID_obj = STID_obj_SS,
                  niche_key = "Niche_host",
                  group_by = "new_cell",
                  loop_id = "LoopAllSamp",
                  k_neighbors = 8,
                  min_agg_size = 10,
                  dist_thres = 1,
                  col = col_lasso)
dev.off()


# 细胞分布图
pdf("./photo/Z_other/20260413_Figure4_result/JEV/JEV_microbe_host_DistLine_cell.pdf",width = 4.5,height = 8)
Plot_DistLine_Ratio(STID_obj = STID_obj_SS,
                    celltypes = c("Dendritic cells","Monocytes","NK cells","Oligodendrocytes"),
                    group_by = "new_cell",
                    # loop_id = "DPI_4_2",
                    col = col_lasso,
                    meta_key = "M2_NicheDetect_STS_STS_JEV_microbe_region",coord_interval_ratio = 15)
Plot_DistLine_Ratio(STID_obj = STID_obj_SS,
                    celltypes = c("Dendritic cells","Monocytes","NK cells","Oligodendrocytes"),
                    group_by = "new_cell",
                    # loop_id = "DPI_4_2",
                    col = col_lasso,
                    meta_key = "M2_NicheDetect_STS_STS_JEV_host_region",coord_interval_ratio = 12)
dev.off()


# 共定位图
tmp <- CalSampCoLoc(STID_obj = STID_obj_SS,
                    niche_key = "Niche_microbe",
                    group_by = "new_cell",
                    # group_use = c(" Hepatocytes","Infla Heps","Fibroblasts","Spp1+ MoMFs","MoKCs","Neutrophils","B/plasma cells"),
                    # loop_id = "DPI_4_2",
                    method = "mistyR",
                    mistyR_params = c(juxtaview_radius = 15, paraview_radius = 10),
                    return_data = TRUE)
tmp <- CalSampCoLoc(STID_obj = STID_obj_SS,
                    niche_key = "Niche_host",
                    group_by = "new_cell",
                    # group_use = c(" Hepatocytes","Infla Heps","Fibroblasts","Spp1+ MoMFs","MoKCs","Neutrophils","B/plasma cells"),
                    # loop_id = "DPI_4_2",
                    method = "mistyR",
                    mistyR_params = c(juxtaview_radius = 15, paraview_radius = 10),
                    return_data = TRUE)

#>
configure_conda("squidpy")
tmp <- CalSampCoLoc(STID_obj = STID_obj_SS,
                    niche_key = "Niche_microbe",
                    group_by = "new_cell",
                    # group_use = c(" Hepatocytes","Infla Heps","Fibroblasts","Spp1+ MoMFs","MoKCs","Neutrophils","B/plasma cells"),
                    # loop_id = "DPI_4_2",
                    method = "squidpy",
                    mistyR_params = c(juxtaview_radius = 15, paraview_radius = 10),
                    return_data = TRUE,
                    grp_nm = "JEV_microbe_Niche_squidpy")
tmp <- CalSampCoLoc(STID_obj = STID_obj_SS,
                    niche_key = "Niche_host",
                    group_by = "new_cell",
                    # group_use = c(" Hepatocytes","Infla Heps","Fibroblasts","Spp1+ MoMFs","MoKCs","Neutrophils","B/plasma cells"),
                    # loop_id = "DPI_4_2",
                    method = "squidpy",
                    mistyR_params = c(juxtaview_radius = 15, paraview_radius = 10),
                    return_data = TRUE,
                    grp_nm = "JEV_host_Niche_squidpy")


# CalSampDEGs
pdf("./photo/Z_other/20260413_Figure4_result/JEV/JEV_microbe_host_DEGs_cell.pdf",width = 5,height = 5)
SS_DEGs_microbe <- CalSampDEGs(STID_obj = STID_obj_SS,
                               niche_key = "Niche_microbe",
                               assay_id = "Spatial",
                               layer_id = "counts",
                               loop_id = "LoopAllSamp",
                               padj_thres = 0.05,
                               logfc_thres = 1, # log2(1.5)不好
                               adjust_method = "BH",
                               col = col_lasso,
                               group_by = "new_cell",
                               group_value = c("Dendritic cells","Monocytes","NK cells","Oligodendrocytes"),
                               remove_genes = c(grep("^Gm",rownames(STID_obj_SS),value = T))
)
dev.off()

pdf("./photo/Z_other/20260413_Figure4_result/JEV/JEV_microbe_host_DEGs_cell.pdf",width = 5,height = 5)
SS_DEGs_host <- CalSampDEGs(STID_obj = STID_obj_SS,
                            niche_key = "Niche_host",
                            assay_id = "Spatial",
                            layer_id = "counts",
                            loop_id = "LoopAllSamp",
                            padj_thres = 0.05,
                            logfc_thres = 1, # log2(1.5)不好
                            adjust_method = "BH",
                            col = col_lasso,
                            group_by = "new_cell",
                            group_value = c("Dendritic cells","Monocytes","NK cells","Oligodendrocytes"),
                            remove_genes = c(grep("^Gm",rownames(STID_obj_SS),value = T))
)
dev.off()

#> 富集图
up_gene <- SS_DEGs_microbe[[1]]$data$Overall_DEGs$gene[SS_DEGs_microbe[[1]]$data$Overall_DEGs$change == "Up"]
down_gene <- SS_DEGs_microbe[[1]]$data$Overall_DEGs$gene[SS_DEGs_microbe[[1]]$data$Overall_DEGs$change == "Down"]
length(up_gene);length(down_gene)
tmp2 <- GeneEnrichment(STID_obj = STID_obj_SS,
                       DEGs = NULL,
                       up_gene = up_gene,
                       down_gene = down_gene,
                       enrichment_method = "GO_KEGG", # GO_KEGG, 
                       go_ont = "BP",
                       return_data = T,
                       grp_nm = "JEV_microbe_DEGs_enrichment")

up_gene <- SS_DEGs_host[[1]]$data$Overall_DEGs$gene[SS_DEGs_host[[1]]$data$Overall_DEGs$change == "Up"]
down_gene <- SS_DEGs_host[[1]]$data$Overall_DEGs$gene[SS_DEGs_host[[1]]$data$Overall_DEGs$change == "Down"]
length(up_gene);length(down_gene)
tmp2 <- GeneEnrichment(STID_obj = STID_obj_SS,
                       DEGs = NULL,
                       up_gene = up_gene,
                       down_gene = down_gene,
                       enrichment_method = "GO_KEGG", # GO_KEGG, 
                       go_ont = "BP",
                       return_data = T,
                       grp_nm = "JEV_host_DEGs_enrichment")

#> cellchat
tmp_microbe <- CalSampCellComm(STID_obj = STID_obj_SS,
                               niche_key = "Niche_microbe",
                               group_by = NULL,
                               assay_id = "Spatial",
                               layer_id = "counts", # ??
                               loop_id = "LoopAllSamp",
                               col = col_lasso,
                               is_Spatial = F, # 是否进行空间cellchat，默认TRUE
                               spatial.factors = NULL, 
                               interaction.range = 250,
                               return_data = TRUE,
                               grp_nm = "JEV_microbe_Niche_CellComm",
                               dir_nm = "M3_CalSampCellComm"
)
tmp_host <- CalSampCellComm(STID_obj = STID_obj_SS,
                            niche_key = "Niche_host",
                            group_by = NULL,
                            assay_id = "Spatial",
                            layer_id = "counts", # ??
                            loop_id = "LoopAllSamp",
                            col = col_lasso,
                            is_Spatial = F, # 是否进行空间cellchat，默认TRUE
                            spatial.factors = NULL, 
                            interaction.range = 250,
                            return_data = TRUE,
                            grp_nm = "JEV_host_Niche_CellComm",
                            dir_nm = "M3_CalSampCellComm"
)


# cellchat plot
cell_types <- c(
  "Adipocytes", "Astrocytes", "Dendritic cells", "Endothelial cells", 
  "Epithelial cells", "Fibroblasts", "Macrophages", "Microglia", 
  "Monocytes", "Neurons", "NK cells", "Oligodendrocytes", "T cells"
)
CellComm_data <- readRDS("F:/Scientific_research/ShiXiaoFeng-InfectiousST/01_Project/2025/20250124_STIDtools/Analysis/outputdata/M3_CalSampCellComm/20260311_103631/CellComm_data.rds")
Plot_NicheCellComm(
  STID_obj = STID_obj_SS,
  CellComm_data = CellComm_data,
  samp_mode = "SS",
  loop_id = "LoopAllSamp",
  sources.use = cell_types,
  targets.use = cell_types,
  # signaling = c("CXCL","CCL","SAA","IL1","COMPLEMENT","TNF","ANNEXIN","ANGPTL"),
  signaling = c("CXCL","CCL","SAA","IL1","COMPLEMENT","TNF"),
  pairLR.use = NULL,
  col = col_lasso
)

CellComm_data <- readRDS("F:/Scientific_research/ShiXiaoFeng-InfectiousST/01_Project/2025/20250124_STIDtools/Analysis/outputdata/M3_CalSampCellComm/20260311_111424/CellComm_data.rds")
Plot_NicheCellComm(
  STID_obj = STID_obj_SS,
  CellComm_data = CellComm_data,
  samp_mode = "SS",
  loop_id = "LoopAllSamp",
  sources.use = cell_types,
  targets.use = cell_types,
  # signaling = c("CXCL","CCL","SAA","IL1","COMPLEMENT","TNF","ANNEXIN","ANGPTL"),
  signaling = c("CXCL","CCL","SAA","IL1","COMPLEMENT","TNF"),
  pairLR.use = NULL,
  col = col_lasso
)



### 七、Figure5 & supp Figure4 ####
# 读取STID_obj_SS_CE.rds
col_lasso <- c(
  "HsPCs" = "#E41A1C",
  "Hepatocytes" = "grey95",
  "Infla Heps" = "#4DAF4A",
  "Fibroblasts" = "#984EA3",
  "Cho/Spp1+ cells" = "#FFFF33",
  "Spp1+ MoMFs" = "#FF7F00",
  "MoKCs" = "#377EB8",
  "Neutrophils" = "#F781BF",
  "B/plasma cells" = "#A65628",
  "Others" = "#8DA0CB"
)
col_lasso2 <- c(
  "HsPCs" = "#E41A1C",
  "Hepatocytes" = "#66C2A5",
  "Infla Heps" = "#4DAF4A",
  "Fibroblasts" = "#984EA3",
  "Cho/Spp1+ cells" = "#FFFF33",
  "Spp1+ MoMFs" = "#FF7F00",
  "MoKCs" = "#377EB8",
  "Neutrophils" = "#F781BF",
  "B/plasma cells" = "#A65628",
  "Others" = "#8DA0CB"
)

### 1.CE ####
STID_obj_SS %>% print()
STID_obj_MS <- CreateMultiSampNiche(
  STID_obj = STID_obj_SS,
  multi_id = NULL, # 可以为NULL，自动生成就挺好
  loop_id = c("DPI_4_2","DPI_79_1"),
  compare_mode = "Comparative",
  niche_key = "Niche",
  description = NULL
)
STID_obj_MS <- CreateMultiSampNiche(
  STID_obj = STID_obj_MS,
  multi_id = NULL, # 可以为NULL，自动生成就挺好
  loop_id = c("DPI_0_1","DPI_4_2"),
  compare_mode = "Comparative",
  niche_key = "Niche",
  description = NULL
)
STID_obj_MS %>% print()


#> CalSampComp
GetMSNicheInfo(STID_obj_MS,
               loop_id = "LoopAllMulti",
               niche_key = "Niche")
pdf("./photo/Z_other/20260421_Figure5_result/CE/CE_multi_sample_NicheComposition_all.pdf",width = 4.5,height = 5)
CalSampComp(STID_obj = STID_obj_MS, samp_mode = "MS", 
            loop_id = "LoopAllMulti", # "LoopAllSamp","Comparative_2_4"
            samp_grp_index = T,
            niche_key = NULL,group_by = "anno",
            col = col_lasso,
            return_data = F
)
dev.off()
pdf("./photo/Z_other/20260421_Figure5_result/CE/CE_multi_sample_NicheComposition_Niche.pdf",width = 7,height = 5)
CalSampComp(STID_obj = STID_obj_MS, samp_mode = "MS", 
            loop_id = "LoopAllMulti", # "LoopAllSamp","Comparative_2_4"
            samp_grp_index = T,
            niche_key = "Niche",group_by = "anno",
            col = col_lasso,
            return_data = F
)
dev.off()


#> 病原比例
tmp <- GetSSNicheCells(STID_obj_SS,niche_key = "Niche")[[1]]
tmp <- GetMetaData(STID_obj_MS, meta_key = "M1_SpotDetect_Gene_CE_correct_after_all_gene_white",add_coord = F)[[1]]
STID_obj_MS <- AddMSNicheCells(
  STID_obj = STID_obj_MS,
  loop_id = "Comparative_2_4",
  meta_key = "M1_SpotDetect_Gene_CE_correct_after_host_gene_white",
  select_colnm = "Label_all_gene_nFeature(sum)",
  niche_key = "Niche"
)
tmp <- GetMSNicheCells(STID_obj_MS,
                       loop_id = "Comparative_2_4",
                       niche_key = "Niche")[[1]]
colnames(tmp)
pdf("./photo/Z_other/20260421_Figure5_result/CE/CE_multi_sample_NicheComposition_Pathogen.pdf",width = 4,height = 5)
CalSampComp(STID_obj = STID_obj_MS, samp_mode = "MS", 
            loop_id = "LoopAllMulti", # "LoopAllSamp","Comparative_2_4"
            samp_grp_index = T,
            meta_key = "M1_SpotDetect_Gene_CE_correct_after_host_gene_white",
            niche_key = NULL,group_by = "Label_all_gene_nFeature(sum)",
            col = rev(c("#E41A1C", "#377EB8")),
            return_data = F
)
dev.off()


#> CalSampCAI
pdf("./photo/Z_other/20260421_Figure5_result/CE/CE_multi_sample_NicheAggIndex.pdf",width = 5,height = 12)
tmp <- CalSampCAI(STID_obj = STID_obj_MS,
                  samp_mode = "MS",
                  loop_id = "LoopAllMulti",
                  samp_grp_index = T,
                  meta_key = NULL,
                  niche_key = "Niche",
                  group_by = "anno",
                  k_neighbors = 8,
                  min_agg_size = 10,
                  dist_thres = 1,
                  col = col_lasso)
dev.off()


#> CalSampDEGs
STID_obj_MS %>% print()
table(STID_obj_MS@meta.data$anno)
MS_DEGs <- CalSampDEGs(STID_obj = STID_obj_MS,
                       samp_mode = "MS",
                       loop_id = "Comparative_2_4", # 
                       samp_grp_index = T,
                       logfc_thres = 2,
                       # niche_key = "Niche", # "Niche"时DEGs很少
                       # meta_key = "coord",
                       group_by = "anno",
                       group_value = c("Neutrophils","Spp1+ MoMFs","Fibroblasts", "B/plasma cells"),
                       assay_id = "Spatial",
                       padj_thres = 0.05,
                       adjust_method = "BH",
                       col = col_lasso,
                       remove_genes = c(grep("^Gm",rownames(STID_obj_SS),value = T),
                                        grep("^EmuJ",rownames(STID_obj_SS),value = T)),
                       grp_nm = "Comparative_2_4_All"
)
tmp <- MS_DEGs[["Comparative_2_4_All"]][["data"]][["Overall_DEGs"]]
tmp <- MS_DEGs[["Comparative_2_4_All"]][["data"]][["Celltype_DEGs"]] %>% 
  filter(group == "Fibroblasts")

tmp$gene %>% head(20)
colnames(MS_DEGs$Comparative_2_4_All$data$Overall_DEGs)
DEGs_gene <- MS_DEGs$Comparative_2_4_All$data$Overall_DEGs %>% 
  filter(change %in% c("Up","Down")) %>%
  # filter(avg_log2FC > 0.25 | avg_log2FC < -0.25) %>%
  filter(p_val_adj<0.05) %>%
  pull(gene) %>% unique()
# DEGs_gene <- SS_DEGs$DPI_4_2$data$Overall_DEGs$gene %>% unique()
length(DEGs_gene)
c("Ccl3","Ccl4","Ccr5") %in% DEGs_gene
up_gene <- MS_DEGs$Comparative_2_4_All$data$Overall_DEGs %>% 
  filter(change == "Up") %>% 
  pull(gene) %>% unique()
down_gene <- MS_DEGs$Comparative_2_4_All$data$Overall_DEGs %>%
  filter(change == "Down") %>% 
  pull(gene) %>% unique()
length(up_gene);length(down_gene)

#> enrichment
tmp2 <- GeneEnrichment(STID_obj = STID_obj_SS,
                       DEGs = NULL,
                       up_gene = up_gene,
                       down_gene = down_gene,
                       enrichment_method = "GO_KEGG", # GO_KEGG, 
                       go_ont = "BP",
                       return_data = T,
                       grp_nm = "CE_Comparative_2_4_All_DEGs_enrichment")


#> Plot_NicheCellComm
# 全部样本cellchat，好像也不需要
if(0){
  tmp <- CalSampCellComm(STID_obj = STID_obj_SS,
                         niche_key = NULL,
                         group_by = NULL,
                         assay_id = "Spatial",
                         loop_id = "LoopAllSamp",
                         col = col_lasso2,
                         is_Spatial = T, # 是否进行空间cellchat，默认TRUE
                         spatial.factors = NULL, 
                         interaction.range = 250,
                         return_data = TRUE,
                         grp_nm = "CE_Comparative_2_4",dir_nm = "M3_CalSampCellComm"
  )
}
CellComm_data <- readRDS(file = "F:/Scientific_research/ShiXiaoFeng-InfectiousST/01_Project/2025/20250124_STID/Analysis/outputdata/M3_CalSampCellComm/20260414_201513/CellComm_data.rds")
pdf("./photo/Z_other/20260421_Figure5_result/CE/CE_multi_sample_NicheCellComm_all.pdf",width = 10,height = 8)
Plot_NicheCellComm(
  STID_obj = STID_obj_MS,
  CellComm_data = CellComm_data,
  samp_mode = "MS",
  loop_id = "Comparative_2_4",
  # sources.use = c("Neutrophils","Spp1+ MoMFs"),
  # targets.use = c("Neutrophils","Spp1+ MoMFs","Fibroblasts", "B/plasma cells"),
  # signaling = c("CXCL","CCL","SAA","IL1","COMPLEMENT","TNF","ANNEXIN","ANGPTL"),
  signaling = c("CXCL","CCL","SAA","SPP1","MIF","VEGF","FGF"),
  pairLR.use = NULL,
  col = col_lasso2
)
dev.off()

pdf("./photo/Z_other/20260421_Figure5_result/CE/CE_multi_sample_NicheCellComm_slect.pdf",width = 10,height = 8)
Plot_NicheCellComm(
  STID_obj = STID_obj_MS,
  CellComm_data = CellComm_data,
  samp_mode = "MS",
  loop_id = "Comparative_2_4",
  sources.use = c("Hepatocytes","Neutrophils","Spp1+ MoMFs"),
  targets.use = c("Neutrophils","Spp1+ MoMFs","Fibroblasts", "B/plasma cells"),
  # signaling = c("CXCL","CCL","SAA","IL1","COMPLEMENT","TNF","ANNEXIN","ANGPTL"),
  signaling = c("CXCL","CCL","SAA","SPP1","MIF","VEGF"),
  pairLR.use = NULL,
  col = col_lasso2
)
dev.off()


### 2.JEV ####
## 2.1 spot识别 ####
pathogen_genes <- c("NS5","C","NS3","NS1","E","Prm","NS4aAlt","NS4bAlt","NS2a","NS2b") # 只有包含至少两个 JEV 基因或三个病毒读数的 bin 才被视为 JEV 阳性
STID_obj <- GetGeneStat(STID_obj = STID_obj, features = pathogen_genes,prefix = "all_gene",func = "sum") %>% 
  AddMetaColumn(STID_obj = STID_obj,
                add_data = ., # data.frame
                meta_key = "raw", # string
                igrnore_rownm = FALSE)
pathogen_genes <- c("NS5") # 节约时间
STID_obj <- SpotDetect_Gene(STID_obj,
                            features = pathogen_genes,
                            feature_colnm = grep("all_gene",colnames(STID_obj@meta.data),value = T),
                            PosThres_prob = 0, PosThres_count = 1,
                            # col = c("#3D3576","#92D74D"),
                            # col = c("#8DBDDC","#FD9F8F"),
                            col = COLOR_DIS_CON,
                            black_bg = F,pt_size = 0.25,
                            # blur_method = "isoblur",
                            blur_method = NULL,
                            blur_n = 1,blur_sigma = 0.5, 
                            plot_method = "single",
                            grp_nm = "JEV_multisamp_microbe_gene")
STID_obj <- SpotDetect_Gene(STID_obj,
                            # features = host_genes,
                            # features = c("Ccl5","Irf7"),
                            features = host_genes,
                            feature_colnm = grep("all_gene",colnames(STID_obj@meta.data),value = T),
                            PosThres_prob = 0, PosThres_count = 4,
                            # col = c("grey95","red"),
                            col = COLOR_DIS_CON,
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#6A1B9A"),
                            # col = c("grey95","#512DA8"),
                            black_bg = F,pt_size = 0.25,
                            blur_method = NULL, blur_n = 1,blur_sigma = 0.5, 
                            plot_method = "single",
                            grp_nm = "JEV_multisamp_host_gene")

# geneset
geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/Mouse_PCD_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
geneset_df <- geneset_df[c(5)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,seed = 10,
                               PosThres_prob = 0.75, PosThres_score = 0, # 最终确定是0.75
                               pt_size = 0.25,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "JEV_multisamp_PCD")

#>
geneset_df <- read.table(file = "./inputdata/Gene_Geneset/Mouse/Geneset/GO/Mouse_GO_BP_Detect_viral_geneset.txt",
                         sep = "\t",row.names = NULL,header = T,na.strings = "")
colnames(geneset_df) <- gsub("GOBP_","",colnames(geneset_df))
geneset_df <- geneset_df[c(11)]
geneset_list <- lapply(geneset_df, function(x) na.omit(x))
names(geneset_list)
STID_obj <- SpotDetect_Geneset(STID_obj,
                               geneset_list = geneset_list,
                               score_method = "AddModuleScore", n_iter = 5, nbin = 24,
                               PosThres_prob = 0.75, PosThres_score = 0, # 最终改为0.75
                               pt_size = 0.25,
                               col = COLOR_DIS_CON,
                               black_bg = F, blur_method = NULL,
                               plot_method = "single",
                               grp_nm = "JEV_multisamp_GO_viral")



## 2.2 Niche识别 ####
# 2.2.1 JEV：microbe ####
STID_obj

#>DBSCAN
STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj, meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene", 
                                   spatial_scale_method = "region", region_detect_method = "convex", update_spots = F,
                                   ROI_size = NULL, density_thres = 1, 
                                   pos_colnm = "Label_all_gene_nFeature(sum)", 
                                   description = NULL,grp_nm = "STS_JEV_multisamp_microbe_region", dir_nm = "M2_NicheDetect_STS")
STID_obj_detect
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_STS_JEV_multisamp_microbe_region",
                           add_coord = F)[[1]]

#> plot
# SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C"  ,"#50C49FFF" ,"#3B95C4FF" ,"#B560D4FF")
# SEVEN_LIGHT <- c("#F88A7E", "#FCC093", "#FCDB9A", "#BBBFA1", "#9DB7C4", "#D1CAD4")
SEVEN_DARK <- c("#50C49FFF" ,"#FC7715FF" ,"#FCB11C"  ,"#F81B02FF" ,"#3B95C4FF" ,"#B560D4FF")
SEVEN_LIGHT <- c("#BBBFA1",  "#FCC093", "#FCDB9A","#F88A7E", "#9DB7C4", "#D1CAD4")

pdf("./photo/Z_other/20260421_Figure5_result/JEV/JEV_microbe_edge_raw.pdf",width = 60,height = 15)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "new_samp",
             datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
# pdf("./photo/Z_other/20260421_Figure5_result/JEV/JEV_microbe_region_raw.pdf",width = 15,height = 15)
# detect_meta <- detect_meta %>% 
#   group_by(new_samp) %>% 
#   arrange(All_ROI_label2) %>% 
#   mutate(new_group = dense_rank(All_ROI_label2)) %>%
#   mutate(new_group = paste0("ROI",new_group))
# Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "new_group",
#              facet_grpnm = "new_samp",
#              datatype = "discrete",
#              col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK,"grey"),con = NULL),
#              pt_size = 1.1,vmin = NULL, vmax = "p99",
#              title = NULL, subtitle = NULL,black_bg = F)
# dev.off()

#> Plot_DistLine_Exp
STID_obj_detect %>% print()
tmp <- GetMetaData(STID_obj_detect,
                   meta_key = list(c("M2_NicheDetect_STS_STS_JEV_microbe_region")))[[1]]
pdf("./photo/Z_other/20260421_Figure5_result/JEV/JEV_microbe_DistLine_Exp_host.pdf",width = 4.8,height = 6)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = c("NS5","Ccl2"), feature_colnm = "all_gene_nFeature(sum)", 
                  loop_id = "D5_1", col = c("#F81B02FF"  ,"#3B95C4FF","#F81B02FF" ) ,
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Gene_JEV_correct_before_all_gene_white",
                                                        "M2_NicheDetect_STS_STS_JEV_microbe_region")))
dev.off()
pdf("./photo/Z_other/20260421_Figure5_result/JEV/JEV_microbe_DistLine_Exp_geneset.pdf",width = 5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("RESPONSE_TO_VIRUS"),
                  loop_id = "D5_1", col = "#3B95C4FF",
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_GO_viral_white",
                                                        "M2_NicheDetect_STS_STS_JEV_microbe_region")))
dev.off()
pdf("./photo/Z_other/20260421_Figure5_result/JEV/JEV_microbe_DistLine_Exp_PCD.pdf",width = 4.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Necroptosis"), 
                  loop_id = "D5_1", col = "#3B95C4FF",
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_PCD_white",
                                                        "M2_NicheDetect_STS_STS_JEV_microbe_region")))
dev.off()


# 2.2.2 JEV：host ####
STID_obj_detect
meta_data <- GetMetaData(
  STID_obj_detect,
  meta_key = "M1_SpotDetect_Geneset_JEV_multisamp_GO_viral"
)[[1]]
STID_obj_detect <- NicheDetect_STS(STID_obj = STID_obj_detect, meta_key = "M1_SpotDetect_Geneset_JEV_multisamp_GO_viral", 
                                   spatial_scale_method = "region", region_detect_method = "convex", update_spots = T,
                                   pos_colnm = "Label_RESPONSE_TO_VIRUS",
                                   ROI_size = NULL,
                                   density_thres = 0.3,
                                   description = NULL,grp_nm = "STS_JEV_multisamp_host_region", dir_nm = "M2_NicheDetect_STS")
STID_obj_detect
detect_meta <- GetMetaData(STID_obj_detect, meta_key = "M2_NicheDetect_STS_STS_JEV_multisamp_host_region",add_coord = F)[[1]]

#> plot
SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C"  ,"#B560D4FF")
SEVEN_LIGHT <- c("#F88A7E", "#FCC093", "#FCDB9A", "#D1CAD4")
# SEVEN_DARK <- c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C"  ,"#50C49FFF" ,"#3B95C4FF" ,"#B560D4FF")
# SEVEN_LIGHT <- c("#F88A7E", "#FCC093", "#FCDB9A", "#BBBFA1", "#9DB7C4", "#D1CAD4")
pdf("./photo/Z_other/20260421_Figure5_result/JEV/JEV_host_edge_raw.pdf",width = 60,height = 15)
Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "ROI_region",
             facet_grpnm = "new_samp", datatype = "discrete",
             col = list(dis = c("grey95","#FFC4E1","#244D7F","#EB1E2C"),con = NULL),
             pt_size = 1,vmin = NULL, vmax = "p99",
             title = NULL, subtitle = NULL,black_bg = F)
dev.off()
# pdf("./photo/Z_other/20260421_Figure5_result/JEV/JEV_host_region_raw.pdf",width = 15,height = 15)
# Plot_Spatial(plot_data = detect_meta,x_colnm = "x",y_colnm = "y",group_by = "All_ROI_label2",
#              facet_grpnm = "new_samp", datatype = "discrete",
#              col = list(dis = c(SEVEN_LIGHT,SEVEN_DARK),con = NULL),
#              pt_size = 1.1,vmin = NULL, vmax = "p99",
#              title = NULL, subtitle = NULL,black_bg = F)
# dev.off()

#> Plot_DistLine_Exp
STID_obj_detect %>% print()
tmp <- GetMetaData(STID_obj_detect,
                   meta_key = list(c("M2_NicheDetect_STS_STS_JEV_host_region")))[[1]]
pdf("./photo/Z_other/20260421_Figure5_result/JEV/JEV_host_DistLine_Exp_host.pdf",width = 4.8,height = 6)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = c("NS5","Ccl2"), feature_colnm = "all_gene_nFeature(sum)", 
                  loop_id = "D5_1", col = c("#F81B02FF"  ,"#3B95C4FF","#F81B02FF" ) ,
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Gene_JEV_correct_before_all_gene_white",
                                                        "M2_NicheDetect_STS_STS_JEV_host_region")))
dev.off()
pdf("./photo/Z_other/20260421_Figure5_result/JEV/JEV_host_DistLine_Exp_geneset.pdf",width = 5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("RESPONSE_TO_VIRUS"),
                  loop_id = "D5_1", col = "#3B95C4FF",
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_GO_viral_white",
                                                        "M2_NicheDetect_STS_STS_JEV_host_region")))
dev.off()
pdf("./photo/Z_other/20260421_Figure5_result/JEV/JEV_host_DistLine_Exp_PCD.pdf",width = 4.5,height = 2.8)
Plot_DistLine_Exp(STID_obj = STID_obj_detect, features = NULL, feature_colnm = c("Necroptosis"), 
                  loop_id = "D5_1", col = "#3B95C4FF",
                  facet_grpnm = "grp",meta_key = list(c("M1_SpotDetect_Geneset_JEV_correct_before_PCD_white",
                                                        "M2_NicheDetect_STS_STS_JEV_host_region")))
dev.off()
saveRDS(STID_obj_detect,file = "./rds/STID_obj_detect_JEV_multisamp.rds")



## 2.3 SS ####
STID_obj_detect %>% print()

#> 组织颜色
col_lasso_tissue <- c('OLF' = '#FFDEAD', 'CTX_HPF' = '#8FBC8F', 'HPF' = '#A0522D', 
               'TH' = '#7FFFAA', 'HY' = '#FFC0CB', 'CNU' = '#FF8C00', 
               'MB' = '#000080', 'HB' = '#9932CC', 'CB' = '#87CEEB', 
               'FB' = '#FFFF00', 'MEN' = '#FF0000', 'CHP' = '#006400', 
               'UK' = '#D2B48C', 'MB_HY' = '#7B68EE', 'MB_HY_HB' = '#EF6FD0', 
               'HB_CB' = '#A9A9A9', 'MB_HB' = '#FF1493', 'CTX' = '#80AD80', 
               'MB_CNU' = '#C1F1A4', 'MB_CNU_HB' = '#B6A7E2', 'TH_HY' = '#F7D89A')
col_lasso_tissue2 <- c('OLF' = '#FFDEAD', 'CTX_HPF' = '#8FBC8F', 'HPF' = '#A0522D', 
                'TH' = '#7FFFAA', 'HY' = '#FFC0CB', 'CNU' = '#FF8C00', 
                'MB' = '#000080', 'HB' = '#9932CC', 'CB' = '#87CEEB', 
                'FB' = '#FFFF00', 'MEN' = '#FF0000', 'CHP' = '#006400', 
                'UK' = '#D2B48C', 'MB_HY' = '#7B68EE', 'MB_HY_HB' = '#EF6FD0', 
                'HB_CB' = '#A9A9A9', 'MB_HB' = '#FF1493', 'CTX' = '#80AD80', 
                'MB_CNU' = '#C1F1A4', 'MB_CNU_HB' = '#B6A7E2', 'TH_HY' = '#F7D89A')


#> 细胞颜色
col_lasso_cell <- c("#E41A1C", "#377EB8", '#000080',"#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#66C2A5",
               "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
names(col_lasso_cell) <- unique(STID_obj_detect@meta.data$new_cell) %>% sort()
col_lasso_cell2 <- c("#E41A1C", "#377EB8",'#000080' ,"#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#F781BF", "#66C2A5",
                "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
names(col_lasso_cell2) <- c(
  "Adipocytes", "Astrocytes", "Dendritic cells", "Endothelial cells", 
  "Epithelial cells", "Fibroblasts", "Macrophages", "Microglia", 
  "Monocytes", "Neurons", "NK cells", "Oligodendrocytes", "T cells"
)


# >>> STID_obj_SS
STID_obj_detect %>% print()
meta_data <- GetMetaData(STID_obj_detect,meta_key = list(c("M2_NicheDetect_STS_STS_JEV_multisamp_microbe_region")))[[1]]
colnames(meta_data)
STID_obj_SS <- CreateSingleSampNiche(
  STID_obj = STID_obj_detect,
  # loop_id = "DPI_4_2",
  niche_key = "Niche_microbe",
  meta_key = list(c("M2_NicheDetect_STS_STS_JEV_multisamp_microbe_region")),
  ROI_type = "ROI",
  pos_colnm = "ROI_label",
  center_colnm = "ROI_center",
  edge_colnm = "ROI_edge",
  all_label_colnm = "All_ROI_label",
  all_dist_colnm = "All_Dist2ROIcenter",
  description = NULL
)
STID_obj_SS <- CreateSingleSampNiche(
  STID_obj = STID_obj_SS,
  # loop_id = "DPI_4_2",
  niche_key = "Niche_host",
  meta_key = list(c("M2_NicheDetect_STS_STS_JEV_multisamp_host_region")),
  ROI_type = "ROI",
  pos_colnm = "ROI_label",
  center_colnm = "ROI_center",
  edge_colnm = "ROI_edge",
  all_label_colnm = "All_ROI_label",
  all_dist_colnm = "All_Dist2ROIcenter",
  description = NULL
)
STID_obj_SS <- AddSSNicheCells(
  STID_obj = STID_obj_SS,
  meta_key = "raw", # string
  select_colnm = "new_tissue", # vector
  niche_key = "Niche_microbe"
)
STID_obj_SS <- AddSSNicheCells(
  STID_obj = STID_obj_SS,
  meta_key = "raw", # string
  select_colnm = "new_tissue", # vector
  niche_key = "Niche_host"
)
STID_obj_SS %>% print()
SS_cells <- GetSSNicheCells(STID_obj_SS,niche_key = "Niche_microbe")[[1]]  
table(SS_cells$new_tissue)


## 2.4 MS ####
STID_obj_SS %>% print()
STID_obj_MS <- CreateMultiSampNiche(
  STID_obj = STID_obj_SS,
  multi_id = NULL, # 可以为NULL，自动生成就挺好
  loop_id = c("D3_1","D5_1","D7_1"),
  compare_mode = "Temporal",
  niche_key = "Niche_microbe",
  description = NULL
)
STID_obj_MS <- CreateMultiSampNiche(
  STID_obj = STID_obj_MS,
  multi_id = NULL, # 可以为NULL，自动生成就挺好
  # loop_id = c("DPI_4_2","DPI_79_1"),
  compare_mode = "Temporal",
  niche_key = "Niche_host",
  description = NULL
)
STID_obj_MS <- CreateMultiSampNiche(
  STID_obj = STID_obj_MS,
  multi_id = NULL, # 可以为NULL，自动生成就挺好
  loop_id = c("D3_1","D5_1","D7_1"),
  compare_mode = "Temporal",
  niche_key = "Niche_host",
  description = NULL
)
STID_obj_MS %>% print()


## 2.5 CalSampPathoTrack ####
STID_obj_MS <- AddMetaColumn(STID_obj_MS,meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",add_data = STID_obj_MS@meta.data["new_tissue"])
tmp <- GetMetaData(STID_obj_MS, meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",add_coord = F)[[1]]
CalSampPathoTrack(STID_obj = STID_obj_MS,
                  loop_id = "Temporal_1_2_3_4", # LoopAllMulti
                  pos_colnm = "Label_all_gene_nFeature(sum)", neg_value = "neg",
                  samp_grp_index = FALSE,
                  meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",
                  niche_key = NULL, # only support one value
                  group_by = NULL,
                  col = col_lasso_cell,
                  return_data = FALSE,
                  grp_nm = "Temporal_1_2_3_4_cell",dir_nm = "M4_CalSampPathoTrack")
CalSampPathoTrack(STID_obj = STID_obj_MS,
                  loop_id = "Temporal_1_2_3_4", # LoopAllMulti
                  pos_colnm = "Label_all_gene_nFeature(sum)", neg_value = "neg",
                  samp_grp_index = FALSE,
                  meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",
                  niche_key = NULL, # only support one value
                  group_by = "new_tissue",
                  col = col_lasso_tissue,
                  return_data = FALSE,
                  grp_nm = "Temporal_1_2_3_4_tissue",dir_nm = "M4_CalSampPathoTrack")


## 2.6 CalSampOSE ####
plan(sequential)
CalSampOSE(STID_obj = STID_obj_MS,
           loop_id = "Temporal_1_2_3_4", # LoopAllMulti
           meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",
           group_by = "new_cell",
           col = col_lasso_cell,
           only_plot = T,
           return_data = FALSE,
           grp_nm = "Temporal_1_2_3_4_tissue",dir_nm = "M4_CalSampOSE")

#> boxplot
file_paths <- list.files(path = c("./outputdata/M4_CalSampOSE/Temporal_1_2_3_4_tissue/"),
                         pattern = "entropy_results_fil",full.names = T,recursive = T)
file_nms <- basename(file_paths) %>% 
  gsub("_entropy_results_fil.txt","",.)
for (i in 1:length(file_paths)) {
  df <- read.table(file_paths[i],sep = "\t",row.names = NULL,header = T)
  df <- df %>%
    mutate(group = file_nms[i])
  if (i == 1) {
    df_merge <- df
  } else {
    df_merge <- rbind(df_merge,df)
  }
}
table(df_merge$group)
colnames(df_merge)
df_merge_bin50 <- df_merge %>% 
  mutate(group = factor(group,levels = c("CTRL_1","D3_1","D5_1","D7_1"))) %>% 
  arrange(group)
comparison <- list(c('CTRL_1','D3_1'), c('D3_1','D5_1'),c('D5_1','D7_1'))
p1 <- ggplot(df_merge_bin50,aes(x=group, y=entropy_adj,color=group))+
  # facet_grid(~group,scales = "free_x",space = "free")+
  geom_violin(width=0.8,position = position_dodge(width=0.9),
              linewidth = 0.4,
              color = "grey30",fill = NA,alpha = 0.8,
              scale = "width")+
  geom_boxplot(
    width=0.3,
    linewidth=0.25,
    outliers = F, # 是否显示异常值，只影响展示，不影响计算
    # outlier.colour = 'white',
    outlier.fill = "white",
    outlier.color = "white",
    outlier.size = 0.5,
    outlier.stroke = 0.3,
    outlier.shape=21,
    position = position_dodge(width=0.9),
    color = "grey50",fill = NA,alpha = 0.3
  )+
  geom_jitter(size=0.75,alpha=0.6,shape = 21,
              position = position_jitterdodge(jitter.width = 0.6,dodge.width = 0.9))+ # width：抖动宽度
  theme_bw()+
  # c("#4DBBD5FF","#F39B7FFF","#3C5488FF", "#E64B35FF"
  # scale_color_manual(values = color_vec)+
  # scale_fill_manual(values = color_vec)+
  # scale_x_discrete(expand = c(0.1,0,0.1,0))+
  scale_y_continuous(expand = c(0,0.05,0,0.1))+
  labs(x=NULL,y="entropy_adj",title=NULL,color="Group",fill="Group")+
  theme(
    legend.position = "none",
    # panel.grid = element_blank(),
    plot.margin = margin(0.4,0.4,0.4,0.4,'cm'),
    plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
    axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
    axis.text.x = element_text(size = 10,face = 'bold',hjust = 1,angle = 45),
    # axis.text.x = element_text(size = 10,face = 'bold',hjust = 0.5),
    axis.text.y  = element_text(size = 10,face = 'bold',hjust = 0.5),
    legend.title = element_text(size = 11,face = 'bold',hjust = 0),
    legend.text = element_text(size = 10,face = 'bold',hjust = 0),
    legend.key.size = unit(20, "pt")
  ) + 
  geom_signif(
    comparisons = comparison, 
    test = wilcox.test, # t.test/wilcox.test
    test.args = c("two.sided"), # two.sided, greater, less
    map_signif_level = function(p){
      sprintf("p = %.3g", p) # p值展示形式，F：仅p值，细节见sprintf
    },
    # map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
    size = 0.2, # 线粗
    textsize = 4, # 文本/***大小
    ##
    # y_position = 6, # 起始位置，default：自动判断最高图形高度
    margin_top = 0, # 实际标记其实位置：y_position + margin_top*y_position
    step_increase = 0.17, # 加一次比较，加总高度的比例 
    tip_length = 0.05,
    vjust = -0.2,  # 文字相对于横线上下移动文本，负值向上
    color="black"
  )
pdf("./photo/Z_other/20260421_Figure5_result/JEV/JEV_microbe_Entropy_boxplot.pdf",width = 3.5,height = 4)
print(p1)
dev.off()


## 2.7 CalSampGeneTrend ####
# all sample
tmp <- CalSampGeneTrend(STID_obj = STID_obj_MS,
                        loop_id = "Temporal_1_2_3_4", # LoopAllMulti
                        samp_grp_index = FALSE,
                        meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",
                        niche_key = NULL, # only support one value
                        group_by = NULL,
                        gene_list = NULL,
                        method = "fitting",
                        col = col_lasso_cell,
                        remove_genes = c(grep("^Gm",rownames(STID_obj_MS),value = T),
                                         grep("Rik$",rownames(STID_obj_MS),value = T)),
                        return_data = T,
                        grp_nm = "Temporal_1_2_3_4_fitting_all",
                        dir_nm = "M4_CalSampGeneTrend")

tmp2 <- CalSampGeneTrend(STID_obj = STID_obj_MS,
                         loop_id = "Temporal_1_2_3_4", # LoopAllMulti
                         samp_grp_index = FALSE,
                         meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",
                         niche_key = NULL, # only support one value
                         group_by = NULL,
                         gene_list = NULL,
                         method = "mfuzz",
                         col = col_lasso_cell,
                         remove_genes = c(grep("^Gm",rownames(STID_obj_MS),value = T),
                                          grep("Rik$",rownames(STID_obj_MS),value = T)),
                         return_data = T,
                         grp_nm = "Temporal_1_2_3_4_mfuzz_all",
                         dir_nm = "M4_CalSampGeneTrend")
up_gene <- tmp[["Temporal_1_2_3_4"]][["data"]][["up_concave_genes"]] %>% 
  dplyr::slice(1:100) %>% pull(gene)
down_gene <- tmp[["Temporal_1_2_3_4"]][["data"]][["down_concave_genes"]] %>% 
  dplyr::slice(1:100) %>% pull(gene)
GeneEnrichment(STID_obj = STID_obj_MS,
               DEGs = NULL,
               up_gene = up_gene,
               down_gene = down_gene,
               method = "GO_KEGG", # GO_KEGG, 
               go_ont = "BP",
               return_data = F,
               grp_nm = "up_down_concave")
up_gene <- tmp[["Temporal_1_2_3_4"]][["data"]][["up_convex_genes"]] %>% 
  dplyr::slice(1:100) %>% pull(gene)
down_gene <- tmp[["Temporal_1_2_3_4"]][["data"]][["down_convex_genes"]] %>% 
  dplyr::slice(1:100) %>% pull(gene)
GeneEnrichment(STID_obj = STID_obj_MS,
               DEGs = NULL,
               up_gene = up_gene,
               down_gene = down_gene,
               method = "GO_KEGG", # GO_KEGG, 
               go_ont = "BP",
               return_data = F,
               grp_nm = "up_down_convex")


#> all sample + group
tmp <- CalSampGeneTrend(STID_obj = STID_obj_MS,
                 loop_id = "Temporal_1_2_3_4", # LoopAllMulti
                 samp_grp_index = FALSE,
                 meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",
                 niche_key = NULL, # only support one value
                 group_by = "new_cell",
                 gene_list = NULL,
                 method = "fitting",
                 col = col_lasso_cell,
                 remove_genes = c(grep("^Gm",rownames(STID_obj_MS),value = T),
                                  grep("Rik$",rownames(STID_obj_MS),value = T)),
                 return_data = T,
                 grp_nm = "Temporal_1_2_3_4_fitting_celtype",
                 dir_nm = "M4_CalSampGeneTrend")

tmp2 <- CalSampGeneTrend(STID_obj = STID_obj_MS,
                 loop_id = "Temporal_1_2_3_4", # LoopAllMulti
                 samp_grp_index = FALSE,
                 meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",
                 niche_key = NULL, # only support one value
                 group_by = "new_cell",
                 gene_list = NULL,
                 method = "mfuzz",
                 col = col_lasso_cell,
                 remove_genes = c(grep("^Gm",rownames(STID_obj_MS),value = T),
                                  grep("Rik$",rownames(STID_obj_MS),value = T)),
                 return_data = T,
                 grp_nm = "Temporal_1_2_3_4_mfuzz_celtype",
                 dir_nm = "M4_CalSampGeneTrend")

#> niche + mfuzz
tmp <- CalSampGeneTrend(STID_obj = STID_obj_MS,
                        loop_id = "Temporal_2_3_4", # LoopAllMulti
                        samp_grp_index = FALSE,
                        meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",
                        niche_key = "Niche_microbe", # only support one value
                        group_by = NULL,
                        gene_list = NULL,
                        method = "mfuzz",
                        col = col_lasso_cell,
                        remove_genes = c(grep("^Gm",rownames(STID_obj_MS),value = T),
                                         grep("Rik$",rownames(STID_obj_MS),value = T)),
                        return_data = T,
                        grp_nm = "Temporal_1_2_3_4_mfuzz_microbe_niche",
                        dir_nm = "M4_CalSampGeneTrend")

tmp2 <- CalSampGeneTrend(STID_obj = STID_obj_MS,
                         loop_id = "Temporal_1_2_3_4", # LoopAllMulti
                         samp_grp_index = FALSE,
                         meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",
                         niche_key = "Niche_host", # only support one value
                         group_by = NULL,
                         gene_list = NULL,
                         method = "mfuzz",
                         col = col_lasso_cell,
                         remove_genes = c(grep("^Gm",rownames(STID_obj_MS),value = T),
                                          grep("Rik$",rownames(STID_obj_MS),value = T)),
                         return_data = T,
                         grp_nm = "Temporal_1_2_3_4_mfuzz_host_niche",
                         dir_nm = "M4_CalSampGeneTrend")

#> niche + group + mfuzz
STID_obj_MS %>% print()
tmp <- CalSampGeneTrend(STID_obj = STID_obj_MS,
                        loop_id = "Temporal_2_3_4", # LoopAllMulti
                        samp_grp_index = FALSE,
                        meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",
                        niche_key = "Niche_microbe", # only support one value
                        group_by = "new_cell",
                        gene_list = NULL,
                        method = "mfuzz",
                        col = col_lasso_cell,
                        remove_genes = c(grep("^Gm",rownames(STID_obj_MS),value = T),
                                         grep("Rik$",rownames(STID_obj_MS),value = T)),
                        return_data = T,
                        grp_nm = "Temporal_1_2_3_4_mfuzz_microbe_niche_celltype",
                        dir_nm = "M4_CalSampGeneTrend")

tmp2 <- CalSampGeneTrend(STID_obj = STID_obj_MS,
                         loop_id = "Temporal_1_2_3_4", # LoopAllMulti
                         samp_grp_index = FALSE,
                         meta_key = "M1_SpotDetect_Gene_JEV_multisamp_microbe_gene",
                         niche_key = "Niche_host", # only support one value
                         group_by = "new_cell",
                         gene_list = NULL,
                         method = "mfuzz",
                         col = col_lasso_cell,
                         remove_genes = c(grep("^Gm",rownames(STID_obj_MS),value = T),
                                          grep("Rik$",rownames(STID_obj_MS),value = T)),
                         return_data = T,
                         grp_nm = "Temporal_1_2_3_4_mfuzz_host_niche_celltype",
                         dir_nm = "M4_CalSampGeneTrend")



