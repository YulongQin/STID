

### ---------------
###
### Create: qyl
### Date: 2026-06-11 19:11:00
### Email: 524583919@qq.com
### Function:
### - Code1_Data_Process.R
###
### Update Log: 2026-06-11
###
### ---------------


if(1){
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
}else{
  library(STID)
}


#### 一、10X ####
# 如果缺少样本信息，有些会保存在GSE***_series_matrix.txt中

### 1.GSE189636: Lang ####
# 1 型 Lang 呼肠孤病毒
if(0){
  # process
  input_path <- "./inputdata/PublicData/Visium/GSE189636/GSE189636_RAW/"
  file_paths <- list.files(input_path, full.names = TRUE)
  file_nms <- list.files(input_path)
  file_nms
  rds_list <- list()
  for(i in 1:length(file_paths)){
    file_path <- file_paths[i]
    file_nm <- file_nms[i]
    print(file_nm)

    object <- Load10X_Spatial(
      data.dir = file_path, # 10X数据文件目录
      filename = "filtered_feature_bc_matrix.h5",  # 必须文件1：.h5；必须文件2：spatial文件夹，且位于同一级
      assay = "Spatial", # assay名字
      slice = "spatial", # 指定加载特定的切片,data.dir目录下至少要有tissue_lowres_image.png？
      filter.matrix = TRUE # 过滤掉低质量的细胞
    )
    object@meta.data$sample <- file_nm
    rds_list[[file_nm]] <- object
  }
  stRNA <- merge(rds_list[[1]],rds_list[2:4],add.cell.ids = file_nms)
  stRNA <- JoinLayers(stRNA)
  row_nm <- rownames(stRNA)
  write.table(row_nm,file = "./inputdata/PublicData/Visium/GSE189636/process/gene_list.txt",
              sep = "\t",quote = F,row.names = F,col.names = F)
  rownames(stRNA) <- gsub("GRCm38-","",rownames(stRNA))
  rownames(stRNA) <- gsub("ReoT1L-","",rownames(stRNA))
  dir.create("./inputdata/PublicData/Visium/GSE189636/process/",showWarnings = F,recursive = T)
  saveRDS(stRNA,file = "./inputdata/PublicData/Visium/GSE189636/process/stRNA.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/Visium/GSE189636/process/stRNA_Lang.rds")
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  stRNA <- NormalizeData(stRNA)
  samps_nm <- c("Mock_D4","Mock_D7","T1L_D4","T1L_D7")
  names(samps_nm) <- c("Mock_Heart_D4PI","Mock_Heart_D7PI","T1L_Heart_D4PI","T1L_Heart_D7PI")
  stRNA@meta.data$sample <- samps_nm[stRNA@meta.data$sample]
  meta_data <- stRNA@meta.data
  table(meta_data$sample)
  STID_obj <- as.STID(stRNA,samp_colnm = "sample",samp_grp_colnm = "sample",celltype_colnm = "sample",
                    host_org = "mouse",pathogen_grp = "virus",pathogen_org = "Lang",
                    data_format = "hex_grid",data_platform = "Visium")
  # pathogen_genes <- grep("T1LReo",rownames(STID_obj),value = T)
  # length(pathogen_genes)
  pathogen_genes <- c("T1LReoS1","T1LReoS2","T1LReoS3",
                     "T1LReoS4","T1LReoM1","T1LReoM2",
                     "T1LReoM3","T1LReoL1","T1LReoL2",
                     "T1LReoL3") # 这些基因都没有表达？？？？？，只能用基因集打分
  host_genes <- c("Bst2", "Isg15", "B2m", "Ifi27l2a", "Iigp1", "H2-K1",
                  "Irf7", "H2-D1", "Ly6e", "Psmb8", "Ifit3", "Ifit1",
                  "Gbp7", "Ly6a", "Rtp4", "Lgals3bp", "H2-Q7", "Xaf1",
                  "H2-T23", "Oasl2", "Rsad2", "Rnf213", "Gbp3", "Psmb9",
                  "Psmb10", "Lgals9", "Ifi203", "Irgm1", "Ifit2")
}


### 2.GSE275201: VEE ####
# 委内瑞拉马脑炎病毒VEE
if(0){
  input_path <- "./inputdata/PublicData/Visium/GSE275201/GSE275201_RAW/"
  file_paths <- list.files(input_path, full.names = TRUE)
  file_nms <- list.files(input_path)
  file_nms
  rds_list <- list()
  for(i in 1:length(file_paths)){
    file_path <- file_paths[i]
    file_nm <- file_nms[i]
    print(file_nm)

    spe2 <- Read10X(paste0(file_path,"/filtered_feature_bc_matrix/"))
    image2 <- Read10X_Image(image.dir = paste0(file_path,"/spatial/"),
                            # image.name = "tissue_hires_image.png", # 默认是tissue_lowres_image.png
                            filter.matrix = TRUE)
    spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")
    image2 <- image2[Cells(x = spe2)]
    DefaultAssay(spe2) <- "Spatial"
    DefaultAssay(image2) <- "Spatial"
    spe2 <- subset(spe2,cells = rownames(image2@coordinates))
    spe2[["slice1"]] <- image2
    for (i in colnames((spe2@images$slice1@coordinates))) {
      spe2@images$slice1@coordinates[[i]] <- as.integer(spe2@images$slice1@coordinates[[i]])
    }
    spe2@meta.data$sample <- file_nm
    rds_list[[file_nm]] <- spe2
  }
  dir.create("./inputdata/PublicData/Visium/GSE275201/process/",showWarnings = F,recursive = T)
  stRNA <- merge(rds_list[[1]],rds_list[2],add.cell.ids = file_nms)
  stRNA <- JoinLayers(stRNA)
  row_nm <- rownames(stRNA)
  write.table(row_nm,file = "./inputdata/PublicData/Visium/GSE275201/process/gene_list.txt",
              sep = "\t",quote = F,row.names = F,col.names = F)
  # rownames(stRNA) <- gsub("GRCm38-","",rownames(stRNA))
  # rownames(stRNA) <- gsub("VEEV-","",rownames(stRNA))
  saveRDS(stRNA,file = "./inputdata/PublicData/Visium/GSE275201/process/stRNA.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/Visium/GSE275201/process/stRNA_VEE.rds")
  stRNS <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  stRNA <- NormalizeData(stRNA)
  meta_data <- stRNA@meta.data
  table(meta_data$sample)
  STID_obj <- as.STID(stRNA,samp_colnm = "sample",samp_grp_colnm ="sample",celltype_colnm = "sample",
                    host_org = "mouse",pathogen_grp = "virus",pathogen_org = "VEEV",
                    data_format = "hex_grid",data_platform = "Visium")
  pathogen_genes <- NULL
  host_genes <- c("Cd72","Ly6c2","Plac8","Ccl2","Ccl7","Ccl4","Cd3e","Nkg7","Ifng",
                  "Gzmb","Prf1","Klrc1")
}


### 3.GSE243367: HBV × ####
# 人-乙肝病毒：缺少病原信息，not useful
if(0){
  input_path <- "./inputdata/PublicData/Visium/GSE243367/GSE243367_RAW/"
  file_paths <- list.files(input_path, full.names = TRUE)
  file_nms <- list.files(input_path)
  file_nms
  rds_list <- list()
  for(i in 1:length(file_paths)){
    file_path <- file_paths[i]
    file_nm <- file_nms[i]
    print(file_nm)

    spe2 <- Read10X(paste0(file_path,"/filtered_feature_bc_matrix/"))
    image2 <- Read10X_Image(image.dir = paste0(file_path,"/spatial/"),
                            # image.name = "tissue_hires_image.png", # 默认是tissue_lowres_image.png
                            filter.matrix = TRUE)
    spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")
    image2 <- image2[Cells(x = spe2)]
    DefaultAssay(spe2) <- "Spatial"
    DefaultAssay(image2) <- "Spatial"
    spe2 <- subset(spe2,cells = rownames(image2@coordinates))
    spe2[["slice1"]] <- image2
    for (i in colnames((spe2@images$slice1@coordinates))) {
      spe2@images$slice1@coordinates[[i]] <- as.integer(spe2@images$slice1@coordinates[[i]])
    }
    spe2@meta.data$sample <- file_nm
    rds_list[[file_nm]] <- spe2
  }
  dir.create("./inputdata/PublicData/Visium/GSE243367/process/",showWarnings = F,recursive = T)
  stRNA <- merge(rds_list[[1]],rds_list[2:9],add.cell.ids = file_nms)
  stRNA <- JoinLayers(stRNA)
  row_nm <- rownames(stRNA)
  write.table(row_nm,file = "./inputdata/PublicData/Visium/GSE243367/process/gene_list.txt",
              sep = "\t",quote = F,row.names = F,col.names = F)
  saveRDS(stRNA,file = "./inputdata/PublicData/Visium/GSE243367/process/stRNA.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/Visium/GSE243367/process/stRNA_HBV.rds")
  stRNA <- NormalizeData(stRNA)
  meta_data <- stRNA@meta.data
  STID_obj <- as.STID(stRNA,samp_colnm = "sample",celltype_colnm = "sample",
                    host_org = "human",pathogen_grp = "virus",pathogen_org = "HBV",
                    data_format = "hex_grid",data_platform = "Visium")
  pathogen_genes <- NULL
  host_genes <- NULL
}

### 4.GSE306575: HBV ####
# 人-乙肝病毒
if(0){
  pre_processed <- readRDS(file = "./inputdata/PublicData/Visium/GSE306575/GSE306575_chb_ST_pre-processed_Seurat.RDS")
  processed <- readRDS(file = "./inputdata/PublicData/Visium/GSE306575/GSE306575_chb_ST_processed_Seurat.RDS")
  spatial_coords <- read.table(file = "./inputdata/PublicData/Visium/GSE306575/spatial_coords.txt",
                               sep = "\t",header = T,row.names = 1)
  tmp <- processed@meta.data
  tmp2 <- tmp %>% rownames_to_column("row_nm") %>% mutate(row_nm,.before = 1) %>%  group_by(orig.ident) %>% dplyr::slice(1)
  names(pre_processed) %in% tmp2$orig.ident 
  for (i in 1:length(pre_processed)) {
    i_nm <- names(pre_processed)[i]
    print(i_nm)
    sub_str <- tmp2$row_nm[tmp2$orig.ident == i_nm]
    sub_str <- gsub(".*-","",sub_str)
    print(sub_str)
    colnames(pre_processed[[i]]) <- gsub("-.*",paste0("-",sub_str),colnames(pre_processed[[i]]))
  }
  seurat_merge_raw <- merge(pre_processed[[1]],pre_processed[2:length(pre_processed)])
  colnames(seurat_merge_raw)
  seurat_merge_raw <- JoinLayers(seurat_merge_raw)
  seurat_merge_raw <- NormalizeData(seurat_merge_raw)
  seurat_merge_raw@meta.data$category <- processed@meta.data$category[match(rownames(seurat_merge_raw@meta.data),
                                                                  rownames(processed@meta.data))]
  table(seurat_merge_raw@meta.data$category,useNA = "always")
  seurat_merge_raw@meta.data$category <- if_else(seurat_merge_raw@meta.data$category %in% c("cccDNA","iDNA"),seurat_merge_raw@meta.data$category,
                                                 "uninf",missing = "uninf")
  seurat_merge_raw <- subset(seurat_merge_raw,cells = rownames(spatial_coords))
  SpatialDimPlot(seurat_merge_raw,group.by = "category",pt.size.factor = 1.5,
                 ncol = 2)
  
  #> 下面这种少了许多点，不太好
  seurat_merge <- subset(seurat_merge_raw,cells = colnames(processed))
  table(rownames(seurat_merge@meta.data) == rownames(processed@meta.data))
  table(rownames(seurat_merge@meta.data) %in% rownames(processed@meta.data))
  seurat_merge@meta.data <- tmp[rownames(seurat_merge@meta.data),]
  table(rownames(seurat_merge@meta.data) %in% rownames(processed@meta.data))
  seurat_merge <- UpdateSeuratObject(seurat_merge)
  seurat_merge <- NormalizeData(seurat_merge)
  
  SpatialDimPlot(seurat_merge,group.by = "category",pt.size.factor = 1.5,
                 ncol = 2)

  saveRDS(seurat_merge_raw,file = "./inputdata/PublicData/Visium/GSE306575/process/stRNA_HBV.rds")

  }else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/Visium/GSE306575/process/stRNA_HBV.rds")
  # stRNA <- NormalizeData(stRNA)
  stRNA@meta.data$sample <- paste0("CHB_",gsub(".*-","",colnames(stRNA)))
  meta_data <- stRNA@meta.data
  pathogen_genes <- NULL
  host_genes <- c("IFIT2","IFI27","CXCL9","CXCL10")
  STID_obj <- as.STID(stRNA,samp_colnm = "sample", samp_grp_colnm = "sample", celltype_colnm = "category",
                      host_org = "human",pathogen_grp = "virus",pathogen_org = "HBV",
                      pathogen_genes = pathogen_genes,
                      # binsize = 50,
                      data_format = "hex_grid",data_platform = "Visium")

}

# HBV_S_CDS	比对或分配至 HBV S coding sequence 的 counts，主要对应表面抗原 HBsAg 的编码区域。
# HBV_Core_CDS	比对至 HBV precore/core coding sequence 的 counts，与 core/precore 转录本有关。
# HBV_X_CDS	比对至 HBV X coding sequence 的 counts，与 HBx 转录区域有关。
# HBV_DR1down	比对至 HBV DR1 下游自定义区间的 counts。它不是一个标准 HBV 基因名称，而更可能是作者在自定义 HBV 注释中定义的基因组区间。

#> table(tmp$category)
# Healthy   uninf  cccDNA    iDNA 
# 11355    6382    6425   25324
# 做病原阳性与阴性的二分类时，可以把 cccDNA 和 iDNA 合并
# Healthy	11,355	22.95%	来自健康对照肝组织
# uninf	6,382	12.90%	来自 CHB 组织，但该 spot 未检测到符合分类要求的 HBV 信号
# cccDNA	6,425	12.98%	HBV 信号模式被判定更符合 cccDNA 转录来源
# iDNA	25,324	51.17%	HBV 信号模式被判定更符合 整合型 HBV DNA（integrated HBV DNA）转录来源


### 5.GSE190225: kp ####
# 肺炎克雷伯菌
if(0){
  load("./inputdata/PublicData/Visium/GSE190225/merge_lung_spatial.RData") # merge.lung
  stRNA <- merge.lung
  DefaultAssay(stRNA)
  meta_data <- stRNA@meta.data
  table(meta_data$orig.ident) # 四张芯片
  DefaultAssay(stRNA) <- "Spatial"
  stRNA[["Spatial"]] <- as(object = stRNA[["Spatial"]], Class = "Assay5") # V4 -> V5
  stRNA[["SCT"]] <- as(object = stRNA[["SCT"]], Class = "Assay5") # V4 -> V5
  stRNA[["predictions"]] <- as(object = stRNA[["predictions"]], Class = "Assay5") # V4 -> V5
  stRNA[["RCTD"]] <- as(object = stRNA[["RCTD"]], Class = "Assay5") # V4 -> V5
  table(colnames(stRNA[["RCTD"]]) == colnames(stRNA))
  stRNA@meta.data$cell_type <- rownames(stRNA[["RCTD"]])[apply(stRNA[["RCTD"]]$data, 2, which.max)]
  dir.create("./inputdata/PublicData/Visium/GSE190225/process/",showWarnings = F,recursive = T)
  saveRDS(stRNA,file = "./inputdata/PublicData/Visium/GSE190225/process/stRNA.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/Visium/GSE190225/process/stRNA_Kp.rds")
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  stRNA <- NormalizeData(stRNA)
  meta_data <- stRNA@meta.data
  table(meta_data$orig.ident)
  STID_obj <- as.STID(stRNA,samp_colnm = "orig.ident",samp_grp_colnm = "orig.ident",celltype_colnm = "cell_type",
                    host_org = "mouse",pathogen_grp = "bacteria",pathogen_org = "Kpneumoniae",
                    data_format = "hex_grid",data_platform = "Visium")
  pathogen_genes <- NULL
  host_genes <- c("Ccl20","Scgb1a1")
}


### 6.GSE276083: MTB × ####
# 结核杆菌，数据不太对劲，不可用
if(0){
  library(data.table)
  # expr_mat <- read.table("F:/Scientific_research/ShiXiaoFeng-InfectiousST/01_Project/2025/20250124_STID/Analysis/inputdata/PublicData/Visium/GSE276083/CosMx_raw_counts.csv",
  #                        header = TRUE, row.names = 1,sep = ",", check.names = F)
  expr_mat <- fread(file = "F:/Scientific_research/ShiXiaoFeng-InfectiousST/01_Project/2025/20250124_STID/Analysis/inputdata/PublicData/Visium/GSE276083/CosMx_raw_counts.csv",
                    header = TRUE,sep = ",", check.names = F)
  expr_mat <- expr_mat %>% as.data.frame() %>% 
    column_to_rownames("V1")
  expr_mat <- as.matrix(expr_mat)
  stRNA <- CreateSeuratObject(counts = expr_mat, assay = "Spatial")
  
  DefaultAssay(stRNA) <- "Spatial"
  meta_data <- stRNA@meta.data
  
  expr_mat_data <- fread(file = "F:/Scientific_research/ShiXiaoFeng-InfectiousST/01_Project/2025/20250124_STID/Analysis/inputdata/PublicData/Visium/GSE276083/CosMx_normalized_matrix.txt",
                    header = TRUE,sep = "\t", check.names = F)
  expr_mat_data <- expr_mat_data %>% as.data.frame() %>% 
    column_to_rownames("V1")
  expr_mat_data <- as.matrix(expr_mat_data)
  norm_sparse <- as(expr_mat_data, "dgCMatrix")
  stRNA@assays[["Spatial"]]["data"] <- norm_sparse
  
  # > generate object @images$slice1
  # 这个版本目前已经不能适配visiumV2了
  generate_stereo_spatial <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
    if (filter.matrix) {
      tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
    }
    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * 
      scale.factors$tissue_lowres_scalef
    spot.radius <- unnormalized.radius / max(dim(x = image))
    return(new(Class = 'VisiumV1',
               image = image,
               scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                            fiducial = scale.factors$fiducial_diameter_fullres,
                                            hires = scale.factors$tissue_hires_scalef,
                                            lowres = scale.factors$tissue_lowres_scalef),
               coordinates = tissue.positions,
               spot.radius = spot.radius))
  }
  
  # 从tissue_positions_list.csv中获取
  tissue_positions_list <- read.csv("F:/Scientific_research/ShiXiaoFeng-InfectiousST/01_Project/2025/20250124_STID/Analysis/inputdata/PublicData/Visium/GSE276083/CosMx_spatial_coordinates.csv",
                                    header = T,row.names = 1, stringsAsFactors = FALSE)
  tissue_positions_list <- tissue_positions_list %>% 
    dplyr::rename(row = Height, col = Width) %>%
    mutate(imagerow = row, imagecol = col) %>%
    mutate(tissue = 1,.before = 1)
  tissue_lowres_image <- matrix(1, max(tissue_positions_list$imagerow), max(tissue_positions_list$imagecol))

  scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1, 
                                   tissue_hires_scalef = 1, 
                                   tissue_lowres_scalef = 1))
  stereo_spatial <- generate_stereo_spatial(image = tissue_lowres_image, # 结果等价于Read10X_Image
                                            scale.factors = fromJSON(scalefactors_json), # 啥作用？
                                            tissue.positions = tissue_positions_list)
  DefaultAssay(stereo_spatial) <- "Spatial"
  stereo_spatial
  
  
  # > 整合到seurat对象中
  stereo_spatial <- stereo_spatial[Cells(x = stRNA)]
  stRNA <- subset(stRNA,cells = rownames(stereo_spatial@coordinates))
  stRNA[["slice1"]] <- stereo_spatial
  if(0){
    # 没有报错，无需下面转化
    for (i in colnames((stRNA@images$slice1@coordinates))) {
      stRNA@images$slice1@coordinates[[i]] <- as.integer(stRNA@images$slice1@coordinates[[i]])
    }
  }
  stRNA
  
  meta_data <- read.table("./inputdata/PublicData/Visium/GSE276083/CosMx_annotations.csv",sep = ",",
                          header = T,row.names = 1, stringsAsFactors = FALSE)
  stRNA@meta.data <- meta_data
  # stRNA <- NormalizeData(stRNA) # CosMX不可
  saveRDS(stRNA,file = "./inputdata/PublicData/Visium/GSE276083/process/stRNA_MTB.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/Visium/GSE276083/process/stRNA_MTB.rds")
  Plot_Spatial(stRNA,x_colnm = "Width",y_colnm = "Height",group_by = "nb_clus",facet_grpnm = "fov",ncol = 6,
               pt_size = 0.5,
               col = list(dis = COLOR_LIST$PALETTE_WHITE_BG,con = NULL))
  # stRNA <- NormalizeData(stRNA)
  meta_data <- stRNA@meta.data
  pathogen_genes <- NULL
  host_genes <- c("SPP1", "FN1", "MARCO", "CD44", "CD68", "CD163", "CD14", "LYZ")
  STID_obj <- as.STID(seurat_merge_raw,samp_colnm = "orig.ident", samp_grp_colnm = "orig.ident", celltype_colnm = "orig.ident",
                      host_org = "human",pathogen_grp = "virus",pathogen_org = "HBV",
                      pathogen_genes = pathogen_genes,
                      # binsize = 50,
                      data_format = "hex_grid",data_platform = "Visium")
}



### 7.GSE268068: Ma ####
# 疟原虫，可以读取原始GEO数据，但是优先使用处理好的数据
if(0){
  input_path <- "./inputdata/PublicData/Visium/GSE268068/GSE268068_RAW/"
  file_paths <- list.files(input_path, full.names = TRUE)
  file_nms <- list.files(input_path)
  file_nms
  rds_list <- list()
  for(i in 1:length(file_paths)){
    file_path <- file_paths[i]
    file_nm <- file_nms[i]
    print(file_nm)

    object <- Load10X_Spatial(
      data.dir = file_path, # 10X数据文件目录
      filename = "filtered_feature_bc_matrix.h5",  # 必须文件1：.h5；必须文件2：spatial文件夹，且位于同一级
      assay = "Spatial", # assay名字
      slice = "spatial", # 指定加载特定的切片,data.dir目录下至少要有tissue_lowres_image.png？
      filter.matrix = TRUE # 过滤掉低质量的细胞
    )
    object@meta.data$sample <- file_nm
    rds_list[[file_nm]] <- object
  }
  stRNA <- merge(rds_list[[1]],rds_list[2:3],add.cell.ids = file_nms)
  stRNA <- JoinLayers(stRNA)
  dir.create("./inputdata/PublicData/Visium/GSE268068/process/",showWarnings = F,recursive = T)
  row_nm <- rownames(stRNA)
  write.table(row_nm,file = "./inputdata/PublicData/Visium/GSE268068/process/gene_list.txt",
              sep = "\t",quote = F,row.names = F,col.names = F)
  # rownames(stRNA) <- gsub("GRCm38-","",rownames(stRNA))
  # rownames(stRNA) <- gsub("ReoT1L-","",rownames(stRNA))

  saveRDS(stRNA,file = "./inputdata/PublicData/Visium/GSE268068/process/stRNA.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/Visium/GSE268068/process/stRNA_Ma.rds")
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  stRNA <- NormalizeData(stRNA)
  meta_data <- stRNA@meta.data
  table(meta_data$sample)
  STID_obj <- as.STID(stRNA,samp_colnm = "sample",samp_grp_colnm = "sample",celltype_colnm = "sample",
                    host_org = "mouse",pathogen_grp = "parasite",pathogen_org = "plasmodium",
                    data_format = "hex_grid",data_platform = "Visium")
  pathogen_genes <- grep("PBANKA",rownames(STID_obj),value = T)
  pathogen_genes2 <- c("HSP70-pb", "HSP90-pb", "LISP2-pb", "H4-pb", "H2B-pb")
  grep("HSP70",rownames(STID_obj),value = T)

  length(pathogen_genes)
  # 干扰素刺激基因（ISG）
  ISGs <- c("Lgals3bp", "Trim30a", "Ly6a", "Ifi44", "Ifit3", "Ifit3b",
            "Rsad2", "Cmpk2", "Bst2", "Irf7", "Ifi27l2a", "Zbp1",
            "Parp14", "Herc6", "Usp18", "Shisa5", "Ifih1", "Fabp5", "Oasl1")
  host_genes <- ISGs
  # 炎症热点”（IHS）
  # IHSs <- c()
}

# 处理后的数据也没有注释结果
if(0){
  MP <- readRDS(file = "./inputdata/PublicData/Visium/GSE268068/se_visium.rds")
  meta_data <- MP@meta.data
  table(meta_data$sample_id)
  meta_data$nCount_pb %>% summary() # pb就是寄生虫的检出
  table(meta_data$timepoint)
  MP2 <- readRDS(file = "./inputdata/PublicData/Visium/GSE268068/final_merged_curated_annotations_270623.RDS")
  MP3 <- readRDS(file = "./inputdata/PublicData/Visium/GSE268068/STUtility_mus_pb_ST.RDS")
}


### 8.GSE200642: Tbb ####
# 布氏锥虫
if(0){
  trypanosome_Naive <- readRDS(file = "./inputdata/PublicData/Visium/GSE200642/Spatial_Naive.RDS")
  trypanosome_25dpi <- readRDS(file = "./inputdata/PublicData/Visium/GSE200642/Spatial_25dpi.RDS")
  trypanosome_45dpi <- readRDS(file = "./inputdata/PublicData/Visium/GSE200642/Spatial_45dpi.RDS")
  stRNA <- merge(trypanosome_Naive,list(trypanosome_25dpi,trypanosome_45dpi),
                     add.cell.ids = c("Naive","Infd25","Infd45"))
  stRNA@assays$Spatial <- as(object = stRNA@assays$Spatial, Class = "Assay5") # V4 -> V5
  stRNA@assays$SCT <- as(object = stRNA@assays$SCT, Class = "Assay5") # V4 -> V5
  table(stRNA@meta.data$seurat_clusters) # 不是这个注释
  dir.create("./inputdata/PublicData/Visium/GSE200642/process/",showWarnings = F,recursive = T)
  row_nm <- rownames(stRNA)
  write.table(row_nm,file = "./inputdata/PublicData/Visium/GSE200642/process/gene_list.txt",
              sep = "\t",quote = F,row.names = F,col.names = F)
  saveRDS(stRNA,file = "./inputdata/PublicData/Visium/GSE200642/process/stRNA.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/Visium/GSE200642/process/stRNA_Tbb.rds")
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  meta_data <- stRNA@meta.data
  STID_obj <- as.STID(stRNA,samp_colnm = "group",celltype_colnm = NULL, samp_grp_colnm = "group",
                    host_org = "mouse",pathogen_grp = "parasite",pathogen_org = "trypanosome",
                    data_format = "hex_grid",data_platform = "Visium")
  # pathogen_genes <- grep("^Tb",rownames(STID_obj),value = T) # 基因异常，没法使用all基因
  pathogen_genes <- c("Tb927.7.5940","Tb927.6.4280") # GADPH, PAD2
  length(pathogen_genes)
  host_genes <- c("Ttr","Cd138","Il10","Il10ra","Aif1","Chil3","Arg1")
}


### 9.GSE290027: Ld ####
# 杜氏利什曼虫
if(0){
  stRNA <- readRDS(file = "./inputdata/PublicData/Visium/GSE290027/with47_48_integrated_nd2117slcl_dims15_res0.9_prepsct.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/Visium/GSE290027/process/stRNA_Ld.rds")
  # stRNA <- NormalizeData(stRNA)
  meta_data <- stRNA@meta.data
  pathogen_genes <- NULL
  host_genes <- c("IDO1", "CD274", "IL32", "CCL18", "CXCL9", "IL24", "IL1B", "IFNGR2", "TNFRSF6B", "LYZ", "GBP1", "IFI30", "FTH1", "TNF", "IFNG")
  STID_obj <- as.STID(stRNA,samp_colnm = "orig.ident", samp_grp_colnm = "orig.ident", celltype_colnm = "morpho",
                      host_org = "human",pathogen_grp = "parasite",pathogen_org = "Ld",
                      pathogen_genes = pathogen_genes,
                      # binsize = 50,
                      data_format = "hex_grid",data_platform = "Visium")
}


### 二、Stereo-Seq ####
### 1.AE ####
if(0){
  # process
  AE <- readRDS(file = "F:/Scientific_research/ShiXiaoFeng-InfectiousST/01_Project/2025/20250124_STIDtools/Analysis/inputdata/PublicData/StereoSeq/AE/hefan/stRNA.rds")
  DefaultAssay(AE) <- "Spatial"
  AE@assays$SCT <- NULL
  table(AE@meta.data$batch)
  stRNA <- subset(AE,subset = batch %in% c("DPI_0_1","DPI_4_2","DPI_8_2","DPI_15_2","DPI_15_3","DPI_79_1"))
  dir.create("./inputdata/PublicData/StereoSeq/AE/process/",showWarnings = F,recursive = T)
  saveRDS(stRNA,file = "./inputdata/PublicData/StereoSeq/AE/process/stRNA.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/StereoSeq/AE/process/stRNA_AE.rds")
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  stRNA <- subset(stRNA,subset = batch %in% c("DPI_0_1","DPI_4_2","DPI_15_3","DPI_79_1")) # 有效样本
  # saveRDS(stRNA,file = "./inputdata/PublicData/StereoSeq/AE/process/stRNA_AE.rds")
  meta_data <- stRNA@meta.data
  table(meta_data$batch)
  # pathogen_genes
  pathogen_genes <- grep("EmuJ-",rownames(stRNA),value = T) # 打分用的估计也是所有的基因
  length(pathogen_genes)
  STID_obj <- as.STID(stRNA,samp_colnm = "batch", samp_grp_colnm = "group", celltype_colnm = "anno",
                    host_org = "mouse",pathogen_grp = "parasite",pathogen_org = "EmuJ",
                    pathogen_genes = pathogen_genes,binsize = 100, # 50
                    data_format = "square_grid",data_platform = "StereoSeq")
  # host_genes <- c("Il1b","Spp1","Ccl3","Ccr1","Casp4","Mpo")
  host_genes <- c("Il1b","Spp1")
  STID_obj %>% print()
}

#> 细胞类型介绍
# [1] "HsPCs"           "Hepatocytes"     "Infla Heps"      "Fibroblasts"     "Cho/Spp1+ cells" "Spp1+ MoMFs"
# [7] "MoKCs"           "Neutrophils"     "B/plasma cells"  "Others"
# MoKCs	Monocyte-derived Kupffer Cells（单核来源库普弗细胞）
# Spp1+ MoMFs	Spp1+ Monocyte-derived Macrophages（单核来源巨噬细胞）
# Cho/Spp1+ cells	Spp1+ Cholangiocytes（表达Spp1的胆管上皮细胞）
# Fibroblasts	成纤维细胞（在肝中主要指活化的肝星状细胞）
# Infla Heps	Inflammatory Hepatocytes（炎症性肝细胞）

# 中性粒 (S100a8/a9+Il1b) → 释放 IL1β 激活下游；
# 调控单核分化：走向MoKC（清除病原体） & Spp1+ MoMF（修复 / 纤维化） 两大分支


### 2.AE_raw ####
if(0){
  # process
  DPI_0_1_raw <- readRDS(file = "./inputdata/PublicData/StereoSeq/AE/DPI_0_1_raw.rds")
  DPI_0_1_raw@meta.data$orig.ident <- "DPI_0_1"
  DPI_4_2_raw <- readRDS(file = "./inputdata/PublicData/StereoSeq/AE/DPI_4_2_raw.rds")
  DPI_4_2_raw@meta.data$orig.ident <- "DPI_4_2"
  DPI_15_3_raw <- readRDS(file = "./inputdata/PublicData/StereoSeq/AE/DPI_15_3_raw.rds")
  DPI_15_3_raw@meta.data$orig.ident <- "DPI_15_3"
  seurat_merge <- merge(DPI_0_1_raw,list(DPI_4_2_raw,DPI_15_3_raw),
                     add.cell.ids = c("DPI_0_1","DPI_4_2","DPI_15_3"))
  seurat_merge <- JoinLayers(seurat_merge)
  seurat_merge <- NormalizeData(seurat_merge)
  meta_data <- seurat_merge@meta.data
  saveRDS(seurat_merge,file = "./inputdata/PublicData/StereoSeq/AE/process/stRNA_AE_raw.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/StereoSeq/AE/process/stRNA_AE_raw.rds")
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  meta_data <- stRNA@meta.data
  table(meta_data$orig.ident)
  # pathogen_genes
  rownames(stRNA) <- gsub("EmuJ_","EmuJ-",rownames(stRNA))
  pathogen_genes <- grep("EmuJ-",rownames(stRNA),value = T) # 打分用的估计也是所有的基因
  length(pathogen_genes)
  STID_obj <- as.STID(stRNA,samp_colnm = "orig.ident", samp_grp_colnm = "orig.ident", celltype_colnm = "orig.ident",
                      host_org = "mouse",pathogen_grp = "parasite",pathogen_org = "EmuJ",
                      pathogen_genes = pathogen_genes,binsize = 100, # 50
                      data_format = "square_grid",data_platform = "StereoSeq")
  # host_genes <- c("Il1b","Spp1","Ccl3","Ccr1","Casp4","Mpo")
  host_genes <- c("Il1b","Spp1")
  STID_obj %>% print()
  rm(stRNA);gc()
}

### 3.JEV ####
## 3.1 bin35 ####
# JEV在云平台转换h5ad -> rds
if(0){
  # process
  CTRL_1 <- readRDS(file = "./inputdata/PublicData/StereoSeq/JEV/THU48_new_annotated.rds") # bin35
  D3_1 <- readRDS(file = "./inputdata/PublicData/StereoSeq/JEV/THU34_new_annotated.rds")
  D5_1 <- readRDS(file = "./inputdata/PublicData/StereoSeq/JEV/THU50_new_annotated.rds")
  # D5_2 <- readRDS(file = "./inputdata/PublicData/StereoSeq/JEV/THU52_new_annotated.rds")
  # D5_3 <- readRDS(file = "./inputdata/PublicData/StereoSeq/JEV/THU54_new_annotated.rds")
  D7_1 <- readRDS(file = "./inputdata/PublicData/StereoSeq/JEV/THU66_new_annotated.rds")
  meta_data <- D3_1@meta.data
  table(meta_data$brain_region)
  # SpatialPlot(D3_1,features = c("Cxcl10"),pt.size.factor = 0.5) + # 用于判断形状？
  #   ggtitle("D5_1_JEV")
  ncol(CTRL_1)
  ncol(D3_1)
  ncol(D5_1)
  ncol(D7_1)
  if(0){
    # 这是旧的注释，不用了，用陈旗最新给的
    CTRL_1_meta <- read.table(file = "./inputdata/PublicData/StereoSeq/JEV/annotation/THU48_singler_annotation_metadata-self.csv",
                              sep = ",",row.names = 1,header = T)
    D3_1_meta <- read.table(file = "./inputdata/PublicData/StereoSeq/JEV/annotation/THU34_singler_annotation_metadata-self.csv",
                            sep = ",",row.names = 1,header = T)
    D5_1_meta <- read.table(file = "./inputdata/PublicData/StereoSeq/JEV/annotation/THU50_singler_annotation_metadata-self.csv",
                            sep = ",",row.names = 1,header = T)
    D7_1_meta <- read.table(file = "./inputdata/PublicData/StereoSeq/JEV/annotation/THU66_singler_annotation_metadata-self.csv",
                            sep = ",",row.names = 1,header = T)
    CTRL_1_meta$grp <- "CTRL_1"
    D3_1_meta$grp <- "D3_1"
    D5_1_meta$grp <- "D5_1"
    D7_1_meta$grp <- "D7_1"
    CTRL_1 <- subset(CTRL_1,cells = rownames(CTRL_1_meta))
    D3_1 <- subset(D3_1,cells = rownames(D3_1_meta))
    D5_1 <- subset(D5_1,cells = rownames(D5_1_meta))
    D7_1 <- subset(D7_1,cells = rownames(D7_1_meta))
    table(rownames(CTRL_1@meta.data) %in% rownames(CTRL_1_meta))
    table(rownames(D3_1@meta.data) %in% rownames(D3_1_meta))
    table(rownames(D5_1@meta.data) %in% rownames(D5_1_meta))
    table(rownames(D7_1@meta.data) %in% rownames(D7_1_meta))
    CTRL_1_meta[c("x","y")] <- CTRL_1@meta.data[,c("x","y")]
    CTRL_1@meta.data <- CTRL_1_meta
    D3_1_meta[c("x","y")] <- D3_1@meta.data[,c("x","y")]
    D3_1@meta.data <- D3_1_meta
    D5_1_meta[c("x","y")] <- D5_1@meta.data[,c("x","y")]
    D5_1@meta.data <- D5_1_meta
    D7_1_meta[c("x","y")] <- D7_1@meta.data[,c("x","y")]
    D7_1@meta.data <- D7_1_meta
  }else{
    meta_data <- read.table(file = "./inputdata/PublicData/StereoSeq/JEV/annotation/251216_all_12_samples_metadata_info.csv",
                            sep = ",",header = T,row.names = 1) # 陈旗后给的meta数据是正确的
    meta_data <- meta_data %>%
      mutate(new_nm = paste0(new_sample,"_",imagecol,"-",imagerow))
    table(meta_data$sample)
    table(meta_data$Region)
    table(meta_data$singleR_main)
    rownames(meta_data) <- meta_data$new_nm
  }

  #>
  stRNA <- merge(CTRL_1,list(D3_1,D5_1,D7_1),add.cell.ids = c("CTRL_1","D3_1","D5_1","D7_1"))
  stRNA <- JoinLayers(stRNA)
  gc()

  stRNA@meta.data$new_samp <- meta_data[rownames(stRNA@meta.data),"new_sample"]
  stRNA@meta.data$new_tissue <- meta_data[rownames(stRNA@meta.data),"Region"]
  stRNA@meta.data$new_cell <- meta_data[rownames(stRNA@meta.data),"singleR_main"]
  stRNA@meta.data$grp <- meta_data[rownames(stRNA@meta.data),"group"]
  table(stRNA@meta.data$new_samp,useNA = "always")
  table(stRNA@meta.data$new_tissue,useNA = "always")
  table(stRNA@meta.data$new_cell,useNA = "always")
  table(stRNA@meta.data$grp,useNA = "always")

  dir.create("./inputdata/PublicData/StereoSeq/JEV/process/",showWarnings = F,recursive = T)
  saveRDS(stRNA,file = "./inputdata/PublicData/StereoSeq/JEV/process/stRNA.rds")

}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/StereoSeq/JEV/process/stRNA_JEV.rds")
  stRNA@meta.data$new_samp %>% unique()
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  # stRNA <- subset(stRNA, subset = new_samp == "D5_1") # D5_1一个样本
  # stRNA <- subset(stRNA, subset = new_samp == "D3_1")
  gc()
  meta_data <- stRNA@meta.data
  STID_obj <- as.STID(stRNA,samp_colnm = "new_samp",celltype_colnm = "new_cell",samp_grp_colnm = "grp",
                    host_org = "mouse",pathogen_grp = "virus",pathogen_org = "JEV",
                    data_format = "square_grid",data_platform = "StereoSeq", binsize = 35)
  # pathogen_genes <- c("NS5","C","NS3","NS1","E","Prm","NS4aAlt","NS4bAlt","NS2a","NS2b") # 只有包含至少两个 JEV 基因或三个病毒读数的 bin 才被视为 JEV 阳性
  pathogen_genes <- c("NS5")
  # host_genes <- c("Cxcl10","Ifitm3","Isg15","Irf7","Ccl5","Ccl2")
  host_genes <- c("Ccl5","Irf7")
  rm(meta_data,stRNA);gc()
  
  STID_obj_bin35 <- readRDS("F:/Scientific_research/ShiXiaoFeng-InfectiousST/01_Project/2025/20250124_STID/Analysis/rds_old/202605_article/STID_obj_detect_JEV_D5_1.rds")
  gc()
}
#> 组织类型介绍
# CB	Cerebellum	小脑	负责运动控制、平衡。
# CHP	Choroid plexus	脉络丛	产生脑脊液的地方，也是血脑屏障的重要组成部分。
# CNU	Cerebral Nuclei	大脑核团	包含纹状体（Striatum）等深部核团，通常是病毒感染的重灾区。
# CTX	Cortex	大脑皮层	负责高级认知、感觉处理。包含分层结构（如L1-L6）。
# FB	Fiber tracts	纤维束	（重点） 白质区域，由神经轴突组成。这就是你之前问的少突胶质细胞聚集的地方。
# HB	Habenula	缰核	连接前脑和中脑的“反奖赏系统”枢纽，与抑郁、压力有关。
# HPF	Hippocampal formation	海马结构	包含海马体（HIP）和齿状回（DG），负责记忆和空间导航。
# HY	Hypothalamus	下丘脑	负责激素调节、体温、睡眠等。
# MB	Midbrain	中脑	负责视觉、听觉、运动控制（如黑质）。
# MB_HB	Midbrain-Hindbrain boundary	中脑-后脑边界	发育上的重要界限区域。
# MEN	Meninges	脑膜	包裹大脑的保护层（硬脑膜、蛛网膜等）。
# OLF	Olfactory areas	嗅觉区	包含嗅球（OB），负责处理气味信息。
# TH	Thalamus	丘脑	感觉信息的中继站，传入大脑皮层前必经之地。
# UK	Unknown / Unassigned	未知/未分类	测序数据中无法明确比对到特定脑区的细胞或点。

# 原文献中介绍
# TH	Thalamus	丘脑
# CTX	Cerebral cortex	大脑皮层
# CNU	Cerebral nuclei	大脑核团 / 基底核
# HB	Hindbrain	后脑：除了小脑后的延髓、脑桥
# FB	Fiber tracts	神经纤维束
# CHP	Choroid plexus	脉络丛
# MB	Midbrain	中脑
# HY	Hypothalamus	下丘脑
# CB	Cerebellum	小脑
# HPF	Hippocampal formation	海马结构
# MEN	Meninges	脑膜
# OLF	Olfactory bulb	嗅球
# UK	Unknown	未知组织


#> 细胞类型介绍
# Astrocytes（星形胶质细胞）
# Microglia（小胶质细胞）
# Oligodendrocytes（少突胶质细胞）
# Adipocytes（脂肪细胞）
# Endothelial cells（内皮细胞）


## 3.2 bin100 ####
# bin100从gem 1转换而来，gem代码在当前目录下，h5ad -> rds代码在此处
if(0){
  # 无
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/StereoSeq/JEV/process/stRNA_JEV_bin100.rds")
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  gc()
  meta_data <- stRNA@meta.data
  STID_obj_bin100 <- as.STID(stRNA,samp_colnm = "orig.ident",celltype_colnm = "orig.ident",samp_grp_colnm = "orig.ident",
                             host_org = "mouse",pathogen_grp = "virus",pathogen_org = "JEV",
                             data_format = "square_grid",data_platform = "StereoSeq", binsize = 100)
  pathogen_genes <- c("NS5")
  host_genes <- c("Ccl5","Irf7")
  rm(meta_data,stRNA);gc()
}


## 3.3 bin200 ####
# bin200从gem 1转换而来，gem代码在当前目录下，h5ad -> rds代码在此处
if(0){
  # 无
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/StereoSeq/JEV/process/stRNA_JEV_bin200.rds")
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  gc()
  meta_data <- stRNA@meta.data
  STID_obj_bin200 <- as.STID(stRNA,samp_colnm = "orig.ident",celltype_colnm = "orig.ident",samp_grp_colnm = "orig.ident",
                      host_org = "mouse",pathogen_grp = "virus",pathogen_org = "JEV",
                      data_format = "square_grid",data_platform = "StereoSeq", binsize = 200)
  pathogen_genes <- c("NS5")
  host_genes <- c("Ccl5","Irf7")
  rm(meta_data,stRNA);gc()
}


## 3.4 cellbin ####
# cellbin直接从云平台下载的h5ad，h5ad -> rds代码在此处
if(0){
  # 无
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/StereoSeq/JEV/process/stRNA_JEV_cellbin.rds")
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  gc()
  meta_data <- stRNA@meta.data
  STID_obj_cellbin <- as.STID(stRNA,samp_colnm = "orig.ident",celltype_colnm = "brain_region",samp_grp_colnm = "orig.ident",
                      host_org = "mouse",pathogen_grp = "virus",pathogen_org = "JEV",
                      data_format = "single_cell",data_platform = "StereoSeq", binsize = 1)
  pathogen_genes <- c("NS5")
  host_genes <- c("Ccl5","Irf7")
  rm(meta_data,stRNA);gc()
  # Plot_Spatial(STID_obj,col = list(dis = COLOR_LIST$PALETTE_WHITE_BG,con = NULL),pt_size = 0.75)
}


### 4.TB ####
## 4.1 Mouse ####
# 原始数据从gem转换过来：F:\Scientific_research\ShiXiaoFeng-InfectiousST\01_Project\2024\20240924_结核项目
if(0){
  # process
  rds_list <- list()
  file_paths <- list.files("./inputdata/PublicData/StereoSeq/TB/",full.names = T)
  file_names <- gsub(".rds","",basename(file_paths))
  for (i in 1:length(file_paths)) {
    i_file_name <- file_names[i]
    i_file_path <- file_paths[i]
    print(i_file_name)
    i_rds <- readRDS(i_file_path)
    i_rds@meta.data$sample <- i_file_name
    rds_list[i_file_name] <- i_rds
  }
  sample_vector <- c(
    "Y00583M8" = "WT_1",
    "Y00583MA" = "WT_2",
    "Y00583N8" = "PI1d_1",
    "Y00583N9" = "PI1d_2",
    "Y00583F1" = "PI4W_1",
    "Y00583F6" = "PI4W_2",
    "Y00583N2" = "PI8W_1",
    "Y00583N6" = "PI8W_2"
  )
  names(rds_list) <- sample_vector[names(rds_list)]
  rds_list <- rds_list[sample_vector]
  stRNA <- merge(rds_list[[1]],rds_list[2:8],add.cell.ids = sample_vector)
  stRNA <- JoinLayers(stRNA)
  stRNA@meta.data$sample <- sample_vector[stRNA@meta.data$sample]
  meta_data <- stRNA@meta.data
  table(stRNA@meta.data$sample )
  dir.create("./inputdata/PublicData/StereoSeq/TB/process/",showWarnings = F,recursive = T)
  saveRDS(stRNA,file = "./inputdata/PublicData/StereoSeq/TB/process/stRNA.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/StereoSeq/TB/process/stRNA_TB.rds")
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  table(stRNA@meta.data$sample)
  # stRNA <- subset(stRNA, subset = sample %in% c("PI1d_1","PI4W_1","PI8W_1"))
  stRNA <- subset(stRNA, subset = sample %in% c("WT_1","PI1d_1","PI4W_1","PI8W_1"))
  meta_data  <- stRNA@meta.data
  table(meta_data$sample)
  STID_obj <- as.STID(stRNA,samp_colnm = "sample",celltype_colnm = "sample",samp_grp_colnm = "sample",
                    host_org = "mouse",pathogen_grp = "bacteria",pathogen_org = "Mtb",
                    data_format = "square_grid",data_platform = "StereoSeq")
  meta_data <- STID_obj@meta.data
  pathogen_genes <- grep("^Rv",rownames(STID_obj),value = T)
  length(pathogen_genes)
  # pathogen_genes <- c()
  # host_genes <- c("Ighm","Ighg2c","Igha","Igkc","Scgb1a1","Ttn","Rvnr01","Rvnr02")
  host_genes <- c("Igkc","Rvnr01","Rvnr02")
}

## 4.2 TB_human ####
# gem -> h5ad -> rds：F:\Scientific_research\ShiXiaoFeng-InfectiousST\01_Project\2024\20240924_结核项目
if(0){
  # process
  rds_list <- list()
  file_paths <- list.files("./inputdata/PublicData/StereoSeq/TB_human/",full.names = T,pattern = ".rds$")
  file_names <- gsub(".rds","",basename(file_paths))
  for (i in 1:length(file_paths)) {
    i_file_name <- file_names[i]
    i_file_path <- file_paths[i]
    print(i_file_name)
    i_rds <- readRDS(i_file_path)
    i_rds@meta.data$sample <- i_file_name
    rds_list[i_file_name] <- i_rds
  }
  sample_vector <- c(
    "H01203F1" = "C_3",
    "H01203D5" = "N_18",
    "H01203C3" = "O_18",
    "H01203C1" = "N_13",
    "H01203A3" = "O_15",
    "G01138F3" = "C_9"
  )
  names(rds_list) <- sample_vector[names(rds_list)]
  rds_list <- rds_list[sample_vector]
  stRNA <- merge(rds_list[[1]],rds_list[2:6],add.cell.ids = sample_vector)
  stRNA <- JoinLayers(stRNA)
  stRNA@meta.data$sample <- sample_vector[stRNA@meta.data$sample]
  meta_data <- stRNA@meta.data
  table(stRNA@meta.data$sample )
  dir.create("./inputdata/PublicData/StereoSeq/TB_human/process/",showWarnings = F,recursive = T)
  saveRDS(stRNA,file = "./inputdata/PublicData/StereoSeq/TB_human/process/stRNA_TB_human.rds")
  
  stRNA <- readRDS(file = "./inputdata/PublicData/StereoSeq/TB_human/process/stRNA_TB_human.rds")
  stRNA <- subset(stRNA, subset = sample %in% c("N_18"))
  saveRDS(stRNA,file = "./inputdata/PublicData/StereoSeq/TB_human/process/stRNA_TB_human.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/StereoSeq/TB_human/process/stRNA_TB_human.rds")
  stRNA <- suppressMessages({
    UpdateSeuratObject(stRNA)
  })
  table(stRNA@meta.data$sample)
  # stRNA <- subset(stRNA, subset = sample %in% c("PI1d_1","PI4W_1","PI8W_1"TB_human))
  # stRNA <- subset(stRNA, subset = sample %in% c("WT_1","PI4W_1","PI8W_1"))
  meta_data  <- stRNA@meta.data
  table(meta_data$sample)
  pathogen_genes <- grep("^Rv",rownames(stRNA),value = T)
  STID_obj <- as.STID(stRNA,samp_colnm = "sample",celltype_colnm = "sample",samp_grp_colnm = "sample",
                      host_org = "human",pathogen_grp = "bacteria",pathogen_org = "Mtb",
                      pathogen_genes = pathogen_genes,
                      data_format = "square_grid",data_platform = "StereoSeq")
  meta_data <- STID_obj@meta.data
  host_genes <- c("IGKC")
  length(pathogen_genes)
  # host_genes <- c(
  #   "IGHM",
  #   "IGHG1", "IGHG2", "IGHG3", "IGHG4",
  #   "IGHA1", "IGHA2",
  #   "IGKC",
  #   "SCGB1A1",
  #   "TTN"
  # )
}


### 5.HPV × ####
# 效果不好：not useful
if(0){
  # process
  TJH11 <- readRDS(file = "./inputdata/PublicData/StereoSeq/HPV/TJH11_seurat.rds") # bin50,包含多个矩阵
  TJH34 <- readRDS(file = "./inputdata/PublicData/StereoSeq/HPV/TJH34_seurat.rds")
  TJH46 <- readRDS(file = "./inputdata/PublicData/StereoSeq/HPV/TJH46_seurat.rds")
  DefaultAssay(TJH11) <- "Spatial"
  DefaultAssay(TJH34) <- "Spatial"
  DefaultAssay(TJH46) <- "Spatial"
  TJH11@assays$SCT <- NULL
  TJH11@assays$microbe <- NULL
  TJH34@assays$SCT <- NULL
  TJH34@assays$microbe <- NULL
  TJH46@assays$SCT <- NULL
  TJH46@assays$microbe <- NULL
  stRNA <- merge(TJH11,list(TJH34,TJH46),add.cell.ids = c("TJH11","TJH34","TJH46"))
  stRNA[["Spatial"]] <- as(object = stRNA[["Spatial"]], Class = "Assay5") # V4 -> V5
  stRNA[["predictions"]] <- as(object = stRNA[["predictions"]], Class = "Assay5") # V4 -> V5
  dir.create("./inputdata/PublicData/StereoSeq/HPV/process/",showWarnings = F,recursive = T)
  saveRDS(stRNA,file = "./inputdata/PublicData/StereoSeq/HPV/process/stRNA_HPV.rds")
}else{
  # read rds
  stRNA <- readRDS(file = "./inputdata/PublicData/StereoSeq/HPV/process/stRNA_HPV.rds")
  meta_data <- stRNA@meta.data
  STID_obj <- as.STID(stRNA,samp_colnm = "orig.ident",celltype_colnm = "orig.ident",samp_grp_colnm = "orig.ident",
                      host_org = "mouse",pathogen_grp = "virus",pathogen_org = "HPV",
                      data_format = "square_grid",data_platform = "StereoSeq")
  STID_obj@meta.data$HPV <- ifelse(is.na(STID_obj@meta.data$HPV),0,STID_obj@meta.data$HPV)
  pathogen_genes <- c("E1","E2","E5","E6","E7","L1","L2")
}


### 三、Gene_Geneset ####
library(vroom)
read_txt_as_df_vroom <- function(path, return_filenames_only = FALSE, ...) {
  entries <- list.files(
    path,
    full.names = TRUE,
    include.dirs = TRUE,
    recursive = FALSE,
    no.. = TRUE  # 排除 . 和 ..
  )

  if (length(entries) == 0) return(list())

  names_entries <- basename(entries)
  result <- list()

  for (i in seq_along(entries)) {
    entry <- entries[i]
    name <- names_entries[i]

    if (dir.exists(entry)) {
      result[[name]] <- read_txt_as_df_vroom(
        entry,
        return_filenames_only = return_filenames_only,
        ...
      )
    } else if (grepl("\\.txt$", entry, ignore.case = TRUE)) {
      clean_name <- tools::file_path_sans_ext(name)

      if (return_filenames_only) {
        result[[clean_name]] <- name
      } else {
        df <- tryCatch({
          vroom(entry, ...)
        }, error = function(e) {
          warning("Failed to read '", entry, "': ", conditionMessage(e))
          NULL
        })
        result[[clean_name]] <- df
      }
    }
  }

  return(result)
}

Gene_Geneset <- read_txt_as_df_vroom(
  path = "./inputdata/Gene_Geneset/",
  delim = "\t"
)
saveRDS(Gene_Geneset,file = "./rds/Gene_Geneset.rds")

file_structure <- read_txt_as_df_vroom(
  path = "./inputdata/Gene_Geneset/",
  delim = "\t",
  return_filenames_only = TRUE
)
str(file_structure)






