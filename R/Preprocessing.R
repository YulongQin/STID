

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Seurat processing and SingleR annotation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Standard Seurat Analysis Pipeline for Spatial Transcriptomics
#'
#' Performs a complete Seurat analysis pipeline for spatial transcriptomics data,
#' including normalization, feature selection, scaling, PCA, clustering, and
#' UMAP/t-SNE dimensionality reduction. Also generates spatial feature plots.
#'
#' @param seurat_obj A Seurat object containing spatial transcriptomics data
#' @param data_format Character, data format - "stereo" (Stereo-seq) or "visium"
#'        (default: "stereo")
#' @param data_type Character, data type - "stRNA" (spatial) or "scRNA" (single-cell)
#'        (default: "stRNA")
#' @param resolution_index Numeric vector, clustering resolutions to test
#'        (default: seq(0.1, 1.3, 0.2))
#' @param runTSNE_index Logical, whether to perform t-SNE (default: FALSE)
#' @param assay_nm Character, assay name (default: NULL, auto-detected)
#'
#' @return A Seurat object with processed data and cluster annotations
#'
#' @import Seurat
#' @import patchwork
#' @importFrom viridis viridis
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run standard pipeline on spatial data
#' seurat_obj <- Seurat_pipeline(
#'   seurat_obj = my_spatial_data,
#'   data_format = "stereo",
#'   resolution_index = seq(0.2, 1.0, 0.2)
#' )
#' }
Seurat_pipeline <- function(seurat_obj = NULL,
                            data_format = "stereo",data_type = "stRNA",
                            resolution_index= seq(0.1,1.3,0.2),
                            runTSNE_index = FALSE,
                            assay_nm = NULL
) {

  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_check()
  if(have_null(seurat_obj)){
    clog_error("The seurat_obj, ref_obj, seurat_colnm and ref_colnm must be provided")
  }
  if(data_type == "stRNA"){
    assay_nm <- "Spatial"
    image_index <- TRUE
  }else if(data_type == "scRNA"){
    assay_nm <- "RNA"
    image_index <- FALSE
  }else{
    clog_error("data_type is not correct")
  }
  # >>> End check

  # >>> Start main pipeline
  clog_normal("SpatialFeaturePlot")
  p1 <- SpatialFeaturePlot(seurat_obj,
                           features = c("nFeature_Spatial"), # "Alb"
                           pt.size.factor = 1.5,
                           # alpha = c(0.1, 1),
                           stroke = 0,
                           max.cutoff = "q95", # 可以是具体的数值
                           min.cutoff = "q5") +
    theme_void() +
    labs(title = "nFeature_Spatial") +
    # scale_fill_gradientn(colours = c("blue", "white", "red")) +
    scale_fill_gradientn(colours = viridis::viridis(20)) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
  p2 <- SpatialFeaturePlot(seurat_obj,
                           features = c("nCount_Spatial"), # "Alb"
                           pt.size.factor = 1.5,
                           # alpha = c(0.1, 1),
                           stroke = 0,
                           max.cutoff = "q95", # 可以是具体的数值
                           min.cutoff = "q5") +
    theme_void() +
    labs(title = "nCount_Spatial") +
    # scale_fill_gradientn(colours = c("blue", "white", "red")) +
    scale_fill_gradientn(colours = viridis::viridis(20)) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
  print(p1 + p2)

  #
  clog_normal("Perform NormalizeData")
  seurat_obj <- NormalizeData( # 归一化+log
    seurat_obj,
    assay = assay_nm, # 对哪个assay标准化
    normalization.method = "LogNormalize", # 'CLR','RC'
    scale.factor = 10000, # 一个细胞的UMI大概这么多
    margin = 1 # 'CLR'使用:1行 2列
  )
  clog_normal("Perform FindVariableFeatures")
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst", # 根据均值、方差关系拟合曲线，对特征值标准化，降序排列取n个基因
    nfeatures = 3000 # HVF数：2000-5000,3000个比较合适，其它gene不进行后续PCA等分析
  )
  clog_normal("Perform ScaleData")
  seurat_obj <- ScaleData(
    seurat_obj
  )
  clog_normal("Perform RunPCA")
  seurat_obj <- RunPCA(seurat_obj,
                       assay = assay_nm, # "RNA"、"SCT"
                       # features = VariableFeatures(seurat_obj), # 默认为所有选取的HVF
                       features = NULL, # 默认为所有选取的HVF
                       npcs = 50, # PCA总数
                       rev.pca = FALSE,
                       weight.by.var = TRUE,
                       verbose = F,
                       ndims.print = 1:5,
                       nfeatures.print = 30,
                       reduction.name = "pca",
                       reduction.key = "PC_",
                       seed.use = 42)
  clog_normal("Perform FindClusters")
  seurat_obj <- seurat_obj %<>%
    FindNeighbors(., dims = 1:10) %>%  # 结果保存在seurat_obj@graph下面，dim可以和umap不同
    FindClusters(., resolution = resolution_index) %>% # 生成的结果在seurat_clusters.*中，但同时也会生成seurat_clusters
    FindClusters(., resolution = 0.5,
                 cluster.name = 'seurat_clusters', # 生成的结果在seurat_clusters.*中，但同时也会生成seurat_clusters
                 algorithm = 1)
  if(runTSNE_index){
    clog_normal("Perform RunTSNE")
    seurat_obj <- RunTSNE(seurat_obj,dims = 1:20,check_duplicates=F)
  }
  clog_normal("Perform RunUMAP")
  seurat_obj <- RunUMAP(seurat_obj,
                        dims = 1:10, # 要使用的维度数。一般选择1:20,NULL为使用全部
                        reduction = "pca",
                        # reduction = "harmony", # 基于哪个降维结果做umap
                        # reduction.key = "UMAP_", # 可以保存多个umap
                        n.neighbors = 30L, # umap计算时局部近似流形结构中使用的相邻点的数量。一般应在5到50之间，越大则cluster之间交集越多
                        min.dist = 0.3, # 调整umap的坐标，点之间小距离，距离越小越往每个cluster的质心聚集，则图更加分散式的聚集。默认0.3， 0.001-0.5，可视化时值大一点更好看？
                        n.components=2, # 生成的umap维度数，如果是三维图需要是3
                        seed.use = 42)
  clog_normal("Perform SpatialPlot")
  p1 <- SpatialPlot(seurat_obj,
                    group.by = "seurat_clusters", # 等价于SpatialDimPlot
                    pt.size.factor = 1,
                    stroke = 0) +
    theme_void() +
    labs(title = "seurat_clusters") +
    scale_fill_manual(values = color_stRNA1) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
  print(p1)
  # >>> End main pipeline

  # >>> End
  clog_end()

  # >>> return
  return(seurat_obj)
}

#' Annotate Cell Types Using SingleR
#'
#' Performs automated cell type annotation using the SingleR package. Uses a
#' reference dataset (e.g., from celldex) to annotate cells based on gene
#' expression profiles.
#'
#' @param seurat_obj A Seurat object to annotate
#' @param ref_obj Reference dataset (matrix or SummarizedExperiment) with
#'        cell type labels
#' @param seurat_colnm Character, column name in Seurat metadata containing
#'        cluster IDs for aggregation
#' @param ref_colnm Character, column name in reference object containing
#'        cell type labels
#' @param species_index Character, species for built-in references (currently unused)
#' @param assay_nm Character, assay to use (default: "Spatial")
#' @param layer_nm Character, layer/slot to use (default: "data")
#'
#' @return A Seurat object with added 'anno_SingleR' column in metadata
#'
#' @import Seurat
#' @importFrom SummarizedExperiment colData
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Annotate using built-in reference
#' library(SingleR)
#' ref <- HumanPrimaryCellAtlasData()
#'
#' seurat_obj <- anno_SingleR(
#'   seurat_obj = my_seurat,
#'   ref_obj = ref,
#'   seurat_colnm = "seurat_clusters",
#'   ref_colnm = "label.main"
#' )
#' }
anno_SingleR <- function(seurat_obj = NULL,
                         ref_obj = NULL,
                         seurat_colnm = NULL,
                         ref_colnm = NULL,
                         species_index = NULL,
                         assay_nm = "Spatial",
                         layer_nm = "data"
) {

  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!requireNamespace("SingleR", quietly = TRUE)) {
    clog_error("Package 'SingleR' is required but is not installed.")
  }
  if(have_null(seurat_obj, ref_obj, seurat_colnm, ref_colnm)){
    clog_error("The seurat_obj, ref_obj, seurat_colnm and ref_colnm must be provided")
  }
  clog_normal(paste0("Your seurat_colnm is ",seurat_colnm))
  clog_normal(paste0("Your ref_colnm is ",ref_colnm))
  if(! (is.matrix(ref_obj) || inherits(ref_obj, "SummarizedExperiment"))){
    clog_error("The ref_obj must be a matrix or a SummarizedExperiment")
  }
  # >>> End check

  # >>> Start main pipeline
  clog_step("Normalize seurat_obj")
  if("data" %in% names(seurat_obj@assays[[assay_nm]]@layers)){
    clog_normal("Find data in the assay")
  }else{
    clog_warn("No data in the assay")
    clog_normal("Perform NormalizeData")
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  obj_exp <- GetAssayData(seurat_obj, slot=layer_nm)

  #
  clog_step("Perform SingleR")
  clog_normal("Your seurat label is as follows:")
  print(table(seurat_obj@meta.data[[seurat_colnm]]))
  clog_normal("Your ref label is as follows:")
  print(table(colData(ref_obj)[,ref_colnm]))
  res_singleR <- SingleR::SingleR(
    test = obj_exp,
    ref = ref_obj,
    labels = colData(ref_obj)[,ref_colnm],
    clusters = seurat_obj@meta.data[[seurat_colnm]]
  )
  res_df <- data.frame(ClusterID=rownames(res_singleR),
                       anno_SingleR=res_singleR$labels,
                       stringsAsFactors = FALSE)
  seurat_obj@meta.data$anno_SingleR <- res_df[match(seurat_obj@meta.data[[seurat_colnm]], res_df$ClusterID), 'anno_SingleR']
  clog_normal("Your SingleR annotation is added in the meta.data with the name 'SingleR'")
  clog_normal("The annotation result is as follows:")
  print(table(seurat_obj@meta.data$anno_SingleR))
  # >>> End main pipeline

  # >>> End
  clog_end()

  # >>> return
  return(seurat_obj)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# paese the gtf
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Parse GTF File into Tidy Data Frame
#'
#' Parses a GTF (Gene Transfer Format) file and converts it into a tidy
#' data frame with key-value attributes properly separated into columns.
#'
#' @param gtf_file Character, path to the GTF file
#' @param fil_label Character, feature type to filter (e.g., "gene", "exon")
#'        (default: "gene")
#'
#' @return A data frame with GTF columns (seqname, source, feature, start, end,
#'         score, strand, frame) plus additional columns for each attribute
#'
#' @importFrom data.table fread
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Parse a GTF file and keep only gene features
#' gtf_data <- parse_gtf(
#'   gtf_file = "path/to/annotation.gtf",
#'   fil_label = "gene"
#' )
#'
#' # View the parsed data
#' head(gtf_data)
#' }
parse_gtf <- function(gtf_file = NULL, fil_label = "gene"){

  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_check()
  clog_normal(paste0("Your fil_label is ",fil_label))
  # >>> End check

  # >>> Start main pipeline
  clog_step(" read gtf file")
  gtf_df <- fread(gtf_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>%
    as.data.frame()
  clog_normal(paste("The number of rows in gtf file is: ", nrow(gtf_df)))
  setnames(gtf_df, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
  if(!is.null(fil_label)){
    gtf_df <- gtf_df %>%
      filter(grepl(fil_label, feature))
  }
  clog_normal(paste("The number of rows in gtf file after filtering is: ", nrow(gtf_df)))
  gtf_df <- gtf_df %>%
    mutate(attribute = gsub("; ",";",attribute)) %>%
    mutate(attribute = gsub(";$","",attribute)) %>%
    mutate(attribute = gsub(" \"","%\"",attribute)) %>%
    mutate(attribute = gsub(" (\\d+);", "%\\1;", attribute, perl = TRUE))

  #
  clog_step("separate_wider_delim")
  res_gtf <- gtf_df[,-c(1:8),drop = FALSE]
  res_gtf <- separate_wider_delim(res_gtf, cols = "attribute", delim = regex("[%;]"),
                                  names_sep = "_A",
                                  too_few = "align_start" #
  )
  if(ncol(res_gtf) %% 2 == 1){
    table(res_gtf[,ncol(res_gtf)-1],useNA = "always") %>% print()
    table(res_gtf[,ncol(res_gtf)],useNA = "always") %>% print()
    error_df <- res_gtf %>%
      filter(!is.na(res_gtf[,ncol(res_gtf)]))
    clog_error("The number of columns is not odd number! please check the delim of separate_wider_delim function!")
  }
  clog_normal(paste("The number of key-value pairs in gtf file is: ", ncol(res_gtf)/2))
  res_gtf <- res_gtf %>%
    mutate(across(where(is.character), ~gsub("\"", "", .))) %>%
    mutate(across(where(is.character), ~ifelse(. == "", NA, .)))

  #
  clog_step("split df")
  res_list <- lapply(1:(ncol(res_gtf)/2), function(i){
    i_df <- res_gtf[, c(2*i-1, 2*i)]
    colnames(i_df) <- c("key", "value")
    i_df <- i_df %>%
      mutate(id = 1:n()) %>%
      filter(rowSums(is.na(i_df)) < 2)
  })

  #
  clog_step("pivot_wider")
  res_df <- do.call(rbind, res_list)
  table(res_df$key) %>% head()
  res_df2 <- res_df %>%
    pivot_wider(names_from = key, values_from = value,
                values_fn = ~paste(.x, collapse = ";")) %>%
    arrange(id) %>%
    select(-id) %>%
    cbind(gtf_df[, 1:8], .)
  # >>> End main pipeline

  # >>> End
  clog_end()

  # >>> return
  return(res_df2)
}

