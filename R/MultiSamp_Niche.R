

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalNicheOSE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate Niche Organization and Spatial Entropy (OSE)
#'
#' Performs comprehensive spatial niche analysis including region segmentation,
#' entropy calculation, and visualization of spatial organization patterns.
#' This function processes spatial transcriptomics data to identify organizational
#' niches and calculate adjusted spatial entropy metrics.
#'
#' @param x Placeholder parameter for future expansion
#'
#' @return NULL (invisible), generates multiple output files and visualizations
#'
#' @details
#' The function performs the following analyses:
#' \enumerate{
#'   \item **Data Preprocessing**: Loads and normalizes spatial transcriptomics data
#'   \item **Spatial Region Segmentation**: Uses iterative k-means-like algorithm to partition tissue into spatial regions
#'   \item **Spatial Entropy Calculation**: Computes adjusted entropy based on cell type diversity in local neighborhoods
#'   \item **Cell Type Proportion Analysis**: Calculates cell type composition for each region
#'   \item **Boundary Detection**: Identifies boundaries between different regions
#'   \item **Visualization**: Generates multiple plots including region maps, cell type distributions, pie charts, and entropy heatmaps
#' }
#'
#' @section Parameters:
#' The function includes extensive parameter settings based on analysis mode:
#' \itemize{
#'   \item \code{mode}: Analysis resolution - "bin50", "bin200", or "pub"
#'   \item \code{binsize}: Spatial bin size (50, 200, or 1 based on mode)
#'   \item \code{window}: Window size for initial region segmentation
#'   \item \code{m}: Weight parameter balancing spatial distance and expression distance
#'   \item \code{nPC}: Number of principal components to use
#'   \item \code{min_spot_num}: Minimum spots per region for inclusion
#' }
#'
#' @section Outputs:
#' Generates files in outputdata/ and photo/ directories:
#' \itemize{
#'   \item \code{*_center_df.txt}: Coordinates of region centers
#'   \item \code{*_data.txt}: Cell-level assignments to regions
#'   \item \code{*_entropy_results.txt}: Raw entropy calculations
#'   \item \code{*_entropy_results_fil.txt}: Filtered entropy results
#'   \item \code{*_celltype_percent.txt}: Cell type proportions per region
#'   \item \code{*_celltype_percent_fil.txt}: Filtered cell type proportions
#'   \item \code{*_segment.txt}: Region boundary coordinates
#'   \item Multiple visualization PDFs and PNGs
#' }
#'
#' @import Seurat
#' @import Matrix
#' @import ggnewscale
#' @import paletteer
#' @import ggforce
#' @import furrr
#' @import future
#' @importFrom scales percent rescale
#' @importFrom proxy dist
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run niche analysis with default settings
#' CalNicheOSE()
#'
#' # Modify parameters within the function before running
#' # Set mode <- "bin200" for lower resolution analysis
#' # Set only_plot <- TRUE to regenerate plots from existing data
#' }
CalNicheOSE <- function(x){

  ### 0.参数设置 ####
  mode <- "bin50" # "bin50"/"bi200"/"pub"
  # mode <- "bin200"
  # mode <- "pub"
  only_plot <- T
  alpha_index <- 0.8

  #
  values <- c(
    "Others" = "#DEEDFF",
    "Medulla" = "#FF3A00",
    "Junction" = "#B7E6A7",
    'Cortex' = "#006FAB"
  )

  if(1){
    options(future.globals.maxSize= 10*1024^3)
    plan(multicore, workers = 16)  # 重新运行一下就会释放占用的内存
    # plan(sequential)
    message("outside>>how many cores can use now: ", nbrOfWorkers())
  }
  if(0){
    sap_nm <- "chky_wite_score_tag" # thymus_sample_8_st
    input_dir <- paste0("/data/work/outdata_h5ad_rds/",mode,"_mingene_0_rawgem/rds/")
    input_txt_dir <- paste0("/data/work/outdata_entropy/",mode,"_mingene_0_rawgem/")
    output_dir <- paste0("/data/work/outdata_entropy/",mode,"_mingene_0_rawgem/")
    photo_dir <- paste0("/data/work/photo_entropy/",mode,"_mingene_0_rawgem/")
    # values <- c(
    #   "Others" = "#e5e7e9",
    #   "Medulla" = "#e74c3c",
    #   "Junction" = "#52be80",
    #   'Cortex' = "#2980b9"
    # )
  }else{
    sap_nm <- "HC-13-C_wite_score_tag" # thymus_sample_2_st,chky_wite_score_tag
    input_dir <- paste0("./inputdata/",mode,"_mingene_0_rawgem/")
    input_txt_dir <- paste0("./inputdata/",mode,"_mingene_0_rawgem/")
    output_dir <- paste0("./outputdata/",mode,"_mingene_0_rawgem/")
    photo_dir <- paste0("./photo/",mode,"_mingene_0_rawgem/")
  }
  dir.create(output_dir, showWarnings = F)
  dir.create(photo_dir, showWarnings = F)

  if(mode =="bin50"){
    binsize <- 50
    window <- 70   # 60是1.5mm, 原来是bin50=90,bin200=25:决定分区个数，目前的90会得到5*5即25的分区（但左上角分区中无组织，故最终为24个）
    m <- 0.2  # bin50=0.2,bin200=0.5:距离权重，越小表达距离的权重越大，建议在0.001-0.01范围内
    min_spot <- 1 # 最小的区域的spot数目
    nPC <- 15
    neighbor_range <- 1
    celltype_col <- 'celltype'
    index_10x <- F
    min_spot_num <- 1000

    p1_size <- 2
    p2_size <- 0.3 # 0.42
    p2_merge_r <- 15
    p3_merge_size <-9
  }else if(mode == "bin200"){
    binsize <- 200
    window <- 70/4   # bin50=90,bin200=25:决定分区个数，目前的90会得到5*5即25的分区（但左上角分区中无组织，故最终为24个）
    m <- 0.95  # bin50=0.2,bin200=0.5:距离权重，越小表达距离的权重越大，建议在0.001-0.01范围内
    min_spot <- 1 # 最小的区域的spot数目
    nPC <- 15
    neighbor_range <- 1
    celltype_col <- 'celltype'
    index_10x <- F
    min_spot_num <- 100

    p1_size = 2.4
    p2_size <- 2.4
    p2_merge_r <- 4
    p3_merge_size <- 10
  }else if(mode == "pub"){
    binsize <- 1
    # 根据面积：0.65*0.65*0.7=0.29575，我们的是1*1*0.9=0.9，所以是0.9/0.29575=3.04，我们的是36个区域，所以大致是一样的
    window <- 20 # 跟bin200有所不同？pub2/4/8，10x的是根据实际的小叶大小来大致分割的，每个片子，大概12-14个？sh
    m <- 0.95
    min_spot <- 1
    nPC <- 15
    neighbor_range <- 1
    celltype_col <- 'celltype'
    index_10x <- T
    min_spot_num <- 100

    p1_size <- 5 # pub2/4/7/8
    # p2_size <- 4.5 # pub2/4/8
    p2_size <- 4.7 # pub2/4/8
    p2_merge_r <- 6.5 # pub2/4/8，饼图
    # 熵点的大小，以no_grp为准，不是no_seg，但是也还要缩放，所以还是以no_seg为准
    # 其它为14，4为12
    p3_merge_size <- 14
  }else{
    print("no mode")
  }

  ### 1.load and preprocess the rds ####
  print("step1: load and preprocess the data")
  print(paste0("sap_nm: ",sap_nm))
  print(paste0("input_dir: ",input_dir))
  print(paste0("output_dir: ",output_dir))
  print(paste0("photo_dir: ",photo_dir))

  if(!only_plot){
    # seurat data
    seurat_obj <- readRDS(file = paste0(input_dir,"/",sap_nm,".rds"))
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj,verbose = F)

    # metadata
    if(index_10x){
      min_x <- min(seurat_obj@meta.data$x)
      min_y <- min(seurat_obj@meta.data$y)
      gap_x <- seurat_obj@meta.data$x %>% unique() %>% sort() %>% diff() %>% table() %>% names() %>% as.numeric() %>% sort() %>% .[2]
      print(paste0("gap_x: ",gap_x))
      gap_y <- seurat_obj@meta.data$y %>% unique() %>% sort() %>% diff() %>% table() %>% names() %>% as.numeric() %>% sort() %>% .[2]
      print(paste0("gap_y: ",gap_y))
      gap_select <- min(gap_x,gap_y)
      print(paste0("gap_x and gap_y ratio: ",gap_x/gap_y))
      seurat_obj@meta.data$x <- round((seurat_obj@meta.data$x-min_x+1)/gap_select) # 除以gap_select，而不是gap_x，因为x和y比例关系是正确的，不能矫正了
      seurat_obj@meta.data$y <- round((seurat_obj@meta.data$y-min_y+1)/gap_select)
    }
    meta_data <- seurat_obj@meta.data
    seurat_obj@meta.data <- data.frame(
      orig.ident = rep("sample", nrow(meta_data)),
      celltype = meta_data$tag,
      row = meta_data$y,
      col = meta_data$x
    )
    row.names(seurat_obj@meta.data) <- rownames(meta_data)
    meta_data <- seurat_obj@meta.data
    if(T){
      meta_data$row <- meta_data$row + binsize # 方式出现0的情况，后面构建稀疏矩阵index1 = T会报错
      meta_data$col <- meta_data$col + binsize # 这个坐标都是增加，不会影响实际的图形的相对位置
    }
    write.table(meta_data, file = paste0(output_dir,"/",sap_nm,"_meta_data.txt"),
                sep = "\t", quote = F, row.names = F)

    ### 2.空间区域 ####
    # 重新迭代从此开始

    ## 2.1 提取坐标和PC数据 ####
    print(paste0("step2: extract pc data ",sap_nm))
    pc_data <- seurat_obj@reductions[["pca"]]@cell.embeddings %>%
      as.data.frame() %>%
      mutate(cellid = rownames(meta_data))
    data <- data.frame(
      row = meta_data$row/binsize, # row和col可能为0
      col = meta_data$col/binsize,
      celltype = meta_data$celltype,
      cellid = rownames(meta_data)
    ) %>%
      left_join(pc_data, by = "cellid") %>%
      select(1:(nPC+4)) # row + col + celltype + cellid  + PC1-15
    row.names(data) <- data$cellid
    data <- data %>%
      mutate(
        label = NA,
        dist = Inf
      )
    min_col <- min(data$col)
    min_row <- min(data$row)

    ## 2.2 生成超像素块区域及中心点 ####
    create_group <- function(x, interval) {
      ceiling((x - min(x)) / interval + 0.0001) # 向上取整，将每个spot分为哪个区域中
    }
    center_df <- data %>%   # 重新迭代前要将这段跑一下，重新得到初始的center_df
      mutate(
        col_group = create_group(col, window),  # 先添加后删除
        row_group = create_group(row, window)
      ) %>%
      group_by(col_group, row_group) %>%
      summarise(
        spot_num = n(),
        across(starts_with("PC_"), mean) # 其余的列都丢弃，每个区域用所有细胞的均值PC值代表
      ) %>%
      ungroup() %>%
      mutate(
        col = as.numeric(col_group) * window + min_col-0.5*window, # 手动定义中心点
        row = as.numeric(row_group) * window + min_row-0.5*window
      ) %>%
      filter(spot_num > min_spot) %>%
      mutate(label = str_c("V", seq(n()))) %>%
      arrange(label) %>%
      select(c(-col_group,-row_group)) # 最终将col_group和row_group转为坐标了

    ## 2.3 迭代更新 ###
    print(paste0("step3: update ",sap_nm))
    res <- Inf
    step <- 1
    data_backup <- data
    if(binsize == 50){
      res_value = 10 # 这个值不变了，就这么定了，16过大了，5过小了
    }else{
      res_value = 1 # bin200和pub用这个都合适
    }
    while (res > res_value) {
      data <- data_backup

      # 计算每个区域与细胞的距离
      dist_s <- proxy::dist(
        data %>% select(col, row), # data源代码中为spots
        center_df %>% select(col, row), # center_df源代码中为centers
        method = "Euclidean"
      ) %>%
        unclass() %>%
        as.data.frame() %>%
        rownames_to_column("cellid") %>%
        pivot_longer(-cellid, names_to = "label", values_to = "d_s") #%>% # value就是距离
      #filter(d_s <= window)
      # print(dim(dist_s)) # 每个点距离每个中心点的实际距离

      dist_exp <- proxy::dist(
        data %>% select(starts_with("PC_")),
        center_df %>% select(starts_with("PC_")),
        method = "Euclidean"
      ) %>%
        unclass() %>%
        as.data.frame() %>%
        rownames_to_column("cellid") %>%
        pivot_longer(-cellid, names_to = "label", values_to = "d_e")
      # print(dim(dist_exp)) # 每个点距离每个中心点的表达距离

      # 综合dist_s和dist_exp
      dists <- dist_s %>%
        inner_join(dist_exp) %>%
        mutate(
          dist = sqrt(m * d_s^2 + (1-m)*d_e^2) # m是距离权重，也只有一个m参数
        ) %>%
        group_by(cellid) %>%
        slice_min(dist) %>% # 选出每个点距离最近的中心点
        select(cellid, label, dist)

      data <- data %>% # data在源代码中有时为df
        select(-c(label, dist)) %>%
        left_join(dists)

      # 根据新的点更新中心点
      center_df_new <- data %>%
        group_by(label) %>% # 前面的是group_by(col_group, row_group)
        summarise(
          spot_num = n(),
          across(c(col, row, starts_with("PC_")), mean)
        ) %>%
        ungroup() %>%
        arrange(label) %>%
        filter(!is.na(label)) # 没有对应的中心点区域

      center_df <- center_df %>%
        filter(label %in% center_df_new$label) %>%
        arrange(label)

      res <- as.matrix(center_df %>% select(-label)) - as.matrix(center_df_new %>% select(-label))
      res <- colSums(res)
      res <- sum(abs(res)) # 计算所有位置的偏差，如果小于10或1就停止
      print(
        paste("Step:", step, # 第几次迭代
              " Segment num:", nrow(center_df_new), # 迭代后后的中心点数，即分区数
              " NA spot: ", sum(is.na(data$label)),
              " Diff:", res) # 迭代后与迭代前的中心点的差异
      )
      center_df <- center_df_new
      step <- step + 1
    }
    write.table(center_df, file = paste0(output_dir,"/",sap_nm,"_center_df.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(data, file = paste0(output_dir,"/",sap_nm,"_data.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    ### 3.区域空间熵 ####
    print(paste0("step4: calculate entropy "))

    ## 计算矫正的系数
    combination <- function(n, k) {
      choose(n + k - 1, k)   # 总的niche可能性
    }
    # all_nicheType <- combination(4, 8)*4 # 如果周围都是点的话
    # 其实只要空白多就会熵高，哪怕不算这些不全的点，熵似乎会更高
    min_neighbor_num <- 3 # 3个点只是限制这是一个niche，实际上空白多似乎不会对熵有很大的影响
    all_nicheType <- sum(sapply(min_neighbor_num:8,function(x) combination(4, x)))*4 # 1920，考虑周围最少3个点的情况，确定了，就定为3个，
    expected_covered_species <- function(N, n) {
      N * (1 - (1 - 1/N)^n) # 指定抽样次数的期望覆盖种类数
    }
    celltype_mat <- Matrix::sparseMatrix(
      data$row, # 不能为0，
      data$col,
      x = as.numeric(factor(data$celltype)), # value值？
    )
    entropy_results <- data %>%
      select(row, col, label, cellid, celltype) %>% # 这个因子顺序不重要，只要是不同的数字，后面确定niche
      mutate(celltype = as.numeric(factor(celltype))) %>% # 行数为点数
      mutate(
        neighbors = map2(row, col, ~ {
          celltype_code <-  celltype_mat[
            # 周围8个点的范围
            max(0, .x - neighbor_range):min(.x + neighbor_range, dim(celltype_mat)[1]), # 不要低于0，不要高于最大值，即范围内
            max(0, .y - neighbor_range):min(.y + neighbor_range, dim(celltype_mat)[2])
          ] %>% as.vector()

          celltype_code <- celltype_code[which(celltype_code > 0)] # 0就是dgc中的值为0，即细胞类型0，即没有细胞
          celltype_code <- celltype_code[-match(celltype_mat[.x, .y], celltype_code)] # 去掉一个自己
          table(celltype_code)
        })) %>%
      group_by(label) %>%
      mutate(N_bin50_raw = n()) %>%
      ungroup() %>%
      # 至少3个邻居的点，这在10x数据中是不能使用的，更改了这个，上面的all_nicheType就要更改
      {if (index_10x) filter(., N_bin50_raw >= min_neighbor_num) else .} %>%
      group_by(label,celltype,neighbors) %>% # 同时考虑了label和celltype和neighbors类型
      summarise(N_bin50_raw = first(N_bin50_raw),f_ij = n()) %>% # 每个区域中每种niche数量：中心+周围均不同
      group_by(label) %>% # orig.ident
      summarise(
        # 下面这些值都是每个区域的值
        N_bin50_raw = first(N_bin50_raw),
        N_niche = n(),
        N_bin50 = sum(f_ij), # 是一个区域所有的点数，因为是sum,不一定是一个区域所有状态的数量
        entropy = -sum(f_ij/N_bin50 * log2(f_ij/N_bin50)), # 每一行的熵，不同中心点+周围细胞类型
        entropy_adj_old = entropy/log2(N_bin50), # 标准化熵或归一化熵，除以已经出现的种类数？有和意义？都矫正了我们需要的差异？处于区间 [0, 1]
        expected_nicheTye = expected_covered_species(all_nicheType, N_bin50), # 期望的niche类型
        entropy_adj = entropy/log2(expected_nicheTye) # log2(expected_nicheTye)是最大的熵
        # 常见的标准化方法是将熵除以其最大可能值
        # 不同分布就用log2(种类可能数)，同一种分布，就用log2(点的总数比较)
      )
    small_v_nm <- entropy_results %>% filter(N_bin50<min_spot_num) %>% pull(label)
    write.table(entropy_results, file = paste0(output_dir,"/",sap_nm,"_entropy_results.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(entropy_results %>% filter(!(label %in% small_v_nm)),
                file = paste0(output_dir,"/",sap_nm,"_entropy_results_fil.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    ### 4.圆环数据 ####
    ## 4.1 celltype_percent ####
    print(paste0("step5: plot pie chart "))
    data <- data %>%
      group_by(label) %>%
      mutate(label_num = n())
    label_df <- data[c(20,22)] %>%
      unique()
    celltype_percent <- data %>%
      select(label,celltype) %>%
      group_by(label,celltype) %>%
      mutate(celltype = str_to_title(celltype)) %>%
      summarise(
        type_counts = n()
      ) %>%
      group_by(label) %>%
      mutate(celltype_ratio = type_counts/sum(type_counts)) %>%
      ungroup() %>%
      select(-type_counts) %>%
      tidyr::pivot_wider(names_from = celltype, values_from = celltype_ratio) %>%
      arrange(label) %>%
      mutate(entropy_adj = entropy_results$entropy_adj) %>%
      mutate(across(where(is.numeric), ~ifelse(is.na(.), 0, .)))

    celltype_percent <- celltype_percent %>%
      mutate(row = center_df_new$row[match(label, center_df_new$label)],
             col = center_df_new$col[match(label, center_df_new$label)])
    write.table(celltype_percent, file = paste0(output_dir,"/",sap_nm,"_celltype_percent.txt"),
                sep = "\t", col.names = T, row.names = F, quote =F)

    celltype_percent <- celltype_percent %>% # 不影响上面的保存
      filter(!(label %in% small_v_nm))
    write.table(celltype_percent, file = paste0(output_dir,"/",sap_nm,"_celltype_percent_fil.txt"),
                sep = "\t", col.names = T, row.names = F, quote =F)

    ## 4.2 celltype_percent_long ####
    # 使用的是celltype_percent_fil
    celltype_percent_long <- celltype_percent %>%
      pivot_longer(cols = starts_with("Cortex") | starts_with("Junction") |
                     starts_with("Medulla") | starts_with("Others"),
                   names_to = "type",
                   values_to = "percent") %>%
      mutate(type = gsub("of_", "", type)) %>%
      mutate(label_num = label_df$label_num[match(label, label_df$label)])
    plot_entropy_data <- distinct(celltype_percent_long, row,col,label, entropy_adj,label_num)

    ### 5.轮廓数据和图 ####
    ## 5.1 segment_df ####
    if(index_10x){
      print(paste0("Step6: 10X skip plot silhouette"))
    }else{
      print(paste0("step6: segment_df"))
      mat <- Matrix::sparseMatrix(
        i = data$row,
        j = data$col,
        x = as.numeric(factor(data$label)))
      dim_mat <- dim(mat)

      segment_df <- map2(data$row, data$col, ~{
        x <- .x
        y <- .y

        res <- tibble()
        if(x+1 <= dim_mat [1] && mat[x + 1,y] > 0 && mat[x, y] != mat[x + 1,y]){
          res <- res %>%
            bind_rows(tibble(x1 = x, y1 = y, x2 = x + 1, y2 = y)) # 返回的x代表的是row，y代表的是col
          index1 <- 1
        } else {
          index1 <- 0
        }
        if(y+1 <= dim_mat [2] && mat[x, y + 1] >0 && mat[x, y] != mat[x, y + 1]){
          res <- res %>%
            bind_rows(tibble(x1 = x, y1 = y, x2 = x, y2 = y + 1))
          index2 <- 1
        }else{
          index2 <- 0
        }

        if(index1 + index2 == 2 && x >1 && y>1 && mat[x, y] != mat[x - 1,y] && mat[x, y] != mat[x,y-1] ){
          res <- res[-(c(-1,0)+nrow(res)),]
        }

        res
      }) %>%
        list_rbind() %>%
        distinct()
      write.table(segment_df, file = paste0(output_dir,"/",sap_nm,"_segment.txt"),
                  sep = "\t", col.names = T, row.names = F, quote =F)
    }
  }

  ######## only_plot ############
  if(only_plot){
    center_df_new <- read.table(file = paste0(input_txt_dir,"/",sap_nm,"_center_df.txt"),
                                header = T, sep = "\t")
    data <- read.table(file = paste0(input_txt_dir,"/",sap_nm,"_data.txt"),
                       header = T, sep = "\t")
    data <- data %>%
      group_by(label) %>%
      mutate(label_num = n())
    label_df <- data[c(20,22)] %>%
      unique()
    celltype_percent <- read.table(file = paste0(input_txt_dir,"/",sap_nm,"_celltype_percent_fil.txt"),
                                   header = T, sep = "\t")
    if(!index_10x){
      segment_df <- read.table(file = paste0(input_txt_dir,"/",sap_nm,"_segment.txt"),
                               header = T, sep = "\t")
    }

    celltype_percent_long <- celltype_percent %>%
      pivot_longer(cols = starts_with("Cortex") | starts_with("Junction") |
                     starts_with("Medulla") | starts_with("Others"),
                   names_to = "type",
                   values_to = "percent") %>%
      mutate(type = gsub("of_", "", type)) %>%
      mutate(label_num = label_df$label_num[match(label, label_df$label)])
    plot_entropy_data <- distinct(celltype_percent_long, row,col,label, entropy_adj,label_num)
  }

  ## 5.2 轮廓图 ####
  if(index_10x){
    print(paste0("Step6: 10X skip plot silhouette"))
  }else{
    print(paste0("step6: plot silhouette "))

    p <- ggplot() +
      geom_segment( # 区域线
        data = segment_df,
        # aes(x = y1-0.25, y = x1-0.25, xend = y2+0.25, yend = x2+0.25),
        aes(x = y1, y = x1, xend = y2, yend = x2),
        alpha = 1,
        color = "black",
        linewidth = 2
      )+
      coord_fixed()+  # 保持坐标轴比例一致
      theme_void() +
      scale_y_reverse()

    ggsave(plot = p,filename = paste0(photo_dir,"/",sap_nm,"_segment.png"),
           width = 10, height = 10,bg = "transparent",dpi = 900)
    pdf(paste0(photo_dir,"/",sap_nm,"_segment.pdf"),
        width = 10, height = 10)
    print(p)
    dev.off()
  }

  ### 6.比例饼图 ####
  print(paste0("step7: plot ratio ",sap_nm))
  tag_ratio <- table(data$celltype) %>%
    as.data.frame() %>%
    mutate(Var1 = str_to_title(Var1)) %>%
    mutate(Var1 = factor(Var1, levels = c("Others", "Medulla", "Junction", "Cortex")))
  tag_ratio$percent <- tag_ratio$Freq/sum(tag_ratio$Freq)
  write.table(tag_ratio, file = paste0(output_dir,"/",sap_nm,"_ratio.txt"),
              sep = "\t", quote = F, row.names = F)

  pdf(paste0(photo_dir,"/",sap_nm,"_ratio.pdf"),
      width = 3, height = 2.5)
  print(
    ggplot(tag_ratio,aes(x=1,y=percent,fill=Var1))+
      geom_bar(width = 1,stat = "identity")+
      coord_polar("y",start = 0,direction = -1)+
      theme_void()+
      theme(legend.position = "right")+
      scale_fill_manual(values = values)+
      labs(title = sap_nm,fill ="Group")+
      # geom_text(aes(label = scales::percent(pct)),position = position_stack(vjust = 0.5))+
      # 圆的外边添加标签
      geom_text(aes(x = 1.7,label = scales::percent(percent)),
                position = position_stack(vjust = 0.5),
                size = 2)+
      theme(
        plot.title = element_text(size = 8,face = 'bold',hjust = 0.5),
        legend.title = element_text(size = 6,face = 'bold',hjust = 0),
        legend.text = element_text(size = 5,face = 'bold',hjust = 0)
      )
  )
  dev.off()

  ### 7.分区和分组图 ####
  print(paste0("step8: plot partition visualization "))
  ## 7.1 分区可视化
  p1 <- ggplot(data, aes(x = col, y = row, color = label)) +
    geom_point(
      size = p1_size, # bin50,相比于绘图点要大一点才行，否则没法重叠，是被区域边界需要
      shape = 16
    ) +
    # scale_color_manual() + # 根据细胞类型自定义颜色
    coord_fixed(ratio = 1) + # 保持坐标轴比例一致
    theme_bw() +
    labs(x = "Column", y = "Row", color = "Cell Type") +
    scale_y_reverse()+
    guides(color = guide_legend(override.aes = list(size = 4)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  pdf(paste0(photo_dir,"/",sap_nm,"_region.pdf"),
      width = 8, height = 8)
  print(p1)
  dev.off()

  ## 7.2 分类可视化-ggplot
  plot_grp_data <- data.frame(
    row = data$row, # row和col可能为0
    col = data$col,
    group = data$celltype
  ) %>%
    mutate(group = str_to_title(group))  # 将首字母

  p2 <- ggplot() +
    geom_point(data = plot_grp_data, aes(x = col, y = row, color = group),
               size = p2_size, # pub2/4/8
               shape = 16) +
    theme_void() +
    coord_fixed()+  # 保持坐标轴比例一致
    scale_color_manual(values = values) +
    scale_y_reverse()+
    labs(x = "Column", y = "Row", color = "Cell Type") +
    guides(color = "none")+
    coord_fixed()+  # 保持坐标轴比例一致
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.background = element_rect(fill = "black")
    )
  pdf(paste0(photo_dir,"/",sap_nm,"_grp.pdf"),
      width = 8, height = 8)
  print(p2)
  dev.off()

  #
  p2_alpha <- ggplot() +
    geom_point(data = plot_grp_data, aes(x = col, y = row, color = group),
               size = p2_size, # pub2/4/8
               shape = 16) +
    theme_void() +
    coord_fixed()+  # 保持坐标轴比例一致
    scale_color_manual(values = alpha(values,alpha_index)) +
    scale_y_reverse()+
    labs(x = "Column", y = "Row", color = "Cell Type") +
    guides(color = "none")+
    coord_fixed()+  # 保持坐标轴比例一致
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.background = element_rect(fill = "black")
    )

  ### 8.合并图 ####
  print(paste0("step9: combine partition and group visualization "))
  # p1_merge: 皮髓分组图
  # p2_merge: 加饼图: 变动的值
  # p3_merge: 加熵图：变动的值
  # p4_merge: 加轮廓图
  p1_merge <- p2_alpha # 分组图
  p2_merge <- p1_merge +
    new_scale_fill()+ # 饼图
    geom_arc_bar(
      data = celltype_percent_long,
      aes(x0 = col, y0 = row,
          r0 = 0,
          r = scales::rescale(label_num/mean(label_num),to = c(0.8, 1.4))*p2_merge_r,
          amount = percent,,
          fill = type),
      # color = "grey90",
      alpha = 1,
      size = 0.5, # linewidth, 旧版本没法识别？
      stat = "pie"
    )+
    coord_fixed()+  # 保持坐标轴比例一致
    scale_fill_manual(values = values) +
    labs(fill = "Cell Type")

  # entropy and no_seg
  p3_merge <- p2_merge+
    new_scale_fill()+ # 熵的饼图
    geom_point(data = plot_entropy_data,  # 确保每个 label 只绘制一个点
               aes(x = col, y = row, fill = entropy_adj),
               size = scales::rescale(plot_entropy_data$label_num/mean(plot_entropy_data$label_num),c(0.7,1.5))*p3_merge_size,
               shape = 21,
               stroke = 1.2,
               color = "black") +  # 形状为圆形，并设置大小和内部填充颜色
    scale_fill_viridis_c(option = "H",
                         limits = c(0, 1) # bin200/bin50_rawgem，pub2/4/7/8
                         # limits = c(0, 0.8) # bin50
    ) +
    # scale_fill_gradientn(colors = c('black',"blue",'skyblue', "yellow",'orange', "red",'darkred')
    #                      ,limits = c(0, 1))+
    theme_void()+
    labs(fill = "Adjusted entropy")+
    guides(size = "none")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.background = element_rect(fill = "black"), # 注释掉就是白底黑字
          # plot.background = element_rect(fill = "black"),
          # legend.text = element_text(color = "white"),
          # legend.title = element_text(color = "white")
    )
  pdf(paste0(photo_dir,"/",sap_nm,"_no_seg.pdf"),
      width = 9, height = 8)
  print(p3_merge)
  dev.off()


  # merge all
  if(index_10x){
    print("10X no merge1 and merge2")
  }else{
    # merge1
    p4_merge <- p3_merge+
      geom_segment( # 区域线
        data = segment_df,
        # aes(x = y1-0.25, y = x1-0.25, xend = y2+0.25, yend = x2+0.25),
        aes(x = y1, y = x1, xend = y2, yend = x2),
        alpha = 1,
        color = "black",
        linewidth = 2
      )
    pdf(paste0(photo_dir,"/",sap_nm,"_merge1.pdf"),
        width = 9, height = 8)
    print(p4_merge)
    dev.off()

    # merge2
    pdf(paste0(photo_dir,"/",sap_nm,"_merge2.pdf"),
        width = 16, height = 8)
    print(p1_merge + p4_merge)
    dev.off()
  }

  # merge_black
  if(!index_10x){
    p3_merge_black <- p2_merge+
      new_scale_fill()+ # 熵
      geom_point(data = plot_entropy_data,  # 确保每个 label 只绘制一个点
                 aes(x = col, y = row, fill = entropy_adj),
                 size =  scales::rescale(plot_entropy_data$label_num/mean(plot_entropy_data$label_num),to = c(0.6,1.5))*p3_merge_size,
                 shape = 21,
                 stroke = 1.5,
                 color = "black") +  # 形状为圆形，并设置大小和内部填充颜色
      scale_fill_viridis_c(option = "H",
                           limits = c(0, 1) # bin200/bin50_rawgem，pub2/4/7/8
                           # limits = c(0, 0.8) # bin50
      ) +
      # scale_fill_gradientn(colors = c('black',"blue",'skyblue', "yellow",'orange', "red",'darkred')
      #                      ,limits = c(0, 1))+
      geom_segment( # 区域线
        data = segment_df,
        # aes(x = y1-0.25, y = x1-0.25, xend = y2+0.25, yend = x2+0.25),
        aes(x = y1, y = x1, xend = y2, yend = x2),
        alpha = 1,
        color = "grey95",
        linewidth = 2
      ) +
      theme_void()+
      labs(fill = "Adjusted entropy")+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "black"), # 注释掉就是白底黑字
            plot.background = element_rect(fill = "black"),
            legend.text = element_text(color = "white"),
            legend.title = element_text(color = "white")
      )
    pdf(paste0(photo_dir,"/",sap_nm,"_merge1_black.pdf"),
        width = 9, height = 8)
    print(p3_merge_black)
    dev.off()
  }

  # 10X: no_grp
  if(index_10x){
    p_no_grp <- ggplot() +
      geom_arc_bar(
        data = celltype_percent_long, # 饼图
        aes(x0 = col, y0 = row,
            r0 = 0,
            r = (label_num/mean(label_num))*p2_merge_r,
            amount = percent,,
            fill = type),
        # color = "grey90",
        alpha = 1,
        size = 0.25, # linewidth, 旧版本没法识别？
        stat = "pie"
      )+
      coord_fixed()+  # 保持坐标轴比例一致
      scale_fill_manual(values = values) +
      labs(fill = "Cell Type") +
      new_scale_fill()+ # 熵
      scale_y_reverse()+
      geom_point(data = plot_entropy_data,  # 确保每个 label 只绘制一个点
                 aes(x = col, y = row, fill = entropy_adj),
                 size = scales::rescale(plot_entropy_data$label_num/mean(plot_entropy_data$label_num),to = c(0.6,1.5))*p3_merge_size,
                 shape = 21,
                 stroke = 1.5,
                 color = "black") +  # 形状为圆形，并设置大小和内部填充颜色
      scale_fill_viridis_c(option = "H",
                           limits = c(0, 1) # bin200/bin50_rawgem，pub2/4/7/8
                           # limits = c(0, 0.8) # bin50
      ) +
      # scale_fill_gradientn(colors = c('black',"blue",'skyblue', "yellow",'orange', "red",'darkred')
      #                      ,limits = c(0, 1))+
      theme_void()+
      labs(fill = "Adjusted entropy")+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.background = element_rect(fill = "black"), # 注释掉就是白底黑字
            # plot.background = element_rect(fill = "black"),
            # legend.text = element_text(color = "white"),
            # legend.title = element_text(color = "white")
      )
    pdf(paste0(photo_dir,"/",sap_nm,"_no_grp.pdf"),
        width = 9, height = 8)
    print(p_no_grp)
    dev.off()
  }

  print("All done!")

}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalNicheGeneTrend
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate Gene Expression Trends Across Niche Regions
#'
#' Placeholder function for analyzing gene expression patterns across
#' different spatial niches. Currently under development.
#'
#' @param x Placeholder parameter for future expansion
#'
#' @return NULL (invisible)
#'
#' @details
#' This function is designed to analyze how gene expression changes across
#' the spatial organization hierarchy identified by CalNicheOSE.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Under development
#' CalNicheGeneTrend()
#' }
CalNicheGeneTrend <- function(x){

}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalNichePseudo
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Generate Pseudo-time Analysis for Niche Development
#'
#' Placeholder function for performing pseudo-time analysis on spatial niches
#' to infer developmental trajectories. Currently under development.
#'
#' @param x Placeholder parameter for future expansion
#'
#' @return NULL (invisible)
#'
#' @details
#' This function is designed to order spatial niches along a developmental
#' trajectory based on gene expression patterns and spatial organization.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Under development
#' CalNichePseudo()
#' }
CalNichePseudo <- function(x){

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalNichePredict
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Predict Niche Types or Properties
#'
#' Placeholder function for predictive modeling of niche properties based on
#' spatial and molecular features. Currently under development.
#'
#' @param x Placeholder parameter for future expansion
#'
#' @return NULL (invisible)
#'
#' @details
#' This function is designed to build predictive models for niche classification
#' or property estimation using machine learning approaches.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Under development
#' CalNichePredict()
#' }
CalNichePredict <- function(x){

}




