

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalSampPathoTrack
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Calculate pathogen tracking across samples
#'
#' Analyzes pathogen distribution across temporal or comparative samples,
#' infers dissemination networks, and visualizes pathogen spread patterns.
#'
#' @param STID_obj A STID object
#' @param loop_id Character, multi-sample analysis identifier
#' @param samp_grp_index Logical, whether to group by sample groups (default: FALSE)
#' @param meta_key Character, metadata key for retrieving cell data
#' @param niche_key Character, niche key for filtering niche cells (optional)
#' @param group_by Character, column name for cell type grouping
#' @param pos_colnm Character, column name containing positive spot labels
#' @param neg_value Character, value indicating negative spots (default: "neg")
#' @param col Character vector, color palette for cell types
#' @param return_data Logical, whether to return results (default: FALSE)
#' @param grp_nm Character, group name for output (default: NULL, uses timestamp)
#' @param dir_nm Character, directory name for output (default: "M4_CalSampPathoTrack")
#'
#' @return If return_data = TRUE, returns a list containing:
#'   \itemize{
#'     \item pathogen_data: Pathogen-positive cell counts and ratios per time point
#'     \item pathogen_wide: Wide format pathogen ratio matrix
#'     \item region_summary: Peak time and total burden per cell type
#'     \item edges_df: Inferred dissemination network edges
#'     \item spatiotemporal_edges: Time-resolved transition edges
#'   }
#'
#' @details
#' The function performs three types of analyses:
#' \enumerate{
#'   \item Temporal trend analysis of pathogen-positive cell type composition
#'   \item Peak-based dissemination network inference (which cell types seed others)
#'   \item Transition-based spatiotemporal mapping of pathogen spread
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Track pathogen spread across time points
#' results <- CalSampPathoTrack(
#'   STID_obj = stid_object,
#'   loop_id = "infection_series",
#'   group_by = "cell_type",
#'   pos_colnm = "Label_pathogen"
#' )
#' }
CalSampPathoTrack <- function(STID_obj = NULL, loop_id = NULL, # must support
                              samp_grp_index = FALSE,
                              meta_key = NULL,
                              niche_key = NULL, # only support one value
                              group_by = NULL,
                              pos_colnm = NULL, neg_value = "neg",
                              col = COLOR_LIST[["PALETTE_WHITE_BG"]],
                              return_data = FALSE,
                              grp_nm = NULL,dir_nm = "M4_CalSampPathoTrack"){

  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(loop_id, pos_colnm, col)
  if(is.null(group_by)){
    clog_warn("group_by is NULL, will use the default celltype_colnm from STID object.")
    group_by <- GetInfo(STID_obj, info_key = "data_info",sub_key = "celltype_colnm")[[1]]
  }
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  results_list <- list()
  clog_normal("Execute multi-sample analysis...")
  loop_multi <- .check_loop_multi(STID_obj = STID_obj, loop_id = loop_id)
  if(length(loop_multi) >1){
    clog_error("Multiple loop_id detected in STID object. Please specify one loop_id for analysis.")
  }

  #> main pipeline
  for(i in seq_along(loop_multi)){
    i_multi <- loop_multi[i]
    clog_loop(paste0("Processing multi_id: ", i_multi, " (", i, "/", length(loop_multi), ")"))
    samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
    samp_id2grp <- data.frame(
      samp_id = STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@samp_info$samp_id,
      samp_grp = STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@samp_info$samp_grp
    )

    #>
    if(!is.null(niche_key)){
      logic_niche <- .check_niche_exist(STID_obj, i_multi, niche_key, samp_mode  = "MS")
      if(!logic_niche){
        clog_warn(paste0("niche_key: ", niche_key, " not found in MultiSampNiche of multi_id: ", i_multi, ", skipping..."))
        next
      }
      clog_normal("Using niche_key for analysis...")
      Niche_cells <- GetMSNicheCells(STID_obj = STID_obj, niche_key = niche_key,loop_id = i_multi)[[1]]
    }else{
      if(is.null(meta_key)){
        clog_warn("Both niche_key and meta_key are NULL, will use meta_key: coord for MS niche plotting.")
        meta_key <- "coord"
      }
      clog_normal("Using meta_key for analysis...")
      Niche_cells <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]] %>%
        filter(!!sym(samp_colnm) %in% samp_id2grp$samp_id) # actually all, not niche
    }
    .check_column_exist(Niche_cells, group_by)
    if(samp_grp_index){
      Niche_cells <- Niche_cells %>%
        mutate(samp_grp = samp_id2grp$samp_grp[match(Niche_cells[[samp_colnm]],samp_id2grp$samp_id)],
               .after = all_of(samp_colnm))
      new_samp_colnm <- "samp_grp"
    }else{
      new_samp_colnm <- samp_colnm
    }
    if(!is.null(niche_key)){
      Niche_cells <- Niche_cells %>%
        mutate(Niche_label = ifelse(is_Niche,"Niche","Bystander"), .after = "is_Niche") %>%
        mutate(Niche_label = factor(Niche_label,levels = c("Bystander","Niche")))
    }else{
      Niche_cells <- Niche_cells %>%
        mutate(Niche_label = "All",.after = all_of(group_by)) %>%
        mutate(Niche_label = factor(Niche_label,levels = c("All")))
    }

    #> pos_meta_data
    .check_column_exist(Niche_cells, pos_colnm)
    pos_meta_data <- Niche_cells %>%
      filter(.,!!sym(pos_colnm) != neg_value)
    len_pos_value <- nrow(pos_meta_data)
    cells_pct <- len_pos_value/nrow(Niche_cells)
    clog_normal(paste0("Positive cells percentage: ",round(cells_pct*100,2),"% (",
                       len_pos_value,"/", nrow(Niche_cells),")"))

    #> pathogen_data
    clog_normal("Calculating pathogen load trends...")
    pathogen_data <- pos_meta_data %>%
      group_by(!!sym(group_by),!!sym(new_samp_colnm)) %>%
      mutate(!!sym(new_samp_colnm) := factor(!!sym(new_samp_colnm), levels = samp_id2grp[,samp_grp_index+1])) %>%
      summarise(pos_count = n()) %>%
      arrange(!!sym(new_samp_colnm),!!sym(group_by)) %>%
      group_by(!!sym(new_samp_colnm)) %>%
      mutate(pos_ratio = pos_count/sum(pos_count)) %>%
      ungroup() %>%
      complete(!!sym(group_by), !!sym(new_samp_colnm),
               fill = list(pos_count = 0, pos_ratio = 0)) %>%
      mutate(time = as.integer(!!sym(new_samp_colnm)))
    pathogen_wide <- pathogen_data %>%
      dplyr::select(-pos_count,-time) %>%
      pivot_wider(names_from = !!sym(new_samp_colnm),
                  values_from = pos_ratio, values_fill = 0) %>%
      column_to_rownames(group_by)
    time_points <- seq_along(samp_id2grp[,samp_grp_index+1])
    names(time_points) <- samp_id2grp[,samp_grp_index+1]
    results_list[[i_multi]]$data[["pathogen_data"]] <- pathogen_data
    results_list[[i_multi]]$data[["pathogen_wide"]] <- pathogen_wide

    #> 方法一
    clog_step("Calculating pathogen-positive site composition trends...")
    p1 <- pathogen_data %>%
      ggplot(aes(x = !!sym(new_samp_colnm), y = pos_ratio,
                 color = !!sym(group_by), group = !!sym(group_by))) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      labs(title = "Pathogen-Positive Site Composition",
           x = "Infection Stage",
           y = "Proportion Infected") +
      theme_bw() +
      scale_color_manual(values = col) +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
    print(p1)
    ggsave(plot = p1, filename = paste0(photo_dir,"/",i_multi,"_pathogen_trend.pdf"),
           width = 6, height = 4)

    #> 方法二，效果一般
    clog_step("Inferring peak-based dissemination network...")
    region_summary <- pathogen_data %>%  # 每个细胞的 peak time 和 total burden
      group_by(!!sym(group_by)) %>%
      summarise(
        total_burden = sum(pos_ratio),
        peak_time = time[which.max(pos_ratio)],
        peak_load = max(pos_ratio),
        .groups = "drop"
      )
    results_list[[i_multi]]$data[["region_summary"]] <- region_summary
    edges_df <- tibble(from = character(), to = character(), score = numeric())
    for (i in seq_len(nrow(region_summary))) { # 每种细胞类型只有一个峰值
      A <- region_summary[[group_by]][i]
      tA <- region_summary$peak_time[i]
      loadA_peak <- region_summary$peak_load[i]
      if (loadA_peak <= 0) next

      for (j in seq_len(nrow(region_summary))) { # 不同细胞类型峰间的传播
        B <- region_summary[[group_by]][j]
        if (A == B) next
        tB <- region_summary$peak_time[j]
        if (tB <= tA) next # 只考虑 B 峰时间晚于 A峰时间 的情况

        # 找到 A 峰值之后的第一个时间点
        later_times <- time_points[time_points > tA]
        if (length(later_times) == 0) next
        t_next <- min(later_times)

        # 获取 B 在 tA 和 t_next 的病原载量，# 只考虑 B 峰值晚于 A 的情况
        load_B_at_tA   <- pathogen_wide[B, names(time_points)[time_points == tA]]
        load_B_at_next <- pathogen_wide[B, names(time_points)[time_points == t_next]]
        delta_B <- load_B_at_next - load_B_at_tA
        if (delta_B <= 0) next  #

        # 公式等同于方法三，但是细胞对有限，只看两个时间点峰细胞的乘法
        score <- delta_B * loadA_peak   # 增加的比例分数，最大值为1，不用归一化
        edges_df <- bind_rows(edges_df, tibble(from = A, to = B, score = score))
      }
    }
    edges_df <- edges_df %>%
      group_by(from,to) %>%
      summarise(score = sum(score)) %>%
      arrange(desc(score))
    results_list[[i_multi]]$data[["edges_df"]] <- edges_df

    #> plot
    g <- graph_from_data_frame(edges_df, directed = TRUE, vertices = region_summary)
    p2 <- ggraph(g, layout = 'fr') +
      geom_edge_fan(
        aes(edge_width = score, edge_alpha = score),
        arrow = arrow(length = unit(2, 'mm')),
        end_cap = circle(4, 'mm')
      ) +
      geom_node_point(
        aes(size = total_burden, color = factor(peak_time)),
        show.legend = TRUE
      ) +
      geom_node_text(
        aes(label = name),
        repel = TRUE,
        size = 3,
        family = "sans"
      ) +
      scale_edge_width(range = c(0.5, 1)) +
      scale_edge_alpha(range = c(0.1, 0.9), guide = "none") +
      scale_size_continuous(range = c(5, 10), breaks = breaks_extended(n = 3)) +
      # scale_color_manual(values = c("red", "orange", "skyblue")) +
      labs(
        title = "Inferred Pathogen Dissemination Network",
        color = "Peak Time (dpi)",
        size = "Pathogen Burden"
      ) +
      theme_graph() +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
    print(p2)
    ggsave(p2, filename = paste0(photo_dir,"/",i_multi,"_dissemination_network.svg"), # pdf会报错
           width = 6, height = 5)

    #> 方法三
    clog_step("Inferring transition-based dissemination map...")
    pathogen_wide <- pathogen_wide %>%
      mutate(!!sym(group_by) := rownames(pathogen_wide),.before = 1)
    regions <- rownames(pathogen_wide)
    edges_list <- tibble() # 连续时间点（t → t+1）转移概率矩阵
    for (i in 1:(length(time_points) - 1)) { # t细胞类型
      t_curr <- names(time_points)[i]
      t_next <- names(time_points)[i + 1]
      merged <- pathogen_wide %>% dplyr::select(!!sym(group_by), all_of(c(t_curr,t_next)))
      for (j in seq_along(regions)) {  # 对每一对源-目标区域
        src <- regions[j]
        src_load <- merged[[t_curr]][merged[[group_by]] == src]
        if(src_load <= 0) next

        for( k in seq_along(regions)) { # t+1细胞类型
          tgt <- regions[k]
          if (src == tgt) next
          tgt_load_curr <- merged[[t_curr]][merged[[group_by]] == tgt]
          tgt_load_next <- merged[[t_next]][merged[[group_by]] == tgt]
          delta <- tgt_load_next - tgt_load_curr
          if (delta <= 0) next  # 目标未上升，不考虑传播

          # 传播得分 = 源载量 × 目标增量，不区分源细胞是否上升下降，只看源细胞占比
          score <- src_load * delta
          y_src <- which(regions == src)
          y_tgt <- which(regions == tgt)

          edges_list <- bind_rows(edges_list, tibble(
            from_region = src,
            to_region = tgt,
            time_from = t_curr,
            time_to = t_next,
            score = score,
            x_start = i,
            y_start = y_src,
            x_end = i + 1,
            y_end = y_tgt
          ))
        }
      }
    }
    n_regions <- length(regions)
    edges_list <- edges_list %>%
      mutate(
        y_start = n_regions - y_start + 1,
        y_end = n_regions - y_end + 1
      )
    results_list[[i_multi]]$data[["spatiotemporal_edges"]] <- edges_list

    #> plot
    pathogen_data <- pathogen_data %>%
      mutate(!!sym(group_by) := factor(!!sym(group_by), levels = rev(regions)))
    p3 <- ggplot(data =  pathogen_data,
                 aes(x = !!sym(new_samp_colnm), y = !!sym(group_by), size = pos_ratio, color = pos_ratio),) +
      geom_segment(
        data =  edges_list,
        aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
            alpha = score, linewidth  = score,color = score),
        arrow = arrow(length = unit(0.2, "mm")),
        lineend = "round",
        inherit.aes = FALSE
      ) +
      geom_point(
        shape = 21,
        fill = "grey95",color = "grey30",
        stroke = 0.75
      ) +
      scale_size_continuous(range = c(1, 7),breaks = breaks_extended(n = 3)) +
      scale_linewidth_continuous(range = c(0.5, 2),breaks = breaks_extended(n = 3)) +
      scale_alpha_continuous(range = c(0.2, 0.9), guide = "none") +
      scale_color_gradientn(
        colors = viridis::viridis(100, option = "H")[5:95],
        guide = "colorbar"
      ) +
      scale_x_discrete(expand = c(0,0.5)) +
      scale_y_discrete(expand = c(0,0.5)) +
      labs(
        title = "Transition-based dissemination map",
        x = "Days Post Infection",
        y = "Types",
        linewidth = "Transmission Score",
        color = "Transmission Score"
      ) +
      theme_test() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 13, face = "bold")
      )
    print(p3)
    ggsave(p3, filename = paste0(photo_dir,"/",i_multi,"_spatiotemporal_dissemination_network.pdf"),
           width = 6, height = 6)
  }
  clog_normal("Saving results to RDS file...")
  saveRDS(results_list, file = paste0(output_dir,"/CalSampPathoTrack_Results_List.rds"))

  # >>> Final
  .save_function_params("CalSampPathoTrack", envir = environment(),
                        file = paste0(output_dir,"Log_function_params_(CalSampPathoTrack).log"))
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(CalSampPathoTrack).log")) %>% invisible()
  if(return_data){
    return(results_list)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalSampOSE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Calculate Spatial Organization Entropy (OSE) for samples
#'
#' Performs spatial organization entropy analysis to quantify tissue spatial
#' heterogeneity based on cell type distribution patterns. This function implements
#' the OSE algorithm which partitions tissue into regions and calculates entropy
#' based on local cell type diversity.
#'
#' @param STID_obj A STID object containing spatial transcriptomics data
#' @param loop_id Character, multi-sample analysis identifier
#' @param meta_key Character, metadata key for retrieving cell data
#' @param group_by Character, column name for cell type grouping
#' @param OSE_dist_m Numeric, distance weight (default: 0.2). Smaller values give
#'        more weight to expression distance in region segmentation
#' @param OSE_PCA_nPC Integer, number of PCs for PCA dimensionality reduction
#'        (default: 15)
#' @param OSE_window Numeric, window size for superpixel generation. If NULL,
#'        automatically set to max(x_range, y_range)/10
#' @param OSE_minSpotNum Integer, minimum spots per region. If NULL, automatically
#'        set to nrow(meta_data)/1000
#' @param only_plot Logical, whether to only regenerate plots from existing results
#'        (default: FALSE)
#' @param plot_params List of plotting parameters:
#'        \itemize{
#'          \item p3_size: Point size for partition plot
#'          \item p4_size: Point size for group plot
#'          \item p2_merge_r: Radius for pie chart in merged plot
#'          \item p3_merge_size: Point size for entropy in merged plot
#'        }
#' @param col Character vector, color palette for cell types
#' @param return_data Logical, whether to return results as a list (default: FALSE)
#' @param grp_nm Character, group name for output directory (default: NULL, uses timestamp)
#' @param dir_nm Character, directory name for output (default: "M4_CalSampOSE")
#'
#' @return If return_data = TRUE, returns a list containing:
#'   \itemize{
#'     \item center_df: Data frame of region centers
#'     \item data: Cell-level OSE data with region assignments
#'     \item entropy_results: Region entropy calculations
#'     \item celltype_percent: Cell type proportions per region
#'     \item segment_df: Region boundary segments
#'   }
#'
#' @importFrom scales rescale
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- CalSampOSE(
#'   STID_obj = stid_object,
#'   loop_id = "comparison_1",
#'   meta_key = "M2_NicheDetect_STS_20240101",
#'   group_by = "cell_type",
#'   OSE_window = 50,
#'   OSE_minSpotNum = 100
#' )
#' }
CalSampOSE <- function(STID_obj = NULL,
                       loop_id = NULL, # must support one
                       meta_key = NULL,
                       group_by = NULL,
                       OSE_dist_m = 0.2,  # 距离权重，越小表达距离的权重越大，建议在0.001-0.01范围内
                       OSE_PCA_nPC = 15,
                       OSE_window = NULL,
                       OSE_minSpotNum = NULL,
                       only_plot = FALSE,
                       plot_params = list(
                         p3_size = 2,
                         p4_size = 0.3, # 0.42
                         p2_merge_r = 15,
                         p3_merge_size = 9
                       ),
                       col = COLOR_LIST[["PALETTE_WHITE_BG"]],
                       return_data = FALSE,
                       grp_nm = NULL,dir_nm = "M4_CalSampOSE"){

  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(loop_id, col)
  if(is.null(group_by)){
    clog_warn("group_by is NULL, will use the default celltype_colnm from STID object.")
    group_by <- GetInfo(STID_obj, info_key = "data_info",sub_key = "celltype_colnm")[[1]]
  }
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  results_list <- list()
  clog_normal("Execute multi-sample analysis...")
  loop_multi <- .check_loop_multi(STID_obj = STID_obj, loop_id = loop_id)
  if(length(loop_multi) >1){
    clog_error("Multiple loop_id detected in STID object. Please specify one loop_id for analysis.")
  }

  #> main pipeline
  for(i in seq_along(loop_multi)){
    i_multi <- loop_multi[i]
    clog_loop(paste0("Processing multi_id: ", i_multi, " (", i, "/", length(loop_multi), ")"))
    samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
    samp_id2grp <- data.frame(
      samp_id = STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@samp_info$samp_id,
      samp_grp = STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@samp_info$samp_grp
    )

    #> no Niche_key, no samp_grp_index
    clog_normal("Using meta_key for analysis...")
    Niche_cells <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]] %>%
      filter(!!sym(samp_colnm) %in% samp_id2grp$samp_id) # actually all, not niche
    .check_column_exist(Niche_cells, group_by)
    new_samp_colnm <- samp_colnm
    Niche_cells <- Niche_cells %>%
      mutate(Niche_label = "All",.after = all_of(group_by)) %>%
      mutate(Niche_label = factor(Niche_label,levels = c("All")))

    #> 单样本遍历
    for(j in seq_along(samp_id2grp$samp_id)){
      j_samp <- samp_id2grp$samp_id[j]
      output_samp_dir <- paste0(output_dir,"/",j_samp)
      photo_samp_dir <- paste0(photo_dir,"/",j_samp)
      dir.create(output_samp_dir, showWarnings = FALSE, recursive = TRUE)
      dir.create(photo_samp_dir, showWarnings = FALSE, recursive = TRUE)
      clog_loop(paste0("Processing ",i_multi,": ", j_samp, " (", j, "/", length(samp_id2grp$samp_id), ")"))
      j_seurat <- subset(STID_obj, subset = !!sym(samp_colnm) == j_samp)
      j_Niche_cells <- Niche_cells %>%
        filter(!!sym(new_samp_colnm) == j_samp)
      j_results <- .OSE_pipeline(STID_obj = j_seurat,
                                samp_nm = j_samp,
                                meta_data = j_Niche_cells,
                                group_by = group_by,
                                OSE_window = OSE_window,
                                OSE_dist_m = OSE_dist_m,
                                OSE_PCA_nPC = OSE_PCA_nPC,
                                OSE_minSpotNum = OSE_minSpotNum,
                                only_plot = only_plot,
                                plot_params = plot_params,
                                col = col,
                                output_dir = output_samp_dir,
                                photo_dir = photo_samp_dir )
      results_list[[i_multi]][[j_samp]]$data <- j_results
    }
  }

  #> save results
  clog_normal("Saving results to RDS file...")
  saveRDS(results_list, file = paste0(output_dir,"/CalSampOSE_Results_List.rds"))

  # >>> Final
  .save_function_params("CalSampOSE", envir = environment(),
                        file = paste0(output_dir,"Log_function_params_(CalSampOSE).log"))
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(CalSampOSE).log")) %>% invisible()
  if(return_data){
    return(results_list)
  }
}


#' Internal OSE pipeline for single sample
#'
#' Performs the core OSE algorithm on a single sample, including PCA,
#' superpixel region segmentation, iterative center update, entropy calculation,
#' and visualization generation.
#'
#' @param STID_obj A Seurat object (converted from STID)
#' @param samp_nm Character, sample name for output files
#' @param meta_data Data frame, cell metadata with coordinates and cell types
#' @param group_by Character, column name for cell type grouping
#' @param OSE_dist_m Numeric, distance weight for region segmentation
#' @param OSE_PCA_nPC Integer, number of PCs for PCA
#' @param OSE_window Numeric, window size for superpixel generation
#' @param OSE_minSpotNum Integer, minimum spots per region
#' @param only_plot Logical, whether to only regenerate plots
#' @param plot_params List of plotting parameters
#' @param col Character vector, color palette
#' @param output_dir Character, output directory for data files
#' @param photo_dir Character, output directory for plots
#'
#' @return List containing OSE results for the sample
#'
#' @details
#' The OSE algorithm works as follows:
#' \enumerate{
#'   \item PCA dimensionality reduction on normalized expression data
#'   \item Initial superpixel grid generation based on spatial coordinates
#'   \item Iterative assignment of spots to nearest regions (weighted by space and expression)
#'   \item Region center update until convergence
#'   \item Spatial entropy calculation based on cell type diversity in local neighborhoods
#'   \item Visualization of region partitions, cell type distributions, and entropy heatmaps
#' }
#'
#' @importFrom scales rescale
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggforce geom_arc_bar
#'
#' @keywords internal
#' @noRd
.OSE_pipeline <- function(STID_obj = NULL,
                         samp_nm = NULL,
                         meta_data = NULL,
                         group_by = NULL,
                         OSE_dist_m = 0.2,  # 距离权重，越小表达距离的权重越大，建议在0.001-0.01范围内
                         OSE_PCA_nPC = 15,
                         OSE_window = NULL,
                         OSE_minSpotNum = NULL,
                         only_plot = FALSE,
                         plot_params = list(
                           p3_size = 2,
                           p4_size = 0.3, # 0.42
                           p2_merge_r = 15,
                           p3_merge_size = 9
                         ),
                         col = COLOR_LIST[["PALETTE_WHITE_BG"]],
                         output_dir = NULL,
                         photo_dir = NULL
){
  results_list <- list()

  #> basic parameter setting
  data_format <- GetInfo(STID_obj, info_key = "data_info", sub_key = "data_format")[[1]]
  data_platform <- GetInfo(STID_obj, info_key = "data_info", sub_key = "data_platform")[[1]] # Visium
  binsize <- GetInfo(STID_obj, info_key = "data_info", sub_key = "binsize")[[1]]
  coord_interval <- GetInfo(STID_obj, info_key = "data_info", sub_key = "coord_interval")[[1]][samp_nm]
  base_unit <- GetInfo(STID_obj, info_key = "data_info", sub_key = "base_unit")[[1]]
  x_max_min <- max(meta_data$x) - min(meta_data$x)
  y_max_min <- max(meta_data$y) - min(meta_data$y)
  if(is.null(OSE_minSpotNum)){
    OSE_minSpotNum <- (nrow(meta_data)/1000) %>% round()
    clog_warn(paste0("OSE_minSpotNum is NULL, will set to ",OSE_minSpotNum))
  }
  if(is.null(OSE_window)){
    OSE_window <- max(x_max_min/10,y_max_min/10) %>% round()
    clog_warn(paste0("OSE_window is NULL, will set to ",OSE_window))
  }
  STID_obj <- as.Seurat(STID_obj)

  #> process the coord
  if(data_platform == "Visium"){
    min_x <- min(STID_obj@meta.data$x)
    min_y <- min(STID_obj@meta.data$y)
    gap_x <- STID_obj@meta.data$x %>% unique() %>% sort() %>% diff() %>% table() %>% names() %>% as.numeric() %>% sort() %>% .[2]
    print(paste0("gap_x: ",gap_x))
    gap_y <- STID_obj@meta.data$y %>% unique() %>% sort() %>% diff() %>% table() %>% names() %>% as.numeric() %>% sort() %>% .[2]
    print(paste0("gap_y: ",gap_y))
    gap_select <- min(gap_x,gap_y)
    print(paste0("gap_x and gap_y ratio: ",gap_x/gap_y))
    STID_obj@meta.data$x <- round((STID_obj@meta.data$x-min_x+1)/gap_select) # 除以gap_select，而不是gap_x，因为x和y比例关系是正确的，不能矫正了
    STID_obj@meta.data$y <- round((STID_obj@meta.data$y-min_y+1)/gap_select)
  }
  meta_data$row <- meta_data$y + coord_interval # 方式出现0的情况，后面构建稀疏矩阵index1 = T会报错
  meta_data$col <- meta_data$x + coord_interval # 这个坐标都是增加，不会影响实际的图形的相对位置

  #> plot parameter setting
  p3_size <- plot_params$p3_size
  p4_size <- plot_params$p4_size
  p2_merge_r <- plot_params$p2_merge_r
  p3_merge_size <- plot_params$p3_merge_size

  #> OSE start
  clog_step(paste0("Start OSE pipeline..."))

  if(!only_plot){
    #> Seurat
    clog_normal("Running PCA...")
    STID_obj <- NormalizeData(STID_obj)
    STID_obj <- FindVariableFeatures(STID_obj)
    STID_obj <- ScaleData(STID_obj)
    STID_obj <- RunPCA(STID_obj,verbose = F,
                       features = VariableFeatures(STID_obj),
                       npcs = OSE_PCA_nPC)

    #> OSE OSE_meta_data
    clog_step("Extract PC OSE_meta_data...")
    pc_data <- STID_obj@reductions[["pca"]]@cell.embeddings %>%
      as.data.frame() %>%
      mutate(cellid = rownames(meta_data))
    OSE_meta_data <- data.frame(
      row = meta_data$row, # row和col可能为0
      col = meta_data$col,
      celltype = meta_data[[group_by]],
      cellid = rownames(meta_data)
    ) %>%
      left_join(pc_data, by = "cellid") %>%
      dplyr::select(1:(OSE_PCA_nPC+4)) %>%  # row + col + celltype + cellid  + PC1-15
      `row.names<-`(.$cellid) %>%
      mutate(label = NA,dist = Inf)
    min_col <- min(OSE_meta_data$col)
    min_row <- min(OSE_meta_data$row)


    #> 生成超像素块区域及中心点
    clog_step(paste0("Generate superpixel blocks and center points..."))
    create_group <- function(x, interval) {
      ceiling((x - min(x)) / interval + 0.0001) # 向上取整，将每个spot分为哪个区域中
    }
    center_df <- OSE_meta_data %>%   # 重新迭代前要将这段跑一下，重新得到初始的center_df
      mutate(
        col_group = create_group(col, OSE_window),  # 先添加后删除
        row_group = create_group(row, OSE_window)
      ) %>%
      group_by(col_group, row_group) %>%
      summarise(
        spot_num = n(),
        across(starts_with("PC_"), mean) # 其余的列都丢弃，每个区域用所有细胞的均值PC值代表
      ) %>%
      ungroup() %>%
      mutate(
        col = as.numeric(col_group) * OSE_window + min_col-0.5*OSE_window, # 手动定义中心点
        row = as.numeric(row_group) * OSE_window + min_row-0.5*OSE_window
      ) %>%
      filter(spot_num > 1) %>% # !!!!!!
      mutate(label = str_c("V", seq(n()))) %>%
      arrange(label) %>%
      dplyr::select(c(-col_group,-row_group)) # 最终将col_group和row_group转为坐标了

    #> 迭代更新区域
    clog_step(paste0("Iteratively update center points..."))
    res <- Inf
    step <- 1
    data_backup <- OSE_meta_data
    res_value <- (nrow(meta_data)/10000) %>% round()
    while (res > res_value) {
      OSE_meta_data <- data_backup

      # 计算每个区域与细胞的距离
      dist_s <- proxy::dist(
        OSE_meta_data %>% dplyr::select(col, row), # OSE_meta_data源代码中为spots
        center_df %>% dplyr::select(col, row), # center_df源代码中为centers
        method = "Euclidean"
      ) %>%
        unclass() %>%
        as.data.frame() %>%
        rownames_to_column("cellid") %>%
        pivot_longer(-cellid, names_to = "label", values_to = "d_s") #%>% # value就是距离
      #filter(d_s <= OSE_window)
      # print(dim(dist_s)) # 每个点距离每个中心点的实际距离

      dist_exp <- proxy::dist(
        OSE_meta_data %>% dplyr::select(starts_with("PC_")),
        center_df %>% dplyr::select(starts_with("PC_")),
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
          dist = sqrt(OSE_dist_m * d_s^2 + (1-OSE_dist_m)*d_e^2) # OSE_dist_m是距离权重，也只有一个OSE_dist_m参数
        ) %>%
        group_by(cellid) %>%
        slice_min(dist) %>% # 选出每个点距离最近的中心点
        dplyr::select(cellid, label, dist)

      OSE_meta_data <- OSE_meta_data %>% # OSE_meta_data在源代码中有时为df
        dplyr::select(-c(label, dist)) %>%
        left_join(dists)

      # 根据新的点更新中心点
      center_df_new <- OSE_meta_data %>%
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

      res <- as.matrix(center_df %>% dplyr::select(-label)) - as.matrix(center_df_new %>% dplyr::select(-label))
      res <- colSums(res)
      res <- sum(abs(res)) # 计算所有位置的偏差，如果小于10或1就停止
      print(
        paste("Step:", step, # 第几次迭代
              " Segment num:", nrow(center_df_new), # 迭代后后的中心点数，即分区数
              " NA spot: ", sum(is.na(OSE_meta_data$label)),
              " Diff:", res) # 迭代后与迭代前的中心点的差异
      )
      center_df <- center_df_new
      step <- step + 1
    }
    results_list$center_df <- center_df
    results_list$data <- OSE_meta_data
    write.table(center_df, file = paste0(output_dir,"/",samp_nm,"_center_df.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(OSE_meta_data, file = paste0(output_dir,"/",samp_nm,"_OSE_meta_data.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    #> 区域空间熵
    clog_step(paste0("Calculate spatial entropy..."))
    clog_normal("it may take a while...")

    # 计算矫正的系数
    combination <- function(n, k) {
      choose(n + k - 1, k)   # 总的niche可能性
    }
    # all_nicheType <- combination(4, 8)*4 # 如果周围都是点的话
    # 其实只要空白多就会熵高，哪怕不算这些不全的点，熵似乎会更高
    min_neighbor_num <- 3 # 3个点只是限制这是一个niche，实际上空白多似乎不会对熵有很大的影响
    all_nicheType <- sum(sapply(min_neighbor_num:8, function(x) combination(4, x)))*4 # 1920，考虑周围最少3个点的情况，确定了，就定为3个，
    expected_covered_species <- function(N, n) {
      N * (1 - (1 - 1/N)^n) # 指定抽样次数的期望覆盖种类数
    }
    celltype_mat <- Matrix::sparseMatrix(
      OSE_meta_data$row, # 不能为0，
      OSE_meta_data$col,
      x = as.numeric(factor(OSE_meta_data$celltype)), # value值？
    )
    entropy_results <- OSE_meta_data %>%
      dplyr::select(row, col, label, cellid, celltype) %>% # 这个因子顺序不重要，只要是不同的数字，后面确定niche
      mutate(celltype = as.numeric(factor(celltype))) %>% # 行数为点数
      mutate(
        neighbors = map2(row, col, ~ {
          celltype_code <-  celltype_mat[
            # 周围8个点的范围
            max(0, .x - coord_interval):min(.x + coord_interval, dim(celltype_mat)[1]), # 不要低于0，不要高于最大值，即范围内
            max(0, .y - coord_interval):min(.y + coord_interval, dim(celltype_mat)[2])
          ] %>% as.vector()

          celltype_code <- celltype_code[which(celltype_code > 0)] # 0就是dgc中的值为0，即细胞类型0，即没有细胞
          celltype_code <- celltype_code[-match(celltype_mat[.x, .y], celltype_code)] # 去掉一个自己
          table(celltype_code)
        })) %>%
      group_by(label) %>%
      mutate(N_Spots_raw = n()) %>%
      ungroup() %>%
      # 至少3个邻居的点，这在10x数据中是不能使用的，更改了这个，上面的all_nicheType就要更改
      # {if (data_platform == "Visium") filter(., N_Spots_raw >= min_neighbor_num) else .} %>% # !!!!!
      group_by(label,celltype,neighbors) %>% # 同时考虑了label和celltype和neighbors类型
      summarise(N_Spots_raw = first(N_Spots_raw),f_ij = n()) %>% # 每个区域中每种niche数量：中心+周围均不同
      group_by(label) %>% # orig.ident
      summarise(
        # 下面这些值都是每个区域的值
        N_Spots_raw = first(N_Spots_raw),
        N_niche = n(),
        N_Spots = sum(f_ij), # 是一个区域所有的点数，因为是sum,不一定是一个区域所有状态的数量
        entropy = -sum(f_ij/N_Spots * log2(f_ij/N_Spots)), # 每一行的熵，不同中心点+周围细胞类型
        entropy_adj_old = entropy/log2(N_Spots), # 标准化熵或归一化熵，除以已经出现的种类数？有和意义？都矫正了我们需要的差异？处于区间 [0, 1]
        expected_nicheTye = expected_covered_species(all_nicheType, N_Spots), # 期望的niche类型
        entropy_adj = entropy/log2(expected_nicheTye) # log2(expected_nicheTye)是最大的熵
        # 常见的标准化方法是将熵除以其最大可能值
        # 不同分布就用log2(种类可能数)，同一种分布，就用log2(点的总数比较)
      )
    small_v_nm <- entropy_results %>% filter(N_Spots < OSE_minSpotNum) %>% pull(label)
    results_list$entropy_results <- entropy_results
    write.table(entropy_results, file = paste0(output_dir,"/",samp_nm,"_entropy_results.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    entropy_results_fil <- entropy_results %>% filter(!(label %in% small_v_nm))
    results_list$entropy_results_fil <- entropy_results_fil
    write.table(entropy_results_fil,
                file = paste0(output_dir,"/",samp_nm,"_entropy_results_fil.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    #> 圆环数据
    clog_step(paste0("Calculate cell type percentage for pie chart..."))

    # celltype_percent
    OSE_meta_data <- OSE_meta_data %>%
      group_by(label) %>%
      mutate(label_num = n())
    label_df <- OSE_meta_data[c(20,22)] %>% unique()
    celltype_percent <- OSE_meta_data %>%
      dplyr::select(label,celltype) %>%
      group_by(label,celltype) %>%
      summarise(type_counts = n()) %>%
      group_by(label) %>%
      mutate(celltype_ratio = type_counts/sum(type_counts)) %>%
      ungroup() %>%
      dplyr::select(-type_counts) %>%
      tidyr::pivot_wider(names_from = celltype, values_from = celltype_ratio) %>%
      arrange(label) %>%
      mutate(entropy_adj = entropy_results$entropy_adj) %>%
      mutate(across(where(is.numeric), ~ifelse(is.na(.), 0, .)))
    celltype_percent <- celltype_percent %>%
      mutate(row = center_df_new$row[match(label, center_df_new$label)],
             col = center_df_new$col[match(label, center_df_new$label)])
    results_list$celltype_percent <- celltype_percent
    write.table(celltype_percent, file = paste0(output_dir,"/",samp_nm,"_celltype_percent.txt"),
                sep = "\t", col.names = T, row.names = F, quote =F)

    celltype_percent <- celltype_percent %>% # celltype_percent已经切换为celltype_percent_fil
      filter(!(label %in% small_v_nm))
    results_list$celltype_percent_fil <- celltype_percent
    write.table(celltype_percent, file = paste0(output_dir,"/",samp_nm,"_celltype_percent_fil.txt"),
                sep = "\t", col.names = T, row.names = F, quote =F)

    # celltype_percent_long
    # 使用的是celltype_percent_fil
    celltype_percent_long <- celltype_percent %>%
      pivot_longer(cols = c(-label,-row,-col,-entropy_adj), # 其他列都是细胞类型的比例
                   names_to = "type",
                   values_to = "percent") %>%
      mutate(type = gsub("of_", "", type)) %>%
      mutate(label_num = label_df$label_num[match(label, label_df$label)])
    results_list$celltype_percent_long <- celltype_percent_long
    write.table(celltype_percent_long, file = paste0(output_dir,"/",samp_nm,"_celltype_percent_long.txt"),
                sep = "\t", col.names = T, row.names = F, quote =F)
    plot_entropy_data <- distinct(celltype_percent_long, row,col,label, entropy_adj,label_num)

    #> 轮廓数据和图
    clog_step(paste0("Generate segment OSE_meta_data..."))
    if(data_platform == "Visium"){
      clog_normal("Visium OSE_meta_data detected, skipping silhouette plot generation...")
    }else{
      clog_normal("Generating segment OSE_meta_data for silhouette plot...")

      #> 轮廓图
      mat <- Matrix::sparseMatrix(
        i = OSE_meta_data$row,
        j = OSE_meta_data$col,
        x = as.numeric(factor(OSE_meta_data$label))
      )
      dim_mat <- dim(mat)

      # 横向邻居
      right <- mat[, -1, drop = FALSE]
      left <- mat[, -dim_mat[2], drop = FALSE]
      mask_right <- (left != 0 & right != 0 & left != right)
      cols <- which(mask_right, arr.ind = TRUE)
      right_segments <- tibble(
        x1 = cols[,1],
        y1 = cols[,2],
        x2 = cols[,1],
        y2 = cols[,2]+1
      )

      # 纵向邻居
      down <- mat[-1, , drop = FALSE]
      up <- mat[-dim_mat[1], , drop = FALSE]
      mask_down <- (up != 0 & down != 0 & up != down)
      rows <- which(mask_down, arr.ind = TRUE)
      down_segments <- tibble(
        x1 = rows[,1],
        y1 = rows[,2],
        x2 = rows[,1]+1,
        y2 = rows[,2]
      )
      segment_df <- bind_rows(right_segments, down_segments) %>% distinct()
      results_list$segment_df <- segment_df
      write.table(segment_df, file = paste0(output_dir,"/",samp_nm,"_segment.txt"),
                  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
  }else{
    clog_normal("only_plot is TRUE, skipping OSE calculation and using existing data for plotting...")
    center_df_new <- read.table(file = paste0(output_dir,"/",samp_nm,"_center_df.txt"),
                                header = T, sep = "\t",check.names = F)
    OSE_meta_data <- read.table(file = paste0(output_dir,"/",samp_nm,"_OSE_meta_data.txt"),
                       header = T, sep = "\t",check.names = F)
    OSE_meta_data <- OSE_meta_data %>%
      group_by(label) %>%
      mutate(label_num = n())
    label_df <- OSE_meta_data[c(20,22)] %>%unique()
    celltype_percent <- read.table(file = paste0(output_dir,"/",samp_nm,"_celltype_percent_fil.txt"), # 去掉了小于OSE_minSpotNum的点
                                   header = T, sep = "\t",check.names = F)
    segment_df <- read.table(file = paste0(output_dir,"/",samp_nm,"_segment.txt"),
                             header = T, sep = "\t",check.names = F)
    celltype_percent_long <- celltype_percent %>%
      pivot_longer(cols = c(-label,-row,-col,-entropy_adj), # 其他列都是细胞类型的比例
                   names_to = "type",
                   values_to = "percent") %>%
      mutate(type = gsub("of_", "", type)) %>%
      mutate(label_num = label_df$label_num[match(label, label_df$label)])
    plot_entropy_data <- distinct(celltype_percent_long, row,col,label, entropy_adj,label_num)
  }

  ### 后面的图都是去掉了小于OSE_minSpotNum的点了
  #> subplots
  clog_step(paste0("plot all subplots..."))
  if(!data_platform == "Visium"){
    #> 轮廓图
    clog_normal("Plotting silhouette plot...")
    p1 <- ggplot() +
      geom_segment( # 区域线
        data =  segment_df,
        # aes(x = y1-0.25, y = x1-0.25, xend = y2+0.25, yend = x2+0.25),
        aes(x = y1, y = x1, xend = y2, yend = x2),
        alpha = 1,
        color = "black",
        linewidth = 2
      )+
      coord_fixed()+  # 保持坐标轴比例一致
      theme_void() +
      scale_y_reverse()
    ggsave(plot = p1,filename = paste0(photo_dir,"/",samp_nm,"_segment.png"),
           width = 10, height = 10,bg = "transparent",dpi = 900)
    ggsave(plot = p1,filename = paste0(photo_dir,"/",samp_nm,"_segment.pdf"),
           width = 10, height = 10,bg = "transparent")
  }

  #> 比例饼图
  clog_normal("Plotting cell type ratio pie chart...")
  celltype_ratio <- table(OSE_meta_data$celltype) %>%
    as.data.frame() %>%
    mutate(Var1 = factor(Var1))
  celltype_ratio$percent <- celltype_ratio$Freq/sum(celltype_ratio$Freq)
  results_list$celltype_ratio <- celltype_ratio
  write.table(celltype_ratio, file = paste0(output_dir,"/",samp_nm,"_celltype_ratio.txt"),
              sep = "\t", quote = F, row.names = F)

  p2 <- ggplot(celltype_ratio,aes(x=1,y=percent,fill=Var1))+
    geom_bar(width = 1,stat = "identity")+
    coord_polar("y",start = 0,direction = -1)+
    theme_void()+
    theme(legend.position = "right")+
    scale_fill_manual(values = col)+
    labs(title = samp_nm,fill ="Group")+
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
  ggsave(plot = p2,filename = paste0(photo_dir,"/",samp_nm,"_celltype_ratio.pdf"),
         width = 4, height = 4)

  #> 分区和分组图
  # 分区可视化
  clog_normal("Plotting partition visualization...")
  p3 <- ggplot(OSE_meta_data, aes(x = col, y = row, color = label)) +
    geom_point(
      size = p3_size, # bin50,相比于绘图点要大一点才行，否则没法重叠，是被区域边界需要
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
  ggsave(plot = p3,filename = paste0(photo_dir,"/",samp_nm,"_region.pdf"),
         width = 8, height = 8)

  # 分类可视化
  clog_normal("Plotting group visualization...")
  plot_grp_data <- data.frame(
    row = OSE_meta_data$row, # row和col可能为0
    col = OSE_meta_data$col,
    group = OSE_meta_data$celltype
  )
  p4 <- ggplot() +
    geom_point(data =  plot_grp_data, aes(x = col, y = row, color = group),
               size = p4_size, # pub2/4/8
               shape = 16) +
    theme_void() +
    coord_fixed()+  # 保持坐标轴比例一致
    scale_color_manual(values = col) +
    scale_y_reverse()+
    labs(x = "Column", y = "Row", color = "Cell Type") +
    guides(color = "none")+
    coord_fixed()+  # 保持坐标轴比例一致
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.background = element_rect(fill = "black")
    )
  ggsave(plot = p4,filename = paste0(photo_dir,"/",samp_nm,"_grp.pdf"),
         width = 8, height = 8)


  #> 合并图
  clog_step(paste0("Plotting merged visualizations..."))
  # p1_merge: 皮髓分组图
  # p2_merge: 加饼图: 变动的值
  # p3_merge: 加熵图：变动的值
  # p4_merge: 加轮廓图

  #
  p1_merge <- ggplot() + # 分组图
    geom_point(data =  plot_grp_data, aes(x = col, y = row, color = group),
               size = p4_size, # pub2/4/8
               shape = 16) +
    theme_void() +
    coord_fixed()+  # 保持坐标轴比例一致
    scale_color_manual(values = alpha(col,0.8)) +
    scale_y_reverse()+
    labs(x = "Column", y = "Row", color = "Cell Type") +
    guides(color = "none")+
    coord_fixed()+  # 保持坐标轴比例一致
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.background = element_rect(fill = "black")
    )
  p2_merge <- p1_merge +
    new_scale_fill()+ # 饼图
    geom_arc_bar(
      data =  celltype_percent_long,
      aes(x0 = col, y0 = row,
          r0 = 0,
          r = scales::rescale(label_num/mean(label_num),to = c(0.8, 1.4))*p2_merge_r,
          amount = percent,,
          fill = type),
      # color = "grey90",
      alpha = 1,
      linewidth = 0.5, # linewidth, 旧版本没法识别？
      stat = "pie"
    )+
    coord_fixed()+  # 保持坐标轴比例一致
    scale_fill_manual(values = col) +
    labs(fill = "Cell Type")

  # entropy and no_seg
  plot_entropy_data <- plot_entropy_data %>%
    mutate(entropy_adj = ifelse(entropy_adj > 1, 1, entropy_adj))
  p3_merge <- p2_merge+
    new_scale_fill()+ # 熵的饼图
    geom_point(data =  plot_entropy_data,  # 确保每个 label 只绘制一个点
               aes(x = col, y = row, fill = entropy_adj),
               size = scales::rescale(plot_entropy_data$label_num/mean(plot_entropy_data$label_num),
                                      c(0.6,1.2))*p3_merge_size,
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
  ggsave(plot = p3_merge,filename = paste0(photo_dir,"/",samp_nm,"_no_seg.png"),
         width = 9, height = 8,bg = "transparent",dpi = 900)

  # merge all
  if(data_platform == "Visium"){
    clog_normal("Visium OSE_meta_data detected, skipping merge1 and merge2 generation...")
  }else{
    p4_merge <- p3_merge+
      geom_segment( # 区域线
        data = segment_df,
        # aes(x = y1-0.25, y = x1-0.25, xend = y2+0.25, yend = x2+0.25),
        aes(x = y1, y = x1, xend = y2, yend = x2),
        alpha = 1,
        color = "black",
        linewidth = 2
      )
    print(p4_merge)
    ggsave(plot = p4_merge,filename = paste0(photo_dir,"/",samp_nm,"_all_merge_white.pdf"),
           width = 9, height = 8)
    pdf(paste0(photo_dir,"/",samp_nm,"_all_merge_compare_white.pdf"),
        width = 16, height = 8)
    print(p1_merge + p4_merge)
    dev.off()

    # merge_black
    p3_merge_black <- p2_merge+
      new_scale_fill()+ # 熵
      geom_point(data =  plot_entropy_data,  # 确保每个 label 只绘制一个点
                 aes(x = col, y = row, fill = entropy_adj),
                 size =  scales::rescale(plot_entropy_data$label_num/mean(plot_entropy_data$label_num),
                                         to = c(0.6,1.2))*p3_merge_size,
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
        data =  segment_df,
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
    print(p3_merge_black)
    ggsave(plot = p3_merge_black,filename = paste0(photo_dir,"/",samp_nm,"_all_merge_black.pdf"),
           width = 9, height = 8)
  }

  # 10X: no_grp
  if(data_platform == "Visium"){
    clog_normal("Plotting merged visualization without group colors for Visium OSE_meta_data...")
    p_no_grp <- ggplot() +
      geom_arc_bar(
        data =  celltype_percent_long, # 饼图
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
      scale_fill_manual(values = col) +
      labs(fill = "Cell Type") +
      new_scale_fill()+ # 熵
      scale_y_reverse()+
      geom_point(data =  plot_entropy_data,  # 确保每个 label 只绘制一个点
                 aes(x = col, y = row, fill = entropy_adj),
                 size = scales::rescale(plot_entropy_data$label_num/mean(plot_entropy_data$label_num),
                                        to = c(0.6,1.5))*p3_merge_size,
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
    ggsave(plot = p_no_grp,filename = paste0(photo_dir,"/",samp_nm,"_all_merge.white.pdf"),
           width = 9, height = 8)
  }
  clog_normal("OSE analysis completed!")

  return(results_list)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalSampGeneTrend
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Calculate gene expression trends across time points
#'
#' Analyzes gene expression patterns across temporal samples using linear
#' regression or Mfuzz clustering to identify genes with different trend patterns.
#'
#' @param STID_obj A STID object
#' @param loop_id Character, multi-sample analysis identifier
#' @param samp_grp_index Logical, whether to group by sample groups instead of
#'        individual sample IDs (default: FALSE)
#' @param meta_key Character, metadata key for retrieving cell data
#' @param niche_key Character, niche key for filtering niche cells (optional)
#' @param group_by Character, column name for cell type grouping
#' @param assay_id Character, assay name (default: "Spatial")
#' @param layer_id Character, layer/slot name (default: "data")
#' @param method Character, trend analysis method - "fitting" (linear/non-linear)
#'        or "mfuzz" (fuzzy clustering) (default: "fitting")
#' @param fitting_paras List of parameters for fitting method:
#'        \itemize{
#'          \item lm_min_slope: Minimum slope threshold for significant linear trend
#'          \item lm_pval_cutoff: P-value cutoff for linear regression significance
#'        }
#' @param mfuzz_paras List of parameters for Mfuzz method:
#'        \itemize{
#'          \item cluster_num: Number of clusters for Mfuzz (default: 10)
#'          \item min_acore: Minimum membership score for gene inclusion (default: 0.3)
#'        }
#' @param gene_list Character vector, specific genes to analyze (default: NULL,
#'        uses all variable features)
#' @param remove_genes Character vector, genes to exclude from analysis
#' @param col Character vector, color palette for visualization
#' @param return_data Logical, whether to return results (default: FALSE)
#' @param grp_nm Character, group name for output (default: NULL, uses timestamp)
#' @param dir_nm Character, directory name for output (default: "M4_CalSampGeneTrend")
#'
#' @return If return_data = TRUE, returns a list of gene trend results per cell type
#'
#' @details
#' The fitting method uses linear regression to identify genes with significant
#' monotonic trends (increasing/decreasing) and quadratic regression to classify
#' concave/convex patterns. The Mfuzz method performs fuzzy clustering to group
#' genes with similar temporal expression patterns.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Linear regression based trend analysis
#' results <- CalSampGeneTrend(
#'   STID_obj = stid_object,
#'   loop_id = "time_series",
#'   group_by = "cell_type",
#'   method = "fitting",
#'   fitting_paras = list(lm_min_slope = 1.5, lm_pval_cutoff = 0.05)
#' )
#'
#' # Mfuzz clustering based trend analysis
#' results <- CalSampGeneTrend(
#'   STID_obj = stid_object,
#'   loop_id = "time_series",
#'   group_by = "cell_type",
#'   method = "mfuzz",
#'   mfuzz_paras = list(cluster_num = 8, min_acore = 0.4)
#' )
#' }
CalSampGeneTrend <- function(STID_obj = NULL, loop_id = NULL, # must support
                             samp_grp_index = FALSE,
                             meta_key = NULL,
                             niche_key = NULL, # only support one value
                             group_by = NULL, # 如果提供了，就是细胞级别的
                             assay_id = "Spatial",
                             layer_id = "data",
                             method = "fitting", # "mfuzz"/ "fitting"
                             fitting_paras = list(
                               lm_min_slope = 2,
                               lm_pval_cutoff = 0.2
                             ),
                             mfuzz_paras = list(
                               cluster_num = 10,
                               min_acore = 0.3
                             ),
                             gene_list = NULL,
                             remove_genes = NULL,
                             col = COLOR_LIST[["PALETTE_WHITE_BG"]],
                             return_data = FALSE,
                             grp_nm = NULL,dir_nm = "M4_CalSampGeneTrend"){

  # refer to CalSampDEGs
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(loop_id, col)
  match.arg(method, choices = c("fitting", "mfuzz"))
  valid_genes <- rownames(STID_obj)[!rownames(STID_obj) %in% remove_genes]
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  results_list <- list()
  clog_normal("Execute multi-sample analysis...")
  loop_multi <- .check_loop_multi(STID_obj = STID_obj, loop_id = loop_id)
  if(length(loop_multi) >1){
    clog_error("Multiple loop_id detected in STID object. Please specify one loop_id for analysis.")
  }

  #> main pipeline
  for(i in seq_along(loop_multi)){
    i_multi <- loop_multi[i]
    clog_loop(paste0("Processing multi_id: ", i_multi, " (", i, "/", length(loop_multi), ")"))
    samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
    samp_id2grp <- data.frame(
      samp_id = STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@samp_info$samp_id,
      samp_grp = STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@samp_info$samp_grp
    )

    #>
    if(!is.null(niche_key)){
      logic_niche <- .check_niche_exist(STID_obj, i_multi, niche_key, samp_mode  = "MS")
      if(!logic_niche){
        clog_warn(paste0("niche_key: ", niche_key, " not found in MultiSampNiche of multi_id: ", i_multi, ", skipping..."))
        next
      }
      clog_normal("Using niche_key for analysis...")
      Niche_cells <- GetMSNicheCells(STID_obj = STID_obj, niche_key = niche_key,loop_id = i_multi)[[1]]
    }else{
      if(is.null(meta_key)){
        clog_warn("Both niche_key and meta_key are NULL, will use meta_key: coord for MS niche plotting.")
        meta_key <- "coord"
      }
      clog_normal("Using meta_key for analysis...")
      Niche_cells <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]] %>%
        filter(!!sym(samp_colnm) %in% samp_id2grp$samp_id) # actually all, not niche
    }
    .check_column_exist(Niche_cells, group_by)
    if(samp_grp_index){
      Niche_cells <- Niche_cells %>%
        mutate(samp_grp = samp_id2grp$samp_grp[match(Niche_cells[[samp_colnm]],samp_id2grp$samp_id)],
               .after = all_of(samp_colnm))
      new_samp_colnm <- "samp_grp"
    }else{
      new_samp_colnm <- samp_colnm
    }
    if(!is.null(niche_key)){
      Niche_cells <- Niche_cells %>%
        mutate(Niche_label = ifelse
               (is_Niche,"Niche","Bystander"), .after = "is_Niche") %>%
        mutate(Niche_label = factor(Niche_label,levels = c("Bystander","Niche")))
    }else{
      Niche_cells <- Niche_cells %>%
        mutate(Niche_label = "All",.after = all_of(group_by)) %>%
        mutate(Niche_label = factor(Niche_label,levels = c("All")))
    }
    if(!is.null(group_by)){
      Niche_cells <- Niche_cells %>%
        mutate(Niche_celltype = paste0(Niche_label,"_",!!sym(group_by)))
    }else{
      Niche_cells <- Niche_cells %>%
        mutate(Niche_celltype = "All")
    }

    #> calculate gene expression
    seurat_obj <- as.Seurat(STID_obj)
    seurat_obj@meta.data <- Niche_cells # !!!
    clog_normal("Subsetting Seurat object for niche cells and valid genes...")
    seurat_obj <- subset(seurat_obj, subset = Niche_label != "Bystander",features = valid_genes)
    celltypes <- unique(seurat_obj@meta.data$Niche_celltype)
    seurat_obj@meta.data <- seurat_obj@meta.data %>%
      mutate(!!sym(new_samp_colnm) := factor(!!sym(new_samp_colnm),levels = samp_id2grp[[samp_grp_index + 1]])) %>%
      arrange(!!sym(new_samp_colnm)) %>%
      mutate(grp_celltype = paste0(!!sym(new_samp_colnm),"~",Niche_celltype))
    clog_normal("Aggregating gene expression by sample groups...")
    exp_list <- AggregateExpression(seurat_obj, group.by = "grp_celltype", features = gene_list,
                                    assays = assay_id, slot = layer_id)
    exp_df <- exp_list[[1]] %>% as.data.frame()
    colnames(exp_df) <- gsub("-","_",colnames(exp_df))

    for(i in seq_along(celltypes)){
      # i = 1
      celltype <- celltypes[i]
      clog_loop(paste0("Processing cell type: ", celltype, " (", i, "/", length(celltypes), ")"))
      output_dir_celltype <- paste0(output_dir,"/",celltype)
      photo_dir_celltype <- paste0(photo_dir,"/",celltype)
      dir.create(output_dir_celltype, showWarnings = FALSE, recursive = TRUE)
      dir.create(photo_dir_celltype, showWarnings = FALSE, recursive = TRUE)
      i_exp_df <- exp_df[,grepl(paste0("~",celltype,"$"),colnames(exp_df)),drop = FALSE]
      if(ncol(i_exp_df) < 3){
        clog_warn(paste0("The number of samples for cell type: ", celltype, " is less than 3, skipping mfuzz analysis..."))
        next
      }
      colnames(i_exp_df) <- gsub(paste0("~",celltype),"",colnames(i_exp_df))

      #> calculate gene trend
      if(method == "fitting"){
        clog_normal("Calculating gene trends using fitting method...")

        i_results <- .fitting_pipeline(exp_df = i_exp_df, fitting_paras = fitting_paras,
                                       output_dir = output_dir_celltype, photo_dir = photo_dir_celltype)
        results_list[[i_multi]]$data[[celltype]] <- i_results

      }else if(method == "mfuzz"){
        i_results <- .mfuzz_pipeline(exp = i_exp_df, grp = colnames(i_exp_df),
                                     cluster_num = mfuzz_paras$cluster_num,
                                     min_acore = mfuzz_paras$min_acore,
                                     output_dir = output_dir_celltype, photo_dir = photo_dir_celltype)
        results_list[[i_multi]]$data[[celltype]] <- i_results
      }
    }
  }

  #> save results
  clog_normal("Saving results to RDS file...")
  saveRDS(results_list, file = paste0(output_dir,"/CalSampGeneTrend_Results_List.rds"))

  # >>> Final
  .save_function_params("CalSampGeneTrend", envir = environment(),
                        file = paste0(output_dir,"Log_function_params_(CalSampGeneTrend).log"))
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(CalSampGeneTrend).log")) %>% invisible()
  if(return_data){
    return(results_list)
  }
}


#' Fitting pipeline for gene trend analysis
#'
#' Performs linear regression (slope) and quadratic regression (a coefficient)
#' on gene expression across time points to classify genes into four trend types:
#' up-concave, up-convex, down-concave, and down-convex.
#'
#' @param exp_df Data frame, expression matrix with genes as rows and time points
#'        as columns
#' @param fitting_paras List containing lm_min_slope and lm_pval_cutoff
#' @param output_dir Character, output directory for saving results
#' @param photo_dir Character, output directory for saving plots
#'
#' @return List containing:
#'   \itemize{
#'     \item all_results: Combined results from linear and quadratic regression
#'     \item up_concave_genes: Genes with increasing concave trends
#'     \item up_convex_genes: Genes with increasing convex trends
#'     \item down_concave_genes: Genes with decreasing concave trends
#'     \item down_convex_genes: Genes with decreasing convex trends
#'   }
#'
#' @importFrom viridis viridis
#' @importFrom ggrepel geom_text_repel
#'
#' @keywords internal
#' @noRd
.fitting_pipeline <- function(exp_df, fitting_paras, output_dir, photo_dir){

  results_list <- list()

  # lm
  clog_normal("Calculating gene trends using linear regression...")
  positions <- seq_len(ncol(exp_df))
  n_rows <- nrow(exp_df)
  pb <- txtProgressBar(min = 0, max = n_rows, style = 3)
  lm_results <- apply(exp_df, 1, function(x){
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    fit <- lm(x ~ positions)
    slope <- coef(fit)[2]
    pval <- summary(fit)$coefficients[2,4]
    return(c(slope = slope, pval = pval))
  })
  close(pb)
  lm_results <- t(lm_results) %>% as.data.frame() %>%
    dplyr::rename(slope = slope.positions)

  # nls
  clog_normal("Calculating gene trends using non-linear regression...")
  n <- nrow(exp_df)
  nls_results <- data.frame(
    gene = rownames(exp_df),
    a = NA_real_,
    b = NA_real_,
    c = NA_real_,
    nls_RSS = NA_real_,
    peak_pos = NA_real_,
    peak_height = NA_real_,
    peak_diff = NA_real_,
    stringsAsFactors = FALSE
  )
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in seq_len(n)) {
    gene <- rownames(exp_df)[i]
    y <- as.numeric(exp_df[gene, ])

    fit_result <- tryCatch({
      nls_fit <- nls(y ~ a*positions^2 + b*positions + c,
                     start = list(a=0, b=0, c=mean(y)),
                     control = nls.control(maxiter = 100, warnOnly = TRUE))

      coefs <- coef(nls_fit)
      peak_pos <- -coefs["b"] / (2 * coefs["a"])
      peak_height <- coefs["a"] * peak_pos^2 + coefs["b"] * peak_pos + coefs["c"]
      baseline <- mean(y[c(1, length(y))])
      peak_diff <- peak_height - baseline
      nls_RSS <- sum(residuals(nls_fit)^2)

      c(coefs, peak_pos = peak_pos, peak_height = peak_height, peak_diff = peak_diff, nls_RSS = nls_RSS)

    }, error = function(e) NULL)

    if (!is.null(fit_result)) {
      nls_results[i, 2:ncol(nls_results)] <- fit_result
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)
  combined_results <- lm_results %>%
    cbind(nls_results) %>%
    relocate(gene,.before = 1)

  #> process the results
  min_slope <- fitting_paras$lm_min_slope
  pval_cutoff <- fitting_paras$lm_pval_cutoff
  up_concave_genes <- combined_results %>%
    filter(slope > min_slope & pval < pval_cutoff) %>%
    filter(a < 0) %>%
    arrange(-slope) %>%
    mutate(top_gene = ifelse(row_number() <= 10, gene, NA))
  up_convex_genes <- combined_results %>%
    filter(slope > min_slope & pval < pval_cutoff) %>%
    filter(a > 0) %>%
    arrange(-slope) %>%
    mutate(top_gene = ifelse(row_number() <= 10, gene, NA))
  down_concave_genes <- combined_results %>%
    filter(slope < -min_slope & pval < pval_cutoff) %>%
    filter(a < 0) %>%
    arrange(slope) %>%
    mutate(top_gene = ifelse(row_number() <= 10, gene, NA))
  down_convex_genes <- combined_results %>%
    filter(slope < -min_slope & pval < pval_cutoff) %>%
    filter(a > 0) %>%
    arrange(slope) %>%
    mutate(top_gene = ifelse(row_number() <= 10, gene, NA))
  clog_normal(paste0("Identified ", nrow(up_concave_genes), " up-concave genes, ",
                     nrow(up_convex_genes), " up-convex genes, ",
                     nrow(down_concave_genes), " down-concave genes, and ",
                     nrow(down_convex_genes), " down-convex genes."))
  results_list <- list(
    all_results = combined_results,
    up_concave_genes = up_concave_genes,
    up_convex_genes = up_convex_genes,
    down_concave_genes = down_concave_genes,
    down_convex_genes = down_convex_genes
  )
  write.table(combined_results, file = paste0(output_dir,"/gene_trend_results.txt"),
              sep = "\t", quote = F, row.names = F)
  write.table(up_concave_genes, file = paste0(output_dir,"/up_concave_genes.txt"),
              sep = "\t", quote = F, row.names = F)
  write.table(up_convex_genes, file = paste0(output_dir,"/up_convex_genes.txt"),
              sep = "\t", quote = F, row.names = F)
  write.table(down_concave_genes, file = paste0(output_dir,"/down_concave_genes.txt"),
              sep = "\t", quote = F, row.names = F)
  write.table(down_convex_genes, file = paste0(output_dir,"/down_convex_genes.txt"),
              sep = "\t", quote = F, row.names = F)

  #> plot
  plot_data1 <- rbind(
    up_concave_genes %>% mutate(trend = "up_concave"),
    up_convex_genes %>% mutate(trend = "up_convex"),
    down_concave_genes %>% mutate(trend = "down_concave"),
    down_convex_genes %>% mutate(trend = "down_convex")
  ) %>%
    mutate(slope = ifelse(slope>quantile(slope, 0.995), quantile(slope, 0.995),
                          ifelse(slope<quantile(slope, 0.005), quantile(slope, 0.005), slope))) %>%
    mutate(a = ifelse(a>quantile(a, 0.995), quantile(a, 0.995),
                      ifelse(a<quantile(a, 0.005), quantile(a, 0.005), a)))
  p1 <- ggplot(plot_data1,aes(x = slope, y = a, color = -log10(pval)))+
    geom_point()+
    theme_bw()+
    labs(title = paste0("Gene trend classification"),
         x = "Slope of linear regression",
         y = "Quadratic coefficient (a) of non-linear regression")+
    scale_color_gradientn(colors = viridis(100, option = "H")[15:85])+
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))+
    geom_text_repel(aes(label=top_gene),  #显示感兴趣基因
                    size = 4,
                    fontface="italic",
                    color="grey50",
                    box.padding=unit(0.35, "lines"), #文本框周边填充
                    point.padding=unit(0.5, "lines"), #点周边填充
                    segment.colour = "grey50",  #连接点与标签的线段的颜色
                    min.segment.length = 0.1,
                    max.overlaps=100000) #最多点的数量
  print(p1)
  ggsave(plot = p1, filename = paste0(photo_dir,"/gene_trend_classification.pdf"),
         width = 6, height = 5)

  return(results_list)
}


#' Internal Mfuzz pipeline for gene trend analysis
#'
#' Performs fuzzy C-means clustering using the Mfuzz package to identify genes
#' with similar temporal expression patterns.
#'
#' @param exp Data frame, expression matrix with genes as rows and time points
#'        as columns
#' @param grp Character vector, time point labels
#' @param cluster_num Integer, number of clusters for Mfuzz (default: 20)
#' @param min_acore Numeric, minimum membership score for filtering genes
#'        (default: 0.3)
#' @param output_dir Character, output directory for saving results
#' @param photo_dir Character, output directory for saving plots
#'
#' @return Data frame with original expression, cluster assignments,
#'         and membership scores for all genes
#'
#' @importFrom Mfuzz filter.NA fill.NA filter.std standardise mestimate mfuzz
#' @importFrom Mfuzz mfuzz.plot2 acore mfuzzColorBar
#' @importFrom patchwork wrap_plots
#' @importFrom paletteer paletteer_c
#' @importFrom viridis viridis
#'
#' @keywords internal
#' @noRd
.mfuzz_pipeline <- function(exp = NULL,grp = NULL,
                           cluster_num = 20,min_acore = 0.3,
                           output_dir = NULL, photo_dir = NULL){

  #> process data
  clog_normal("Processing data for Mfuzz clustering...")
  colnames(exp) <- colnames(exp) %>% gsub("-","_",.) %>%
    make.names()
  grp <- grp %>% gsub("-","_",.)
  col_nm <- colnames(exp)
  exp_mean <- rowsum(t(exp), grp) / as.numeric(table(grp))
  exp_mean <- t(exp_mean)
  colnames(exp_mean) <- colnames(exp)

  #> Mfuzz ExpressionSet
  clog_normal("Performing Mfuzz clustering...")
  protein <- as.matrix(exp_mean)
  mfuzz_class <- new('ExpressionSet',exprs = protein)
  exprs <- Biobase::exprs
  # assign("exprs", Biobase::exprs, envir = .GlobalEnv)
  {
    mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
    mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
    mfuzz_class <- filter.std(mfuzz_class, min.std = 0.3,visu=F) # 推荐0.3
    mfuzz_class <- standardise(mfuzz_class)
  }
  set.seed(123)
  # cluster_num <- 20
  m <- mestimate(mfuzz_class)
  mfuzz_cluster <- mfuzz(mfuzz_class,
                         centers = cluster_num,
                         m = m)
  #> Mfuzz plot
  clog_normal("Plotting Mfuzz clustering results...")
  ifelse(cluster_num<=10,mfow_index <- c(5, 2),
         ifelse(cluster_num<=20,mfow_index <- c(5, 4),
                ifelse(cluster_num<=30,mfow_index <- c(5, 6),mfow_index <- c(5, 8))))
  width_index <- ifelse(cluster_num<=10,6,
                        ifelse(cluster_num<=20,12,
                               ifelse(cluster_num<=30,18,24)))
  pdf(file = paste0(photo_dir,"/mfuzz_merge_all_",cluster_num,".pdf"),
      width = width_index,height = 15)
  mfuzz.plot2(mfuzz_class,
              cl = mfuzz_cluster,
              mfrow = mfow_index,
              time.labels = colnames(protein),
              Xwidth=12,Xheight=10,
              min.mem = 0,
              x11=F)
  dev.off()
  pdf(file = paste0(photo_dir,"/mfuzz_merge_fil_",cluster_num,".pdf"),
      width = width_index,height = 15)
  mfuzz.plot2(mfuzz_class,
              cl = mfuzz_cluster,
              mfrow = mfow_index,
              time.labels = colnames(protein),
              Xwidth=12,Xheight=10,
              min.mem = 0.3,
              x11=F)
  dev.off()
  for (i in 1:cluster_num) {
    pdf(file = paste0(photo_dir,"/mfuzz_fil_",i,".pdf"),
        width = 10,height = 10)
    mat <- matrix(1:2,ncol=2,nrow=1,byrow=TRUE)
    layout(mat,widths=c(5,1))

    mfuzz.plot2(mfuzz_class,
                cl=mfuzz_cluster,
                mfrow=NA,
                colo="fancy",
                single=i,
                min.mem = 0.3,
                x11=F)
    mfuzzColorBar(col = "fancy",
                  main = "Membership",
                  cex.main = 1) # 添加颜色图例
    dev.off()
  }

  #> cluster of genes
  clog_normal("Extracting gene clusters from Mfuzz results...")
  {
    cluster_size <- mfuzz_cluster$size
    names(cluster_size) <- paste0("cluster_",1:cluster_num)
    cluster_size # 未过滤的
    write.table(cluster_size,
                file = paste0(output_dir,"/cluster_size_all.txt"),
                sep = "\t",col.names = F,row.names = T,quote = F)

    protein_mem <- mfuzz_cluster$membership %>%
      as.data.frame() %>%
      dplyr::rename_with(~paste0("cluster_",.x),everything())
  }
  {
    protein_cluster <- mfuzz_cluster$cluster
    protein_oriexp_cluster <- cbind(protein[names(protein_cluster), ],
                                    protein_cluster,
                                    protein_mem)
    write.table(protein_oriexp_cluster,
                file = paste0(output_dir,"/protein_oriexp_cluster.txt"),
                sep = '\t', col.names = NA, quote = FALSE)
  }
  {
    protein_cluster <- mfuzz_cluster$cluster
    protein_standard <- mfuzz_class@assayData$exprs
    protein_standard_cluster <- cbind(protein_standard[names(protein_cluster), ],
                                      protein_cluster,
                                      protein_mem) %>%
      as.data.frame()
    write.table(protein_standard_cluster,
                file = paste0(output_dir,"/protein_standard_cluster.txt"),
                sep = '\t', col.names = NA, quote = FALSE)
  }

  #> min_acore filtering
  clog_normal("Filtering gene clusters based on minimum membership score...")
  acore_fil <- acore(
    mfuzz_class,mfuzz_cluster,
    min.acore = min_acore # 根据membership值过滤基因
  )
  saveRDS(acore_fil,file = paste0(output_dir,"/acore_genes_fil.Rds"))
  names(acore_fil) <- paste0("cluster_",1:cluster_num)
  cluster_size_fil <- sapply(1:length(acore_fil),function(i){
    a <- nrow(acore_fil[[i]])
    names(a) <- names(acore_fil[i])
    return(a)
  })
  write.table(cluster_size_fil,
              file = paste0(output_dir,"/cluster_size_fil.txt"),
              sep = "\t",col.names = F,row.names = T,quote = F)

  #>
  clog_normal("Plotting gene expression trends for each cluster after filtering...")
  plot_list <- list()
  for (i_cluster in 1:cluster_num) {
    # i_cluster <- 1
    cluster_id <- paste0("cluster_",i_cluster)
    cat(paste0(cluster_id," start at ",Sys.time(),"\n"))

    dir.create(paste0(output_dir,"/",cluster_id))
    # dir.create(paste0(photo_dir,"/",cluster_id))
    acore_cluster <- acore_fil[[cluster_id]]$NAME
    fil_gene_cluster <- acore_fil[[cluster_id]] %>%
      arrange(desc(MEM.SHIP))
    write.table(fil_gene_cluster,
                file = paste0(output_dir,"/",cluster_id,"/fil_gene_",cluster_id,".txt"),
                sep = "\t",col.names = NA,row.names = T,quote = F)
    df_cluster <- protein_standard_cluster[acore_cluster,]
    write.table(df_cluster,
                file = paste0(output_dir,"/",cluster_id,"/protein_standard_cluster_",cluster_id,".txt"),
                sep = "\t",col.names = NA,row.names = T,quote = F)

    # 使用ggplot2绘制单个cluster的趋势图
    df_cluster_plot <- df_cluster[c(col_nm,cluster_id)] %>%
      dplyr::rename(cluster=all_of(cluster_id)) %>%
      mutate(gene=rownames(.)) %>%
      pivot_longer(
        cols = -c(cluster,gene),
        names_to = "group",
        values_to = "values") %>%
      mutate(group = factor(group,levels = col_nm)) %>%
      arrange(cluster) %>%
      mutate(gene = factor(gene,levels = unique(gene)))
    colors <- rev(paletteer_c("ggthemes::Red-Blue-White Diverging", 30)[-c(15,16)])
    # colors <- rev(paletteer_c("ggthemes::Red-Blue-White Diverging", 30))
    p1 <- ggplot(data = df_cluster_plot,aes(x=group,y=values,group=gene,color=cluster))+
      geom_line()+
      theme_classic()+
      scale_color_gradientn(colors = colorRampPalette(colors)(50),
                            limits=c(0,1),guide = "colorbar")+
      # scale_color_viridis_c(option = "H",
      #                       limits = c(0, 1) # bin200/bin50_rawgem，pub2/4/7/8
      #                       # limits = c(0, 0.8) # bin50
      # ) +
      scale_x_discrete(expand = c(0,0.1,0,0.1))+
      labs(x=NULL,y="Expression changes",title=paste0(cluster_id,": ",length(unique(df_cluster_plot$gene))),color="Membership",
           subtitle = paste0(rev(unique(df_cluster_plot$gene))[1:3],collapse = ","))+
      theme(#legend.position = "none",
        # panel.grid = element_blank(),
        plot.margin = margin(0.4,0.4,0.4,0.4,'cm'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(size = 10,face = 'bold',hjust = 1,angle = 45),
        # axis.text.x = element_text(size = 10,face = 'bold',hjust = 0.5),
        axis.text.y  = element_text(size = 10,face = 'bold',hjust = 0.5),
        legend.title = element_text(size = 11,face = 'bold',hjust = 0),
        legend.text = element_text(size = 10,face = 'bold',hjust = 0))
    plot_list[[i_cluster]] <- p1
    saveRDS(plot_list,file = paste0(output_dir,"/mfuzz_merge_fil_",cluster_num,"_ggplot2.Rds"))
  }
  pdf(file = paste0(photo_dir,"/mfuzz_merge_fil_",cluster_num,"_ggplot2.pdf"),
      width = width_index*1.2,height = 15)
  p1 <- wrap_plots(plotlist = plot_list,guides = "collect",
                   nrow = mfow_index[1],ncol = mfow_index[2])
  print(p1)
  dev.off()

  return(protein_oriexp_cluster)
}


