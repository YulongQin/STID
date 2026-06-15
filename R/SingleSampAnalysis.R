
#' @include RPyCall.R
NULL


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalSampComp: 第一种模式SS/MS，不需要多线程
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate Niche Cell Type Composition
#'
#' Analyzes and visualizes the cell type composition within and outside niches.
#' Generates stacked barplots and faceted barplots showing cell type proportions
#' for niche vs. bystander regions.
#'
#' @param STID_obj An STID object containing niche analysis results
#' @param samp_mode  Character, sample type - "SS" (single sample) or "MS" (multi-sample)
#'        (default: "SS")
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param samp_grp_index Logical, whether to group by sample groups in MS mode
#'        (default: FALSE)
#' @param meta_key Character, metadata key for MS mode when niche_key is NULL
#' @param niche_key Character, niche key to analyze (only one value supported)
#' @param group_by Character, column name for cell type grouping
#' @param col Color palette for visualization (default: COLOR_LIST$PALETTE_WHITE_BG)
#' @param return_data Logical, whether to return the plot data (default: FALSE)
#'
#' @return If return_data = TRUE, returns a list of plots per sample; otherwise NULL
#'
#' @importFrom ggh4x facet_grid2
#' @import dplyr
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Single-sample niche composition
#' CalSampComp(
#'   STID_obj = STID_obj,
#'   samp_mode  = "SS",
#'   niche_key = "niche_virulence",
#'   group_by = "cell_type"
#' )
#'
#' # Multi-sample niche composition with sample grouping
#' CalSampComp(
#'   STID_obj = STID_obj,
#'   samp_mode  = "MS",
#'   loop_id = "LoopAllMulti",
#'   niche_key = "niche_virulence",
#'   samp_grp_index = TRUE,
#'   group_by = "cell_type"
#' )
#' }
CalSampComp <- function(STID_obj = NULL,
                        samp_mode  = "SS",
                        loop_id = "LoopAllSamp", # LoopAllMulti
                        samp_grp_index = FALSE,
                        meta_key = NULL,
                        niche_key = NULL, # only support one value
                        group_by = NULL,
                        col = COLOR_LIST[["PALETTE_WHITE_BG"]],
                        return_data = FALSE
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(loop_id,col)
  match.arg(samp_mode, .STID_globals$samp_modes)
  if(is.null(group_by)){
    clog_warn("group_by is NULL, will use the default celltype_colnm from STID object.")
    group_by <- GetInfo(STID_obj, info_key = "data_info",sub_key = "celltype_colnm")[[1]]
  }
  # >>> End check

  # >>> Start main pipeline
  results_list <- list()
  if(samp_mode  == "SS"){
    clog_normal("Execute single-sample analysis...")
    .check_null_args(niche_key);.check_one_arg(niche_key)
    loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)

    #> main pipeline
    for(i in seq_along(loop_single)){
      i_single <- loop_single[i]
      clog_loop(paste0("Processing samp_id: ", i_single, " (", i, "/", length(loop_single), ")"))
      logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
      if(!logic_niche){
        clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
        next
      }
      Niche_cells <- GetSSNicheCells(STID_obj = STID_obj, niche_key = niche_key,loop_id = i_single)[[1]]
      .check_column_exist(Niche_cells, group_by)
      Niche_cells <- Niche_cells %>%
        mutate(Niche_label = ifelse(is_Niche,"Niche","Bystander"),.after = "is_Niche") %>%
        mutate(Niche_label = factor(Niche_label,levels = c("Bystander","Niche")))

      clog_normal("Plotting the stacked barplot...")
      p1 <- ggplot(Niche_cells,aes( x=Niche_label	,fill=!!sym(group_by)))+
        geom_bar(position="fill")+
        scale_color_manual(values=col)+
        scale_fill_manual(values=col)+
        labs(y="Percentage of groups",x=NULL,title=i_single,color=NULL,fill=NULL)+
        theme_common()
      print(p1)

      # >
      clog_normal("Plotting the facet barplot...")
      Niche_cells_summary <- Niche_cells %>%
        group_by(Niche_label) %>%
        mutate(cellnum = n()) %>%
        ungroup() %>%
        group_by(!!sym(group_by),Niche_label) %>%
        summarise(celltype_ratio = n()/dplyr::first(cellnum)) %>%
        as.data.frame()
      p2 <- ggplot(Niche_cells_summary,aes( x=Niche_label,y=celltype_ratio, fill=!!sym(group_by)))+
        facet_wrap(vars(!!sym(group_by)), ncol = 3,scales = "free_y")+
        geom_col(position = "identity")+
        scale_y_continuous(
          # expand = c(expand = c(0,0,0,0.1)),
          labels = scales::percent)+
        scale_color_manual(values=col)+
        scale_fill_manual(values=col)+
        labs(y="Percentage of groups",x=NULL,title=i_single,color=NULL,fill=NULL)+
        theme_common()
      print(p2)

      results_list[[i_single]] <- list(p1 = p1,p2 = p2)
    }
  }else if(samp_mode  == "MS"){
    clog_normal("Execute multi-sample analysis...")
    loop_multi <- .check_loop_multi(STID_obj = STID_obj, loop_id = loop_id)

    #> main pipeline
    for(i in seq_along(loop_multi)){
      i_multi <- loop_multi[i]
      clog_loop(paste0("Processing samp_id: ", i_multi, " (", i, "/", length(loop_multi), ")"))
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
          clog_warn("Both niche_key and meta_key are NULL, will use meta_key: coord for analysis...")
          meta_key <- "coord"
        }
        clog_normal("Using meta_key for analysis...")
        Niche_cells <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]] %>%
          filter(!!sym(samp_colnm) %in% samp_id2grp$samp_id) # actually all, not niche
      }
      .check_column_exist(Niche_cells, group_by)
      if(samp_grp_index){
        Niche_cells <- Niche_cells %>%
          mutate(samp_grp = samp_id2grp$samp_grp[match(Niche_cells[[samp_colnm]],samp_id2grp$samp_id)],.after = all_of(samp_colnm))
        new_samp_colnm <- "samp_grp"
      }else{
        new_samp_colnm <- samp_colnm
      }
      if(!is.null(niche_key)){
        Niche_cells <- Niche_cells %>%
          mutate(Niche_label = ifelse(is_Niche,"Niche","Bystander"),.after = "is_Niche") %>%
          mutate(Niche_label = factor(Niche_label,levels = c("Bystander","Niche")))
      }else{
        Niche_cells <- Niche_cells %>%
          mutate(Niche_label = "All",.after = all_of(group_by)) %>%
          mutate(Niche_label = factor(Niche_label,levels = c("All")))
      }

      #>
      clog_normal("Plotting the stacked barplot...")
      p1 <- ggplot(Niche_cells,aes( x=!!sym(new_samp_colnm)	,fill=!!sym(group_by)))+
        geom_bar(position="fill")+
        # scale_color_manual(values=col)+
        scale_fill_manual(values=col)+
        labs(y="Percentage of groups",x=NULL,title=i_multi,color=NULL,fill=NULL)+
        theme_common()
      if(!is.null(niche_key)){
        p1 <- p1 + facet_grid2(cols = vars(Niche_label),
                               scales = "free", independent = "y")
      }
      print(p1)

      # >
      clog_normal("Plotting the facet barplot...")
      if(!is.null(niche_key)){
        Niche_cells_summary <- Niche_cells %>%
          group_by(!!sym(new_samp_colnm),Niche_label) %>%
          mutate(cellnum = n()) %>%
          ungroup() %>%
          group_by(!!sym(new_samp_colnm),!!sym(group_by),Niche_label) %>%
          summarise(celltype_ratio = n()/dplyr::first(cellnum)) %>%
          as.data.frame()
      }else{
        Niche_cells_summary <- Niche_cells %>%
          group_by(!!sym(new_samp_colnm)) %>%
          mutate(cellnum = n()) %>%
          ungroup() %>%
          group_by(!!sym(new_samp_colnm),!!sym(group_by)) %>%
          summarise(celltype_ratio = n()/dplyr::first(cellnum)) %>%
          as.data.frame()
      }
      p2 <- ggplot(Niche_cells_summary,aes( x=!!sym(new_samp_colnm),y=celltype_ratio, fill=!!sym(group_by)))+
        geom_col(position = "identity")+
        scale_y_continuous(
          # expand = c(expand = c(0,0,0,0.1)),
          labels = scales::percent)+
        # scale_color_manual(values=col)+
        scale_fill_manual(values=col)+
        labs(y="Percentage of groups",x=NULL,title=i_multi,color=NULL,fill=NULL)+
        theme_common()
      if(!is.null(niche_key)){
        p2 <- p2 + facet_grid2(rows = vars(!!sym(group_by)),cols = vars(Niche_label),
                               scales = "free", independent = "y")
      }else{
        p2 <- p2 + facet_grid2(rows = vars(!!sym(group_by)),
                               scales = "free", independent = "y")
      }
      print(p2)

      results_list[[i_multi]] <- list(p1 = p1,p2 = p2)
    }
  }

  # >>> Final
  clog_end()
  if(return_data){
    return(results_list)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalSampCAI
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate Niche Aggregation Index
#'
#' Computes spatial aggregation metrics for cell types within niches using
#' k-nearest neighbors and graph-based clustering. Calculates aggregation index
#' and aggregation region ratios.
#'
#' @param STID_obj An STID object containing niche analysis results
#' @param samp_mode  Character, sample type - "SS" (single sample) or "MS" (multi-sample)
#'        (default: "SS")
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param samp_grp_index Logical, whether to group by sample groups in MS mode
#'        (default: FALSE)
#' @param meta_key Character, metadata key for MS mode when niche_key is NULL
#' @param niche_key Character, niche key to analyze
#' @param group_by Character, column name for cell type grouping
#' @param dist_thres Numeric, distance threshold for graph construction
#'        (default: NULL, uses median of 1st neighbor distances)
#' @param k_neighbors Integer, number of nearest neighbors to consider (default: 8)
#' @param min_agg_size Integer, minimum aggregation cluster size (default: 2)
#' @param col Color palette for visualization (default: \code{COLOR_LIST$PALETTE_WHITE_BG[-1]})
#' @param return_data Logical, whether to return the results list (default: FALSE)
#'
#' @return If return_data = TRUE, returns a list of results per sample; otherwise NULL
#'
#' @import dplyr
#' @import FNN
#' @import dbscan
#' @importFrom igraph graph_from_adjacency_matrix V E
#' @import ggsignif
#' @importFrom ggh4x facet_wrap2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate aggregation index for niche cells
#' CalSampCAI(
#'   STID_obj = STID_obj,
#'   samp_mode  = "SS",
#'   niche_key = "niche_virulence",
#'   group_by = "cell_type",
#'   k_neighbors = 8,
#'   min_agg_size = 3
#' )
#' }
CalSampCAI <- function(STID_obj = NULL,
                       samp_mode  = "SS",
                       loop_id = "LoopAllSamp",
                       samp_grp_index = FALSE,
                       meta_key = NULL,
                       niche_key = NULL,
                       group_by = NULL,
                       dist_thres = NULL,
                       k_neighbors = 8,
                       min_agg_size = 2,
                       col = COLOR_LIST[["PALETTE_WHITE_BG"]][-1],
                       return_data = FALSE
) {
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(loop_id,col)
  match.arg(samp_mode, .STID_globals$samp_modes)
  if(is.null(group_by)){
    clog_warn("group_by is NULL, will use the default celltype_colnm from STID object.")
    group_by <- GetInfo(STID_obj, info_key = "data_info",sub_key = "celltype_colnm")[[1]]
  }
  # >>> End check

  # >>> Start main pipeline
  clog_step("Calculating the cell aggregation index...")
  if(samp_mode  == "SS"){
    clog_normal("Execute single-sample analysis...")
    .check_null_args(niche_key)
    .check_one_arg(niche_key)
    loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)
  }else if(samp_mode  == "MS"){
    clog_normal("Execute multi-sample analysis...")
    loop_single <- lapply(STID_obj_MS@STID_analysis@MultiSampNiche,function(x){ x@samp_info$samp_id}) %>%
      unlist(use.names = F) %>% unique()
  }

  # >>> plot_data
  all_metrics <- data.frame()
  for (i in seq_along(loop_single)) {
    i_single <- loop_single[i]
    clog_loop(paste0("Processing samp_id: ", i_single, " (", i, "/", length(loop_single), ")"))

    #>
    if(!is.null(niche_key)){
      logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
      if(!logic_niche){
        clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
        next
      }
      clog_normal("Using niche_key for SS/analysis...")
      Niche_cells <- GetSSNicheCells(STID_obj = STID_obj, niche_key = niche_key, loop_id = i_single)[[1]]
    }else{
      if(is.null(meta_key)){
        clog_warn("Both niche_key and meta_key are NULL, will use meta_key: coord for analysis...")
        meta_key <- "coord"
      }
      clog_normal("Using meta_key for analysis...")
      Niche_cells <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]] %>%
        filter(!!sym(samp_colnm) == i_single) # actually all, not niche
    }
    .check_column_exist(Niche_cells, group_by)
    if(!is.null(niche_key)){
      Niche_cells <- Niche_cells %>%
        mutate(Niche_label = ifelse(is_Niche,"Niche","Bystander"),.after = "is_Niche") %>%
        mutate(Niche_label = factor(Niche_label,levels = c("Bystander","Niche")))
    }else{
      Niche_cells <- Niche_cells %>%
        mutate(Niche_label = "All",.after = group_by) %>%
        mutate(Niche_label = factor(Niche_label,levels = c("All")))
    }
    i_Niche_cells <- Niche_cells %>% filter(!is.na(.data[[group_by]]))
    if (nrow(i_Niche_cells) == 0) {
      clog_warn(paste("No valid cells in", group_by, "for sample", i_single)) # ???
      next
    }

    #> Loop group
    unique_groups <- unique(i_Niche_cells[[group_by]]) %>% sort() %>% as.character()
    for ( j in seq_along(unique_groups)) {
      j_group <- unique_groups[j]
      clog_normal(paste0("Processing group: ", j_group, " (", j, "/", length(unique_groups), ")"))
      j_Niche_cells <- i_Niche_cells %>% filter(.data[[group_by]] == j_group)

      #> Loop region
      unique_regions <- unique(j_Niche_cells$Niche_label) %>% sort() %>% as.character()
      for (k in seq_along(unique_regions)) {
        k_region <- unique_regions[k]
        sub_cells <- j_Niche_cells %>% filter(Niche_label == k_region)
        clog_normal(paste0("Spots in ", k_region, ": ", nrow(sub_cells)))
        if (nrow(sub_cells) < k_neighbors + 1) {
          clog_warn(paste("Not enough cells for group =", j_group,"in", k_region,", skipping..."))
          next
        }
        coords <- sub_cells[, c("x", "y"), drop = FALSE]
        coords_mat <- as.matrix(coords)

        #> Cell Density
        knn_result <- get.knn(coords_mat, k = min(k_neighbors, nrow(coords) - 1))
        avg_distances <- rowMeans(knn_result$nn.dist)
        avg_distances_exp <- sapply(avg_distances, function(x) exp(-x + 1)) # Gaussian kernel
        density_metric <- data.frame(
          samp_id = i_single,
          group_id = j_group,
          Niche_label = k_region,
          spot_dist = avg_distances,
          spot_agg_index = avg_distances_exp,
          stringsAsFactors = FALSE,
          row.names = rownames(coords)
        )

        #> Aggregation via Graph
        if(is.null(dist_thres)){
          dist_thres <- quantile(knn_result$nn.dist[, 1], 0.5)
        }
        n_cells <- nrow(coords)
        adj_matrix <- matrix(FALSE, n_cells, n_cells)
        if(0){
          for (idx in 1:n_cells) {
            neighbors_idx <- knn_result$nn.index[idx, ]
            neighbors_idx <- neighbors_idx[knn_result$nn.dist[idx, ] <= dist_thres]
            adj_matrix[idx, neighbors_idx] <- TRUE
          }
        }else{
          valid_mask <- knn_result$nn.dist <= dist_thres
          rows <- row(knn_result$nn.index)[valid_mask]
          cols <- knn_result$nn.index[valid_mask]
          adj_matrix[cbind(rows, cols)] <- TRUE
        }
        adj_matrix <- adj_matrix | t(adj_matrix) #
        diag(adj_matrix) <- FALSE
        adj_matrix <- as(adj_matrix, "dgCMatrix")
        g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected") # !!! time-consuming step
        comps <- igraph::components(g)$membership
        agg_metric <- data.frame(
          components = comps,
          stringsAsFactors = FALSE,
          row.names = rownames(coords)
        )
        niche_metrics <- bind_cols(density_metric, agg_metric)
        all_metrics <- bind_rows(all_metrics, niche_metrics)
      }
    }
  }

  # >>> plot_data
  clog_step("Plotting the aggregation index and aggregation region ratio...")
  results_list <- list()
  if(samp_mode  == "SS"){
    clog_normal("Execute single-sample analysis...")
    for(i in seq_along(loop_single)){
      i_single <- loop_single[i]
      logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
      clog_loop(paste0("Plotting for samp_id: ", i_single, " (", i, "/", length(loop_single), ")"))
      if(!logic_niche){ next }
      plot_data1 <- all_metrics %>%
        filter(samp_id == i_single)
      plot_data2 <- plot_data1 %>%
        group_by(Niche_label, group_id) %>%
        mutate(all_agg_size = n()) %>%
        group_by(Niche_label, group_id, components) %>%
        mutate(i_agg_size = n()) %>%
        filter(i_agg_size >= min_agg_size) %>%
        group_by(Niche_label, group_id) %>%
        summarise(agg_ratio = n()/dplyr::first(all_agg_size))

      #> plot
      clog_normal("Plotting the aggregation index boxplot...")
      p1 <- ggplot(plot_data1,aes(x=Niche_label, y=spot_agg_index))+
        facet_wrap( ~ group_id, ncol = 3,scales = "free_y")+
        stat_boxplot(
          mapping = aes(color=group_id),
          geom = 'errorbar',
          width=0.7,
          linewidth=0.5,
          position = position_dodge(width=0.8) # 控制组间宽度
        )+ # 1.5 IQR横线
        geom_boxplot(
          mapping = aes(color=group_id),
          width=0.7,
          linewidth=0.5,
          # outliers = F, # 是否显示异常值，只影响展示，不影响计算
          outlier.fill = "white",
          outlier.color = "white",
          outlier.size = 0.5,
          outlier.stroke = 0.3,
          outlier.shape=21,
          position = position_dodge(width=0.8),
        )+
        scale_color_manual(values=col)+
        scale_fill_manual(values=col)+
        scale_y_continuous(expand = c(0,0.07,0,0.15)) +
        labs(y= "Average Aggregation Index",x=NULL,title=paste0("Average Aggregation Index (", i_single,")"),
             color=NULL,fill=NULL)+
        theme_common() +
        geom_signif(comparisons = list(c("Bystander", "Niche")),
                    test = wilcox.test, # t.test/wilcox.test
                    test.args = c("two.sided"), # two.sided, greater, less
                    map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
                    size = 0.2, # 线粗
                    textsize = 3, # 文本/***大小
                    # y_position = 6, # 起始位置，default：自动判断最高图形高度
                    margin_top = 0.2, # 实际标记其实位置：y_position + margin_top*y_position
                    step_increase = 0.1, # 加一次比较，加总高度的比例
                    tip_length = 0.05,
                    vjust = -0.2,  # 文字相对于横线上下移动文本，负值向上
                    color="black")
      print(p1)

      clog_normal("Plotting the aggregation ratio barplot...")
      p2 <- ggplot(plot_data2,aes( x=Niche_label,y=agg_ratio, fill=group_id))+
        facet_wrap( ~ group_id, ncol = 3,scales = "free_y")+
        geom_col(position = "identity")+
        scale_color_manual(values=col)+
        scale_fill_manual(values=col)+
        scale_y_continuous(labels = scales::percent)+
        labs(y= "Aggregation Region Ratio",x=NULL,title=paste0("Aggregation Region Ratio (", i_single,")"),
             color=NULL,fill=NULL)+
        theme_common()
      print(p2)

      results_list[[i_single]] <- list(data = list(all_metrics = all_metrics),
                                       plot = list(plot_spot_dit = p1, Plot_agg_ratio = p2))
    }
  }else if(samp_mode  == "MS"){
    clog_normal("Execute multi-sample analysis...")
    loop_multi <- .check_loop_multi(STID_obj = STID_obj, loop_id = loop_id)
    for(i in seq_along(loop_multi)){
      i_multi <- loop_multi[i]
      clog_loop(paste0("Processing samp_id: ", i_multi, " (", i, "/", length(loop_multi), ")"))
      samp_id2grp <- data.frame(
        samp_id = STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@samp_info$samp_id,
        samp_grp = STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@samp_info$samp_grp
      )
      i_all_metrics <- all_metrics %>% filter(samp_id %in% samp_id2grp$samp_id)
      if(samp_grp_index){
        i_all_metrics <- i_all_metrics %>%
          mutate(samp_grp = samp_id2grp$samp_grp[match(i_all_metrics[["samp_id"]],samp_id2grp$samp_id)],.after = "samp_id")
        new_samp_colnm <- "samp_grp"
      }else{
        new_samp_colnm <- "samp_id"
      }

      #>
      plot_data1 <- i_all_metrics
      if(!is.null(niche_key)){
        logic_niche <- .check_niche_exist(STID_obj, i_multi, niche_key, samp_mode  = "MS")
        if(!logic_niche){
          clog_warn(paste0("niche_key: ", niche_key, " not found in MultiSampNiche of multi_id: ", i_multi, ", skipping..."))
          next
        }
        plot_data2 <- plot_data1 %>%
          group_by(!!sym(new_samp_colnm),Niche_label, group_id) %>%
          mutate(all_agg_size = n()) %>%
          group_by(!!sym(new_samp_colnm),Niche_label, group_id, components) %>%
          mutate(i_agg_size = n()) %>%
          filter(i_agg_size >= min_agg_size) %>%
          group_by(!!sym(new_samp_colnm),Niche_label, group_id) %>%
          summarise(agg_ratio = n()/dplyr::first(all_agg_size))
      }else{
        plot_data2 <- plot_data1 %>%
          group_by(!!sym(new_samp_colnm), group_id) %>%
          mutate(all_agg_size = n()) %>%
          group_by(!!sym(new_samp_colnm), group_id, components) %>%
          mutate(i_agg_size = n()) %>%
          filter(i_agg_size >= min_agg_size) %>%
          group_by(!!sym(new_samp_colnm), group_id) %>%
          summarise(agg_ratio = n()/dplyr::first(all_agg_size))
      }

      #> plot
      clog_normal("Plotting the aggregation index boxplot...")
      p1 <- ggplot(plot_data1,aes(x=!!sym(new_samp_colnm), y=spot_agg_index))+
        stat_boxplot(
          mapping = aes(color=group_id),
          geom = 'errorbar',
          width=0.7,
          linewidth=0.5,
          position = position_dodge(width=0.8) # 控制组间宽度
        )+ # 1.5 IQR横线
        geom_boxplot(
          mapping = aes(color=group_id),
          width=0.7,
          linewidth=0.5,
          # outliers = F, # 是否显示异常值，只影响展示，不影响计算
          outlier.fill = "white",
          outlier.color = "white",
          outlier.size = 0.5,
          outlier.stroke = 0.3,
          outlier.shape=21,
          position = position_dodge(width=0.8),
        )+
        scale_color_manual(values=col)+
        # scale_fill_manual(values=col)+
        scale_y_continuous(expand = c(0,0.07,0,0.15)) +
        labs(y= "Average Aggregation Index",x=NULL,title=paste0("Average Aggregation Index (", i_multi,")"),
             color=NULL,fill=NULL)+
        theme_common()
      if(!is.null(niche_key)){
        p1 <- p1 + facet_grid2(rows = vars(group_id),cols = vars(Niche_label),
                               scales = "free", independent = "y")
      }else{
        p1 <- p1 + facet_grid2(rows = vars(group_id),
                               scales = "free", independent = "y")
      }
      print(p1)

      clog_normal("Plotting the aggregation ratio barplot...")
      p2 <- ggplot(plot_data2,aes( x=!!sym(new_samp_colnm),y=agg_ratio, fill=group_id))+
        geom_col(position = "identity")+
        # scale_color_manual(values=col)+
        scale_fill_manual(values=col)+
        scale_y_continuous(labels = scales::percent)+
        labs(y= "Aggregation Region Ratio",x=NULL,title=paste0("Aggregation Region Ratio (", i_multi,")"),
             color=NULL,fill=NULL)+
        theme_common()
      if(!is.null(niche_key)){
        p2 <- p2 + facet_grid2(rows = vars(group_id),cols = vars(Niche_label),
                               scales = "free", independent = "y")
      }else{
        p2 <- p2 + facet_grid2(rows = vars(group_id),
                               scales = "free", independent = "y")
      }
      print(p2)

      results_list[[i_multi]] <- list(data = list(all_metrics = all_metrics),
                                      plot = list(plot_spot_dit = p1, Plot_agg_ratio = p2))
    }
  }
  clog_normal("Saving results to RDS file...")
  saveRDS(results_list, file = file.path(output_dir, paste0("NicheAggIndex_results_", loop_id, ".rds")))

  # >>> Final
  clog_end()
  if (return_data) {
    return(results_list)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalSampCoLoc：细胞间/基因间/基因分数间, noMS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate Niche Co-localization
#'
#' Performs spatial co-localization analysis using the mistyR or Squidpy method to
#' assess the spatial relationships between cell types, gene expression,
#' or gene scores within niches. Generates co-localization metrics and visualizations for each sample.
#'
#' @param STID_obj An STID object containing niche analysis results
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param niche_key Character, niche key to analyze
#' @param meta_key Character, metadata key for analysis when niche_key is NULL
#' @param group_by Character, column name for cell type grouping
#' @param group_use Character vector, specific groups to include
#' @param features Character vector, gene names to analyze
#' @param feature_colnm Character vector, metadata column names to analyze
#' @param method Character, method for co-localization analysis - "mistyR" or "squidpy" (default: "mistyR")
#' @param mistyR_params Numeric vector, parameters for MSTIDyR analysis (default: c(juxtaview_radius = 15, paraview_radius = 10))
#'       - juxtaview_radius: radius for juxtaview in spatial units
#'       - paraview_radius: radius for paraview in spatial units
#' @param return_data Logical, whether to return the results list (default: TRUE)
#' @param grp_nm Character, group name for output organization (default: NULL)
#' @param dir_nm Character, directory name for output (default: "M3_CalSampCoLoc")
#'
#' @return If return_data = TRUE, returns a list of MSTIDyR results per sample
#'
#' @import distances
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze co-localization of cell types within niches
#' results <- CalSampCoLoc(
#'   STID_obj = STID_obj,
#'   niche_key = "niche_virulence",
#'   group_by = "cell_type",
#'   mistyR_params = c(15,10)
#' )
#' }
CalSampCoLoc <- function(STID_obj = NULL,
                         loop_id = "LoopAllSamp",
                         meta_key = NULL,
                         niche_key = NULL,
                         group_by = NULL,
                         group_use = NULL,
                         features = NULL,
                         feature_colnm = NULL,
                         method = "mistyR",
                         mistyR_params = c(juxtaview_radius = 15,
                                           paraview_radius = 10),
                         return_data = TRUE,
                         grp_nm = NULL,dir_nm = "M3_CalSampCoLoc"
){
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
  # .check_null_args(niche_key)
  match.arg(method, c("mistyR", "squidpy"))
  if(.all_null(group_by,features,feature_colnm)){
    clog_warn("group_by, features and feature_colnm are all NULL, will use
              the default celltype_colnm from STID object as group_by.")
    group_by <- GetInfo(STID_obj, info_key = "data_info",sub_key = "celltype_colnm")[[1]]
  }
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  coord_interval <- GetInfo(STID_obj, info_key = "data_info",sub_key = "coord_interval")[[1]]
  if(method == "mistyR"){
    clog_step("Calculating co-localization using MSTIDyR...")
    if( requireNamespace("mistyR", quietly = TRUE) == FALSE){
      clog_error("Package 'mistyR' is required but not installed. Please install it first.")
    }
    results_list <- .mistyR_pipeline(STID_obj = STID_obj,
                                     loop_single = loop_single,
                                     niche_key = niche_key,
                                     meta_key = meta_key,
                                     group_by = group_by,
                                     group_use = group_use,
                                     features = features,
                                     feature_colnm = feature_colnm,
                                     mistyR_params = mistyR_params,
                                     coord_interval = coord_interval,
                                     output_dir = output_dir,
                                     photo_dir = photo_dir)
    clog_normal("MSTIDyR analysis completed.")
  }else if(method == "squidpy"){
    clog_step("Calculating co-localization using Squidpy...")
    if(!.all_null(features,feature_colnm)){
      clog_warn("features and feature_colnm are not NULL, but Squidpy method does not support them, will ignore these parameters.")
    }
    results_list <- .RPyCall_squidpy(STID_obj = STID_obj,
                                     loop_single = loop_single,
                                     niche_key = niche_key,
                                     meta_key = meta_key,
                                     group_by = group_by,
                                     group_use = group_use,
                                     coord_interval = coord_interval,
                                     output_dir = output_dir,
                                     photo_dir = photo_dir)
    clog_normal("Squidpy analysis completed.")
  }
  clog_normal("Saving results to RDS file...")
  saveRDS(results_list, file = file.path(output_dir, paste0("NicheCoLoc_results_", loop_id, ".rds")))

  # >>> Final
  .save_function_params("CalSampCoLoc", envir = environment(), file = paste0(output_dir,"Log_function_params_(CalSampCoLoc).log"))
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(CalSampCoLoc).log")) %>% invisible()
  if(return_data){
    return(results_list)
  }
}


#' Internal pipeline for MSTIDyR spatial interaction analysis
#'
#' @param STID_obj STID object
#' @param loop_single Sample identifiers
#' @param niche_key Niche key (optional)
#' @param meta_key Metadata key (default: "coord")
#' @param group_by Cell type grouping column
#' @param group_use Specific groups to include
#' @param features Gene names to analyze
#' @param feature_colnm Metadata columns to analyze
#' @param mistyR_params Numeric vector of length 2: c(juxtaview_radius, paraview_radius)
#' @param coord_interval Spatial scaling coord_intervals
#' @param tmp_dir Temporary directory
#' @param output_dir Output directory
#' @param photo_dir Plot directory
#'
#' @return List of MSTIDyR results per sample
#'
#' @keywords internal
#' @noRd
.mistyR_pipeline <- function( STID_obj = NULL,
                              loop_single = NULL,
                              meta_key = NULL,
                              niche_key = NULL,
                              group_by = NULL,
                              group_use = NULL,
                              features = NULL,
                              feature_colnm = NULL,
                              mistyR_params = NULL,
                              coord_interval = NULL,
                              tmp_dir = NULL,
                              output_dir = NULL,
                              photo_dir = NULL
){
  if(length(mistyR_params) != 2) {
    clog_error("mistyR_params should be a numeric vector of length 2: c(juxtaview_radius, paraview_radius).")
  }
  results_list <- list()
  tmp_dir <- tempdir()
  for(i in seq_along(loop_single)){
    i_single <- loop_single[i]
    clog_loop(paste0("Processing samp_id: ", i_single, " (", i, "/", length(loop_single), ")"))
    logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
    i_coord_interval <- coord_interval[i]

    #> check
    if(!is.null(niche_key)){
      if(!logic_niche){
        clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
        next
      }
      clog_normal("Using niche_key ...")
      Niche_cells <- GetSSNicheCells(STID_obj = STID_obj, niche_key = niche_key, loop_id = i_single)[[1]]
      .check_column_exist(Niche_cells, "is_Niche")
      Niche_cells <- Niche_cells %>%
        filter(is_Niche) # only niche, not all
    }else{
      if(is.null(meta_key)){
        clog_warn("Both niche_key and meta_key are NULL, will use meta_key: coord for analysis...")
        meta_key <- "coord"
      }
      clog_normal("Using meta_key ...")
      Niche_cells <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]] %>%
        filter(!!sym(samp_colnm) == i_single) # actually all, not niche
    }
    if(nrow(Niche_cells) >5000){
      clog_warn("Number of niche cells is greater than 5000, will randomly sample
                5000 cells for CoLocalization calculation to save time.")
      Niche_cells <- Niche_cells %>% slice_sample(n = 5000, replace = FALSE)
    }
    if(!is.null(group_by)){
      .check_column_exist(Niche_cells, group_by)
      if(!is.null(group_use)){
        Niche_cells <- Niche_cells %>%
          mutate(!!sym(group_by) := as.character(!!sym(group_by))) %>%
          filter(!!sym(group_by) %in% group_use)
      }
      expr <- Niche_cells[,group_by,drop=F]
      expr <- model.matrix(as.formula(paste0("~",group_by," - 1")),data = expr) %>%
        as.data.frame()
    }else if(!is.null(features)){
      valid_features <- .check_features_exist(STID_obj, features)
      expr <- FetchData(STID_obj, vars = valid_features) %>% #
        as.data.frame()
      expr <- expr[rownames(Niche_cells), , drop=F]
      if(!is.null(feature_colnm)){
        .check_column_exist(Niche_cells, feature_colnm)
        expr2 <- Niche_cells %>%
          dplyr::select(all_of(feature_colnm))
        expr <- cbind(expr, expr2)
      }
    }else if(!is.null(feature_colnm)){
      .check_column_exist(Niche_cells, feature_colnm)
      expr <- Niche_cells %>%
        dplyr::select(all_of(feature_colnm))
    }

    # >
    colnames(expr) <- gsub(paste0("^", group_by), "", colnames(expr))
    clog_normal("Make syntactically valid names out of column names...")
    colnames(expr) <- make.names(colnames(expr))
    pos <- Niche_cells[,c("x","y")]

    #> mistyR
    clog_normal("Create mistyR intraview, juxtaview and paraview. It may take a while...")
    misty.intra <- mistyR::create_initial_view(expr)
    misty.views <- misty.intra %>%
      mistyR::add_juxtaview(pos, neighbor.thr = mistyR_params[1]*i_coord_interval) %>%  # add_juxtaview，阈值和半径过大会比较慢
      mistyR::add_paraview(pos, l = mistyR_params[2]*i_coord_interval)
    view_nms <- c("intra", paste0("juxta.", mistyR_params[1]), paste0("para.", mistyR_params[2]))
    results_list[[i_single]][["data"]][["misty_views"]] <- misty.views

    #>
    clog_normal("Run mistyR model, it may take a while...")
    misty.results <- misty.views %>%
      mistyR::run_misty(results.folder = paste0(tmp_dir,"/results"),
                        seed = 42,cv.folds = 10) %>%
      mistyR::collect_results()
    results_list[[i_single]][["data"]][["misty_results"]] <- misty.results

    #>
    clog_normal("Plot improvement stats ...")
    pdf(paste0(photo_dir,"/",i_single,"_MistyR_Improvement_Stats.pdf"),
        width = 6, height = 6)
    misty.results %>%
      mistyR::plot_improvement_stats("gain.R2") %>%
      mistyR::plot_improvement_stats("gain.RMSE")
    dev.off()

    #>
    clog_normal("Plot view contributions ...")
    pdf(paste0(photo_dir,"/",i_single,"_MistyR_View_Contributions.pdf"),
        width = 6, height = 6)
    misty.results %>%
      mistyR::plot_view_contributions()
    dev.off()

    #>
    clog_normal("Plot interaction heatmaps ...")
    pdf(paste0(photo_dir,"/",i_single,"_MistyR_Interaction_Heatmaps.pdf"),
        width = 6, height = 6)
    for(i in seq_along(view_nms)){
      view_nm <- view_nms[i]
      misty.results %>%
        .plot_interaction_heatmap(view_nm, cutoff = 0) # heatmap_cutoff = 0
    }
    dev.off()

    #> compare
    clog_normal("Plot interaction comparison heatmaps ...")
    pdf(paste0(photo_dir,"/",i_single,"_MistyR_Interaction_Comparison_Heatmaps.pdf"),
        width = 6, height = 6)
    for (i in 2:length(view_nms)) {
      view_nm1 <- view_nms[1]
      view_nm2 <- view_nms[i]
      misty.results %>%
        .plot_contrast_heatmap(view_nm1, view_nm2, cutoff = 0) # heatmap_cutoff = 0
    }
    dev.off()

    #> network
    clog_normal("Plot interaction communities ...")
    pdf(paste0(photo_dir,"/",i_single,"_MistyR_Interaction_Communities.pdf"),
        width = 6, height = 6)
    for(i in seq_along(view_nms)){
      view_nm <- view_nms[i]
      misty.results %>%
        mistyR::plot_interaction_communities(view_nm, cutoff = 1) # comm_cutoff = 1
    }
    dev.off()
  }
  return(results_list)
}


#' Plot Interaction Heatmap
#'
#' Internal function to create interaction heatmaps from MSTIDyR results.
#'
#' @param misty.results MSTIDyR results object
#' @param view Character, view name
#' @param cutoff Numeric, importance cutoff
#' @param trim Numeric, trim threshold
#' @param trim.measure Character, measure to trim by
#' @param clean Logical, whether to clean data
#'
#' @return Modified misty.results (invisible)
#'
#' @keywords internal
#'
#' @noRd
.plot_interaction_heatmap <- function(misty.results, view, cutoff = 1, trim = -Inf,
                                      trim.measure = c("gain.R2", "multi.R2", "intra.R2",
                                                       "gain.RMSE", "multi.RMSE", "intra.RMSE"),
                                      clean = FALSE) {
  trim.measure.type <- match.arg(trim.measure)
  assertthat::assert_that(
    "importances.aggregated" %in% names(misty.results),
    msg = "The provided result list is malformed. Consider using collect_results()."
  )
  assertthat::assert_that(
    "improvements.stats" %in% names(misty.results),
    msg = "The provided result list is malformed. Consider using collect_results()."
  )
  available.views <- misty.results$importances.aggregated %>%
    dplyr::pull(.data$view) %>%
    unique()
  assertthat::assert_that(
    view %in% available.views,
    msg = "The selected view cannot be found in the results table."
  )
  inv <- sign(
    (stringr::str_detect(trim.measure.type, "gain") |
       stringr::str_detect(trim.measure.type, "RMSE", negate = TRUE)) - 0.5
  )
  targets <- misty.results$improvements.stats %>%
    dplyr::filter(
      .data$measure == trim.measure.type,
      inv * .data$mean >= inv * trim
    ) %>%
    dplyr::pull(.data$target)
  plot.data <- misty.results$importances.aggregated %>%
    dplyr::filter(
      .data$view == !!view,
      .data$Target %in% targets
    )
  if (clean) {
    clean.predictors <- plot.data %>%
      dplyr::mutate(
        Importance = .data$Importance * (.data$Importance >= cutoff)
      ) %>%
      dplyr::group_by(.data$Predictor) %>%
      dplyr::summarize(
        total = sum(.data$Importance, na.rm = TRUE)
      ) %>%
      dplyr::filter(.data$total > 0) %>%
      dplyr::pull(.data$Predictor)
    clean.targets <- plot.data %>%
      dplyr::mutate(
        Importance = .data$Importance * (.data$Importance >= cutoff)
      ) %>%
      dplyr::group_by(.data$Target) %>%
      dplyr::summarize(
        total = sum(.data$Importance, na.rm = TRUE)
      ) %>%
      dplyr::filter(.data$total > 0) %>%
      dplyr::pull(.data$Target)
    plot.data.clean <- plot.data %>%
      dplyr::filter(
        .data$Predictor %in% clean.predictors,
        .data$Target %in% clean.targets
      )
  } else {
    plot.data.clean <- plot.data
  }
  set2.blue <- "#1D91C0"
  results.plot <- ggplot2::ggplot(
    plot.data.clean,
    ggplot2::aes(x = .data$Predictor, y = .data$Target)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data$Importance),
      color = "gray50", linewidth = 0.2
    ) +
    ggplot2::scale_fill_gradient2(
      low = "white",
      mid = "white",
      high = set2.blue,
      na.value = "white",
      midpoint = cutoff
    ) +
    ggplot2::theme_test() +  # 使用更简洁的主题
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      panel.border = ggplot2::element_blank(),  # 移除面板边框
      plot.background = ggplot2::element_blank(), # 移除绘图区域背景边框
      plot.title = ggplot2::element_text(hjust = 0.5)  # 标题居中
    ) +
    ggplot2::coord_equal() +
    ggplot2::ggtitle(view)
  print(results.plot)
  invisible(misty.results)
}

#' Plot Contrast Heatmap
#'
#' Internal function to create contrast heatmaps comparing two views from MSTIDyR results.
#'
#' @param misty.results MSTIDyR results object
#' @param from.view Character, source view name
#' @param to.view Character, target view name
#' @param cutoff Numeric, importance cutoff
#' @param trim Numeric, trim threshold
#' @param trim.measure Character, measure to trim by
#'
#' @return Modified misty.results (invisible)
#'
#' @keywords internal
#'
#' @noRd
.plot_contrast_heatmap <- function(misty.results, from.view, to.view, cutoff = 1, trim = -Inf,
                                   trim.measure = c("gain.R2", "multi.R2", "intra.R2",
                                                    "gain.RMSE", "multi.RMSE", "intra.RMSE")) {
  trim.measure.type <- match.arg(trim.measure)
  assertthat::assert_that(
    "importances.aggregated" %in% names(misty.results),
    msg = "The provided result list is malformed. Consider using collect_results()."
  )
  assertthat::assert_that(
    "improvements.stats" %in% names(misty.results),
    msg = "The provided result list is malformed. Consider using collect_results()."
  )
  available.views <- misty.results$importances.aggregated %>%
    dplyr::pull(.data$view) %>% unique()

  assertthat::assert_that(
    from.view %in% available.views,
    msg = "The selected from.view cannot be found in the results table."
  )
  assertthat::assert_that(
    to.view %in% available.views,
    msg = "The selected to.view cannot be found in the results table."
  )
  inv <- sign(
    (stringr::str_detect(trim.measure.type, "gain") |
       stringr::str_detect(trim.measure.type, "RMSE", negate = TRUE)) - 0.5
  )
  targets <- misty.results$improvements.stats %>%
    dplyr::filter(
      .data$measure == trim.measure.type,
      inv * .data$mean >= inv * trim
    ) %>%
    dplyr::pull(.data$target)
  from.view.wide <- misty.results$importances.aggregated %>%
    dplyr::filter(
      .data$view == from.view,
      .data$Target %in% targets
    ) %>%
    tidyr::pivot_wider(
      names_from = "Target",
      values_from = "Importance",
      -c(.data$view, .data$nsamples)
    )
  to.view.wide <- misty.results$importances.aggregated %>%
    dplyr::filter(
      .data$view == to.view,
      .data$Target %in% targets
    ) %>%
    tidyr::pivot_wider(
      names_from = "Target",
      values_from = "Importance",
      -c(.data$view, .data$nsamples)
    )
  mask <- ((from.view.wide %>% dplyr::select(-.data$Predictor)) < cutoff) &
    ((to.view.wide %>% dplyr::select(-.data$Predictor)) >= cutoff)
  masked <- (to.view.wide %>% tibble::column_to_rownames("Predictor")) * mask
  plot.data <- masked %>%
    dplyr::slice(which(masked %>% rowSums(na.rm = TRUE) > 0)) %>%
    dplyr::select(which(masked %>% colSums(na.rm = TRUE) > 0)) %>%
    tibble::rownames_to_column("Predictor") %>%
    tidyr::pivot_longer(
      names_to = "Target",
      values_to = "Importance",
      -.data$Predictor
    )
  set2.blue <- "#1D91C0"
  results.plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x = .data$Predictor, y = .data$Target)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$Importance),
                       color = "gray50",linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      low = "white",
      mid = "white",
      high = set2.blue,
      na.value = "white",
      midpoint = cutoff
    ) +
    ggplot2::theme_test() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      panel.border = ggplot2::element_blank(),  # 移除面板边框
      plot.background = ggplot2::element_blank(), # 移除绘图区域背景边框
      plot.title = ggplot2::element_text(hjust = 0.5)  # 标题居中
    ) +
    ggplot2::coord_equal() +
    ggplot2::ggtitle(paste0(to.view, " - ", from.view))
  print(results.plot)
  invisible(misty.results)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot_SpatialCoLoc: noMS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Plot Spatial Co-localization
#'
#' Visualizes spatial co-localization of two cell groups, genes, or metadata
#' features in selected samples or niches.
#'
#' @param STID_obj An STID object.
#' @param loop_id Character, sample loop identifier
#' (default: \code{"LoopAllSamp"}).
#' @param meta_key Character, metadata key used when \code{niche_key} is
#' \code{NULL}.
#' @param niche_key Character, niche key used to retrieve niche cells.
#' @param group_by Character, metadata column for group-based co-localization.
#' @param group_use Character vector, groups to keep from \code{group_by}.
#' @param features Character vector, feature or gene names for co-localization.
#' @param feature_colnm Character vector, metadata columns used as features.
#' @param exp_thres Numeric, threshold below which feature values are set to zero
#' (default: \code{1}).
#' @param col Character vector, color palette
#' (default: \code{COLOR_LIST$PALETTE_WHITE_BG}).
#' @param pt_size Numeric, point size for spatial plots
#' (default: \code{1.5}).
#'
#' @return Invisibly returns \code{NULL}; spatial plots are printed.
#'
#' @importFrom Seurat FetchData
#' @importFrom dplyr filter mutate if_else select
#' @importFrom rlang sym
#' @importFrom tidyselect all_of
#' @importFrom stats as.formula model.matrix
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Plot_SpatialCoLoc(
#' STID_obj = STID_obj,
#' niche_key = "niche_virulence",
#' group_by = "cell_type",
#' group_use = c("Macrophage", "T cell")
#' )
#'
#' Plot_SpatialCoLoc(
#' STID_obj = STID_obj,
#' niche_key = "niche_virulence",
#' features = c("Cxcl10", "Isg15")
#' )
#' }
Plot_SpatialCoLoc <- function(STID_obj = NULL,
                              loop_id = "LoopAllSamp",
                              meta_key = NULL,
                              niche_key = NULL,
                              group_by = NULL,
                              group_use = NULL,
                              features = NULL,
                              feature_colnm = NULL,
                              exp_thres = 1,
                              col = COLOR_LIST$PALETTE_WHITE_BG,
                              pt_size = 1.5
){
  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(niche_key,col)
  if(.all_null(group_by,features,feature_colnm)){
    clog_warn("group_by, features and feature_colnm are all NULL, will use the default celltype_colnm from STID object as group_by.")
    group_by <- GetInfo(STID_obj, info_key = "data_info",sub_key = "celltype_colnm")[[1]]
  }
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)

  # >>> Start main pipeline
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  for(i in seq_along(loop_single)){
    i_single <- loop_single[i]
    clog_loop(paste0("Processing samp_id: ", i_single, " (", i, "/", length(loop_single), ")"))

    #> check
    if(!is.null(niche_key)){
      logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
      if(!logic_niche){
        clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
        next
      }
      clog_normal("Using niche_key...")
      Niche_cells <- GetSSNicheCells(STID_obj = STID_obj, niche_key = niche_key, loop_id = i_single)[[1]]
      .check_column_exist(Niche_cells, "is_Niche")
      Niche_cells <- Niche_cells %>%
        filter(is_Niche) # only niche, not all
    }else{
      if(is.null(meta_key)){
        clog_warn("Both niche_key and meta_key are NULL, will use meta_key: coord for analysis...")
        meta_key <- "coord"
      }
      clog_normal("Using meta_key...")
      Niche_cells <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]] %>%
        filter(!!sym(samp_colnm) == i_single) # actually all, not niche
    }
    if(!is.null(group_by)){
      .check_column_exist(Niche_cells, group_by)
      if(!is.null(group_use)){
        Niche_cells <- Niche_cells %>%
          # filter(!!sym(group_by) %in% group_use) %>%
          mutate(!!sym(group_by) := if_else(!!sym(group_by) %in% group_use, as.character(!!sym(group_by)), NA_character_))

      }
      expr <- Niche_cells[,group_by,drop=F]
      expr <- model.matrix(as.formula(paste0("~",group_by," - 1")),data = expr) %>%
        as.data.frame()
    }else if(!is.null(features)){
      valid_features <- .check_features_exist(STID_obj, features)
      expr <- FetchData(STID_obj, vars = valid_features) %>% #
        as.data.frame()
      expr <- expr[rownames(Niche_cells), , drop=F]
      if(!is.null(feature_colnm)){
        .check_column_exist(Niche_cells, feature_colnm)
        expr2 <- Niche_cells %>%
          dplyr::select(all_of(feature_colnm))
        expr <- cbind(expr, expr2)
      }
      expr[which(expr < exp_thres, arr.ind = TRUE)] <- 0
    }else if(!is.null(feature_colnm)){
      .check_column_exist(Niche_cells, feature_colnm)
      expr <- Niche_cells %>%
        dplyr::select(all_of(feature_colnm))
      expr[which(expr < exp_thres, arr.ind = TRUE)] <- 0
    }
    if(ncol(expr) >2){
      clog_error("The variables number for co-localization analysis should be 2.")
    }

    # >
    if(!is.null(group_by)){
      clog_normal("Plot group co-localization...")
      plot_data <- Niche_cells[c("x","y", group_by)]
      clog_normal("The number of cells in each group:")
      print(table(plot_data[[group_by]]))
      p1 <- Plot_Spatial(plot_data = plot_data, x_colnm = "x", y_colnm = "y", group_by = group_by,
                         facet_grpnm = NULL, datatype = "discrete",
                         col = col,
                         pt_size = pt_size,vmin = NULL, vmax = "p99",
                         title = i_single, subtitle = NULL,black_bg = F)
      print(p1)
    }else{
      clog_normal("Plot feature co-localization...")
      plot_data <- data.frame(expr, x = Niche_cells$x, y = Niche_cells$y)
      plot_data$celltype <- if_else(expr[,1] >0 & expr[,2] >0, "Both",
                                    if_else(expr[,1] >0, colnames(expr)[1],
                                            if_else(expr[,2] >0, colnames(expr)[2], NA_character_)))
      clog_normal("The number of cells in each group:")
      print(table(plot_data$celltype))
      p1 <- Plot_Spatial(plot_data = plot_data, x_colnm = "x", y_colnm = "y", group_by = "celltype",
                         facet_grpnm = NULL, datatype = "discrete",
                         col = col,
                         pt_size = pt_size,vmin = NULL, vmax = "p99",
                         title = i_single, subtitle = NULL,black_bg = F)
      print(p1)
    }
  }

  clog_end()
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalSampDEGs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate Niche Differential Expression Genes
#'
#' Identifies differentially expressed genes between niche and bystander cells,
#' both globally and per cell type. Generates volcano plots and faceted volcano
#' plots for visualization.
#'
#' @param STID_obj An STID object containing niche analysis results
#' @param samp_mode  Character, sample type - "SS" (single sample) or "MS" (multi-sample)
#'        (default: "SS")
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param samp_grp_index Logical, whether to group by sample groups in MS mode
#'        (default: FALSE)
#' @param meta_key Character, metadata key for MS mode when niche_key is NULL
#' @param niche_key Character, niche key to analyze
#' @param group_by Character, column name for cell type grouping
#' @param group_value Character vector, specific cell types to analyze
#' @param assay_id Character, assay name (default: "Spatial")
#' @param test_method Character, statistical test to use (default: "wilcox")
#' @param logfc_thres Numeric, log2 fold change threshold (default: 1)
#' @param min_pct Numeric, minimum percentage of cells expressing gene (default: 0.01)
#' @param padj_thres Numeric, adjusted p-value threshold (default: 0.05)
#' @param adjust_method Character, p-value adjustment method (default: "BH")
#' @param topGeneN Integer, number of top genes to label (default: 3)
#' @param col Color palette for visualization (default: COLOR_LIST$PALETTE_WHITE_BG)
#' @param remove_genes Character vector, genes to exclude from analysis
#' @param return_data Logical, whether to return results list (default: TRUE)
#' @param grp_nm Character, group name for output organization (default: NULL)
#' @param dir_nm Character, directory name for output (default: "M3_CalSampDEGs")
#' @param ... Additional parameters for FindAllMarkers and FindMarkers functions
#'
#' @return If return_data = TRUE, returns a list of DEG results per sample
#'
#' @import ggrepel
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Find DEGs between niche and bystander cells
#' results <- CalSampDEGs(
#'   STID_obj = STID_obj,
#'   niche_key = "niche_virulence",
#'   group_by = "cell_type",
#'   logfc_thres = 0.5,
#'   padj_thres = 0.05,
#'   topGeneN = 5
#' )
#' }
CalSampDEGs <- function(STID_obj = NULL,
                        samp_mode  = "SS",
                        loop_id = "LoopAllSamp",
                        samp_grp_index = FALSE,
                        meta_key = NULL,
                        niche_key = NULL,
                        group_by = NULL,
                        group_value = NULL,
                        assay_id = "Spatial",
                        test_method = "wilcox",
                        logfc_thres = 1,
                        min_pct = 0.01,
                        padj_thres = 0.05,
                        adjust_method = "BH",
                        topGeneN = 5,
                        remove_genes = NULL,
                        col = COLOR_LIST[["PALETTE_WHITE_BG"]],
                        return_data = TRUE,
                        grp_nm = NULL, dir_nm = "M3_CalSampDEGs",
                        ...
){
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
  .check_null_args(col)
  match.arg(samp_mode, .STID_globals$samp_modes)
  if(is.null(group_by)){
    clog_warn("group_by is NULL, will use the default celltype_colnm from STID object.")
    group_by <- GetInfo(STID_obj, info_key = "data_info",sub_key = "celltype_colnm")[[1]]
  }
  valid_genes <- rownames(STID_obj)[!rownames(STID_obj) %in% remove_genes]
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  results_list <- list()
  clog_step("Start calculating DEGs...")
  if(samp_mode  == "SS"){
    clog_normal("Execute single-sample analysis...")
    .check_null_args(niche_key);.check_one_arg(niche_key)
    loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)

    #>
    for(i in seq_along(loop_single)){
      i_single <- loop_single[i]
      clog_loop(paste0("Processing samp_id: ", i_single, " (", i, "/", length(loop_single), ")"))
      logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
      if(!logic_niche){
        clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
        next
      }

      #>
      Niche_cells <- GetSSNicheCells(STID_obj = STID_obj, niche_key = niche_key,loop_id = i_single)[[1]]
      .check_column_exist(Niche_cells, "is_Niche")
      Niche_cells <- Niche_cells %>%
        mutate(Niche_label = factor(ifelse(is_Niche,"Niche","Bystander"),levels = c("Bystander","Niche"))) %>%
        mutate(Niche_celltype = paste0(Niche_label,"_",!!sym(group_by)))

      #>
      samp_seurat <- subset(STID_obj, subset = !!sym(samp_colnm) == i_single,features = valid_genes)
      if(all(colnames(samp_seurat) == rownames(Niche_cells))){
        samp_seurat@meta.data <- Niche_cells
      }
      valid_Niche_celltypes <- Niche_cells$Niche_celltype %>% na.omit() %>% as.character() %>% unique() %>% sort()

      #> DEGs overall
      clog_title("Calculating DEGs between Niche and Bystander cells, it may take a while...")
      degs_niche <- FindMarkers(object = samp_seurat,
                                ident.1 = "Niche",
                                ident.2 = "Bystander",
                                group.by = "Niche_label",
                                assay = assay_id,
                                min.pct = min_pct,
                                test.use = test_method,
                                only.pos = FALSE,
                                ...)
      degs_niche_list <- .Process_DEGs(degs_niche, grp_id = "Overall", logfc_thres = logfc_thres,
                                       padj_thres = padj_thres, topGeneN = topGeneN)
      degs_niche <- degs_niche_list$degs
      volcano_title <- degs_niche_list$volcano_title
      print(table(degs_niche["change"]))
      write.table(degs_niche, file = paste0(output_dir,"/",i_single,"_Niche_vs_Bystander_DEGs.txt"),
                  sep = "\t", quote = F, row.names = T)


      #> DEGs per celltype
      clog_title("Calculating DEGs between Niche and Bystander cells for each celltype...")
      valid_celltypes <- Niche_cells[[group_by]] %>% na.omit() %>% as.character() %>% unique() %>% sort()
      if(!is.null(group_value)){
        valid_celltypes <- valid_celltypes[na.omit(match(group_value, valid_celltypes))]
      }
      degs_niche_celltypes <- list()
      for(j in seq_along(valid_celltypes)){
        j_celltype <- valid_celltypes[j]
        ident_1 <- paste0("Niche_",j_celltype)
        ident_2 <- paste0("Bystander_",j_celltype)
        clog_normal(paste0("Calculating DEGs for celltype: ", j_celltype, " (", j, "/", length(valid_celltypes), ")"))
        min_value <- table(Niche_cells$Niche_celltype[Niche_cells$Niche_celltype %in% c(ident_1,ident_2)]) %>% min()
        if(min_value < 3){
          clog_warn(paste0("Skipping celltype ", j_celltype, " due to insufficient cell numbers (<3)."))
          next
        }
        if(!all(c(ident_1,ident_2) %in% valid_Niche_celltypes)){
          clog_warn(paste0("Skipping celltype ", j_celltype, " due to lack of Niche or Bystander cells."))
          next
        }
        degs_niche_ct <- FindMarkers(object = samp_seurat,
                                     ident.1 = paste0("Niche_",j_celltype),
                                     ident.2 = paste0("Bystander_",j_celltype),
                                     group.by = "Niche_celltype",
                                     assay = assay_id,
                                     min.pct = min_pct,
                                     test.use = test_method,
                                     only.pos = FALSE,
                                     ...)
        degs_niche_ct <- .Process_DEGs(degs_niche_ct, grp_id = j_celltype, logfc_thres = logfc_thres,
                                       padj_thres = padj_thres, topGeneN = topGeneN)$degs
        degs_niche_celltypes[[j_celltype]] <- degs_niche_ct
      }
      degs_niche_ct <- bind_rows(degs_niche_celltypes)
      print(table(degs_niche_ct[c("group","change")]))
      write.table(degs_niche_ct, file = paste0(output_dir,"/",i_single,"_Niche_vs_Bystander_DEGs_per_Celltype.txt"),
                  sep = "\t", quote = F, row.names = T)

      #> plot
      clog_title("Plotting the overall niche DEGs volcano plot...")
      p1 <- ggplot(data = degs_niche,aes(x = avg_log2FC_cap, y = -log10(p_val_adj_cap),color = change)) +
        geom_point(alpha = 0.6, size = 1.5,shape = 16,stroke = 0) +
        theme_bw() +
        scale_y_continuous(expand = c(0.06, 0,0.15,0)) +
        geom_hline(yintercept=-log10(padj_thres) ,linetype=2,color = "grey30") +
        geom_vline(xintercept=c(logfc_thres,-logfc_thres) ,linetype=2 ,color = "grey30") +
        labs(x="log2(Fold Change)",y='-log10(adj P value)',title = paste0("Niche DEGs Volcano Plot (", i_single,")"),
             subtitle=volcano_title,color='Group')+
        scale_color_manual(values = c('#a121f0','grey80','#ffad21'))+
        geom_text_repel(aes(label=sign),  #显示感兴趣基因
                        size = 4,
                        fontface="italic",
                        color="grey50",
                        box.padding=unit(0.35, "lines"), #文本框周边填充
                        point.padding=unit(0.5, "lines"), #点周边填充
                        segment.colour = "grey50",  #连接点与标签的线段的颜色
                        min.segment.length = 0.1,
                        max.overlaps=100000)+ #最多点的数量
        guides(color=guide_legend(override.aes=list(size=4)))+
        theme_common()
      ggsave(filename = paste0(photo_dir,"/",i_single,"_Niche_vs_Bystander_DEGs_Volcano.pdf"),
             plot = p1, width = 6, height = 5)
      print(p1)

      # per celltype
      clog_title("Plotting the niche DEGs volcano plot for each celltype...")
      p2 <- .Plot_scVolcano(
        diffData = degs_niche_ct,
        group_by = "group",
        topGeneN = 3, col = col,
        orderBy = "avg_log2FC",
        plotTitle = paste0("Niche DEGs per Celltype (", i_single,")"),
        log2FC_thres = logfc_thres,
        backH = 0.1,
        myMarkers = NULL,
        clusterOrder = valid_celltypes
      )
      ggsave(filename = paste0(photo_dir,"/",i_single,"_Niche_vs_Bystander_DEGs_per_Celltype_Volcano.pdf"),
             plot = p2, width = 6, height = 5)
      print(p2)
      results_list[[i_single]] <- list(
        data = list(Overall_DEGs = degs_niche,
                    Celltype_DEGs = degs_niche_ct),
        plots = list(Overall_Volcano = p1,
                     Celltype_Volcano = p2)
      )
    }
  }else if(samp_mode  == "MS"){
    clog_normal("Execute multi-sample analysis...")
    loop_multi <- .check_loop_multi(STID_obj = STID_obj, loop_id = loop_id)
    for(i in seq_along(loop_multi)){
      i_multi <- loop_multi[i]
      clog_loop(paste0("Processing samp_id: ", i_multi, " (", i, "/", length(loop_multi), ")"))
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
          clog_warn("Both niche_key and meta_key are NULL, will use meta_key: coord for analysis...")
          meta_key <- "coord"
        }
        clog_normal("Using meta_key for analysis...")
        Niche_cells <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]] %>%
          filter(!!sym(samp_colnm) %in% samp_id2grp$samp_id) # actually all, not niche
      }
      .check_column_exist(Niche_cells, group_by)
      if(samp_grp_index){
        Niche_cells <- Niche_cells %>%
          mutate(samp_grp = samp_id2grp$samp_grp[match(Niche_cells[[samp_colnm]],samp_id2grp$samp_id)],.after = samp_colnm)
        new_samp_colnm <- "samp_grp"
      }else{
        new_samp_colnm <- samp_colnm
      }
      if(!is.null(niche_key)){
        Niche_cells <- Niche_cells %>%
          mutate(Niche_label = ifelse(is_Niche,"Niche","Bystander"),.after = "is_Niche") %>%
          mutate(Niche_label = factor(Niche_label,levels = c("Bystander","Niche")))
      }else{
        Niche_cells <- Niche_cells %>%
          mutate(Niche_label = "All",.after = all_of(group_by)) %>%
          mutate(Niche_label = factor(Niche_label,levels = c("All")))
      }
      Niche_cells <- Niche_cells %>%
        mutate(Niche_celltype = paste0(Niche_label,"_",!!sym(group_by)))

      #>
      samp_seurat <- subset(STID_obj, subset = !!sym(samp_colnm) %in% samp_id2grp$samp_id,features = valid_genes)
      if(all(colnames(samp_seurat) == rownames(Niche_cells))){
        samp_seurat@meta.data <- Niche_cells
      }
      valid_Niche_celltypes <- Niche_cells$Niche_celltype %>% na.omit() %>% as.character() %>% unique() %>% sort()

      #>
      valid_regions <- Niche_cells$Niche_label %>% unique() %>% sort() %>% as.character()
      valid_samps <- Niche_cells[[new_samp_colnm]] %>% unique() %>% sort() %>% as.character()
      for (j in seq_along(valid_regions)) {
        j_region <- valid_regions[j]
        clog_normal(paste0("Processing region: ", j_region, " (", j, "/", length(valid_regions), ")"))
        j_samp_seurat <- subset(samp_seurat, subset = Niche_label == j_region)

        #> DEGs overall
        clog_title("Calculating DEGs between comparative groups, it may take a while...")
        clog_normal(paste0("Comparing groups: ", valid_samps[2], " vs ", valid_samps[1]))
        degs_niche <- FindMarkers(object = j_samp_seurat,
                                  ident.1 = valid_samps[2], # case
                                  ident.2 = valid_samps[1],
                                  group.by = new_samp_colnm,
                                  assay = assay_id,
                                  min.pct = min_pct,
                                  test.use = test_method,
                                  only.pos = FALSE,
                                  ...)
        degs_niche_list <- .Process_DEGs(degs_niche, grp_id = j_region, logfc_thres = logfc_thres,
                                         padj_thres = padj_thres, topGeneN = topGeneN)
        degs_niche <- degs_niche_list$degs
        volcano_title <- degs_niche_list$volcano_title
        print(table(degs_niche["change"]))
        write.table(degs_niche, file = paste0(output_dir,"/",i_multi,"_",j_region,"_DEGs.txt"),
                    sep = "\t", quote = F, row.names = T)

        #> DEGs per celltype
        clog_title("Calculating DEGs between comparative groups for each celltype...")
        valid_celltypes <- Niche_cells[[group_by]] %>% na.omit() %>% as.character() %>% unique() %>% sort()
        if(!is.null(group_value)){
          valid_celltypes <- valid_celltypes[na.omit(match(group_value, valid_celltypes))]
        }
        degs_niche_celltypes <- list()
        for(k in seq_along(valid_celltypes)){
          k_celltype <- valid_celltypes[k]
          clog_normal(paste0("Calculating DEGs for celltype: ", k_celltype, " (", k, "/", length(valid_celltypes), ")"))
          k_samp_seurat <- subset(j_samp_seurat, subset = !!sym(group_by) == k_celltype)
          samp_table <- table(k_samp_seurat@meta.data[[samp_colnm]]) %>% as.data.frame()
          if(min(samp_table$Freq) < 3 | nrow(samp_table) < 2){
            clog_warn(paste0("Skipping celltype ", k_celltype, " due to insufficient cell numbers (<3)."))
            next
          }
          degs_niche_ct <- FindMarkers(object = k_samp_seurat,
                                       ident.1 = valid_samps[2], # case
                                       ident.2 = valid_samps[1],
                                       group.by = new_samp_colnm,
                                       assay = assay_id,
                                       min.pct = min_pct,
                                       test.use = test_method,
                                       only.pos = FALSE,
                                       ...)
          degs_niche_ct <- .Process_DEGs(degs_niche_ct, grp_id = k_celltype, logfc_thres = logfc_thres,
                                         padj_thres = padj_thres, topGeneN = topGeneN)$degs
          degs_niche_celltypes[[k_celltype]] <- degs_niche_ct
        }
        degs_niche_ct <- bind_rows(degs_niche_celltypes)
        print(table(degs_niche_ct[c("group","change")]))
        write.table(degs_niche_ct, file = paste0(output_dir,"/",i_multi,"_",j_region,"_DEGs_per_Celltype.txt"),
                    sep = "\t", quote = F, row.names = T)

        #> plot
        clog_title("Plotting the overall niche DEGs volcano plot...")
        plot_id <- paste0(i_multi,"_",j_region)
        p1 <- ggplot(data = degs_niche,aes(x = avg_log2FC_cap, y = -log10(p_val_adj_cap),color = change)) +
          geom_point(alpha = 0.6, size = 1.5,shape = 16,stroke = 0) +
          theme_bw() +
          scale_y_continuous(expand = c(0.06, 0,0.15,0)) +
          geom_hline(yintercept=-log10(padj_thres) ,linetype=2,color = "grey30") +
          geom_vline(xintercept=c(logfc_thres,-logfc_thres) ,linetype=2 ,color = "grey30") +
          labs(x="log2(Fold Change)",y='-log10(adj P value)',title = paste0("Niche DEGs Volcano Plot (",plot_id,")"),
               subtitle=volcano_title,color='Group')+
          scale_color_manual(values = c('#a121f0','grey80','#ffad21'))+
          geom_text_repel(aes(label=sign),  #显示感兴趣基因
                          size = 4,
                          fontface="italic",
                          color="grey50",
                          box.padding=unit(0.35, "lines"), #文本框周边填充
                          point.padding=unit(0.5, "lines"), #点周边填充
                          segment.colour = "grey50",  #连接点与标签的线段的颜色
                          min.segment.length = 0.1,
                          max.overlaps=100000)+ #最多点的数量
          guides(color=guide_legend(override.aes=list(size=4)))+
          theme_common()
        ggsave(filename = paste0(photo_dir,"/",plot_id,"_DEGs_Volcano.pdf"),
               plot = p1, width = 6, height = 5)
        print(p1)

        # per celltype
        clog_title("Plotting the niche DEGs volcano plot for each celltype...")
        p2 <- .Plot_scVolcano(
          diffData = degs_niche_ct,
          group_by = "group",
          topGeneN = 3, col = col,
          orderBy = "avg_log2FC",
          plotTitle = paste0("Niche DEGs per Celltype (",plot_id,")"),
          log2FC_thres = logfc_thres,
          backH = 0.1,
          myMarkers = NULL,
          clusterOrder = valid_celltypes
        )
        ggsave(filename = paste0(photo_dir,"/",plot_id,"_DEGs_per_Celltype_Volcano.pdf"),
               plot = p2, width = 6, height = 5)
        print(p2)

        results_list[[plot_id]] <- list(
          # plots = list(Overall_Volcano = p1, # 消耗存储
          #              Celltype_Volcano = p2),
          data = list(Overall_DEGs = degs_niche,
                      Celltype_DEGs = degs_niche_ct)
        )
      }
    }
  }
  clog_normal("Saving results to RDS file...")
  saveRDS(results_list, file = paste0(output_dir,"/CalSampDEGs_Results_List.rds"))

  # >>> Final
  .save_function_params("CalSampDEGs", envir = environment(), file = paste0(output_dir,"Log_function_params_(CalSampDEGs).log"))
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(CalSampDEGs).log")) %>% invisible()
  if(return_data){
    return(results_list)
  }
}

#' Process Differential Expression Results
#'
#' Internal function to process DEG results, add significance labels,
#' and prepare for visualization.
#'
#' @param degs Data frame, differential expression results
#' @param grp_id Character, group identifier
#' @param padj_thres Numeric, adjusted p-value threshold
#' @param logfc_thres Numeric, log2 fold change threshold
#' @param topGeneN Integer, number of top genes to label
#'
#' @return List containing processed DEGs and volcano plot title
#'
#' @keywords internal
#'
#' @noRd
.Process_DEGs <- function(degs,grp_id, padj_thres = 0.05, logfc_thres = 0.25, topGeneN = 10){
  degs <- degs %>%
    mutate(.,gene = rownames(.),
           change = ifelse(.$p_val_adj < padj_thres & abs(.$avg_log2FC) > logfc_thres,
                           ifelse(.$avg_log2FC > logfc_thres, 'Up','Down'), 'No-Change'),
           sign = ifelse(.$p_val_adj < padj_thres & abs(.$avg_log2FC) > 2,rownames(.),NA),
           pct_mean = (pct.1 + pct.2)/2,
           pct_mean = ifelse(avg_log2FC<0,-pct_mean,pct_mean)) %>%
    arrange(desc(avg_log2FC)) %>%
    mutate(group = grp_id) %>%
    `rownames<-`(paste0(grp_id,"_",rownames(.)))

  degs$sign <- ""
  degs[rownames(degs)[1:topGeneN],"sign"] <- degs$gene[1:topGeneN]
  tmp_index <- nrow(degs)-topGeneN + 1
  degs[rownames(degs)[(tmp_index):nrow(degs)],"sign"] <- degs$gene[(tmp_index):nrow(degs)]
  degs <- degs %>%
    mutate(p_val_adj_cap = ifelse(p_val_adj < 10^-100, 10^-100, p_val_adj)) %>%
    mutate(avg_log2FC_cap = ifelse(avg_log2FC > 10, 10, ifelse(avg_log2FC < -10, -10, avg_log2FC)))

  volcano_title <- paste0('Cutoff for log2(FC) is ',round(logfc_thres,2),
                          '\nThe number of up gene is ',nrow(degs[degs$change == 'Up',]),
                          ' ,down gene is ',nrow(degs[degs$change == 'Down',]))
  return(list(degs = degs, volcano_title = volcano_title))
}


#' Plot Single-cell Volcano Plot
#'
#' Internal function to create faceted volcano plots for multiple groups.
#'
#' @param diffData Data frame, differential expression results
#' @param group_by Character, column name for grouping
#' @param topGeneN Integer, number of top genes to label
#' @param col Color palette
#' @param orderBy Character, column to order by
#' @param plotTitle Character, plot title
#' @param log2FC_thres Numeric, log2 fold change threshold
#' @param backH Numeric, background height
#' @param myMarkers Character vector, specific markers to highlight
#' @param clusterOrder Character vector, order of clusters
#'
#' @return A ggplot object
#'
#' @import ggprism
#'
#' @keywords internal
#'
#' @noRd
.Plot_scVolcano <- function(diffData, group_by = NULL,
                            topGeneN = 3, col = NULL,
                            orderBy="avg_log2FC_cap", plotTitle = NULL,
                            log2FC_thres = 4, backH=0.1,
                            myMarkers=NULL, clusterOrder=NULL
){

  diff.marker <- diffData %>%
    filter(abs(avg_log2FC_cap)>=log2FC_thres ) %>%
    filter(p_val_adj_cap < 0.05) %>%
    mutate(type = ifelse(avg_log2FC_cap>0, "Up", "Down")) %>%
    dplyr::rename(cluster = !!sym(group_by))

  cluster_values <- diff.marker$cluster
  if(is.numeric(cluster_values)){
    diff.marker <- diff.marker %>%
      mutate(cluster = paste0("Cluster", cluster_values)) %>%
      mutate(cluster = factor(cluster,levels = paste0("Cluster", sort(unique(cluster_values)))))
  }

  if (!is.null(clusterOrder)) {
    diff.marker$cluster <- factor(diff.marker$cluster, levels = clusterOrder)
  }

  back.data <- purrr::map_df(unique(diff.marker$cluster), function(x){
    tmp <- diff.marker %>% filter(cluster == x)
    new.tmp <- data.frame(
      cluster = x,
      min = min(tmp$avg_log2FC_cap)*1.1 - backH,
      max = max(tmp$avg_log2FC_cap)*1.1 + backH
    )
    return(new.tmp)
  })

  # >
  top.marker.tmp <- diff.marker %>% group_by(cluster)
  top.marker.max <- top.marker.tmp %>%
    filter(avg_log2FC_cap > 0) %>%
    arrange(desc(!!sym(orderBy))) %>%
    slice_head(n = topGeneN)
  top.marker.min <- top.marker.tmp %>%
    filter(avg_log2FC_cap < 0) %>%
    arrange(!!sym(orderBy)) %>%
    slice_head(n = topGeneN)
  top.marker <- rbind(top.marker.max, top.marker.min)

  diff.marker$topN <- ifelse(
    paste0(diff.marker$cluster, diff.marker$gene) %in% paste0(top.marker$cluster, top.marker$gene), diff.marker$gene, "")
  diff.marker$topNShape <- ifelse(
    paste0(diff.marker$cluster, diff.marker$gene) %in% paste0(top.marker$cluster, top.marker$gene), 1, "")
  diff.marker$topNCol <- ifelse(
    paste0(diff.marker$cluster, diff.marker$gene) %in% paste0(top.marker$cluster, top.marker$gene), "grey50", "white")
  diff.marker$topNAlpha <- ifelse(
    paste0(diff.marker$cluster, diff.marker$gene) %in% paste0(top.marker$cluster, top.marker$gene), 1, 0)

  if (!is.null(myMarkers)) {
    top.marker <- diff.marker %>%
      filter(gene %in% myMarkers)
  }else{
    top.marker <- top.marker
  }
  diff.marker.xj <- diff.marker
  diff.marker.xj$xj <- jitter(as.numeric(factor(diff.marker$cluster)), amount = 0.4)

  tile.data <- data.frame(cluster = unique(diff.marker$cluster))

  #> plot
  p1 <- ggplot(data = diff.marker.xj, aes(x = xj, y = avg_log2FC_cap)) +
    geom_col(data = back.data, aes(x = cluster, y = min), fill = "grey98") +
    geom_col(data = back.data, aes(x = cluster, y = max), fill = "grey98") +
    geom_point(aes(color = type)) +
    geom_point( # highlight topN genes
      aes(x = xj, y = avg_log2FC_cap, shape = topNShape),
      size = 2.5, stroke = 1,
      color = diff.marker.xj$topNCol,
      alpha = diff.marker.xj$topNAlpha,
      show.legend = FALSE
    ) +
    geom_tile(
      data = tile.data,
      aes(x = cluster, y = 0, fill = cluster),
      color = "white",
      height = ifelse(log2FC_thres > 1, log2FC_thres*1.8, 1.8),
      alpha = 0.8,
      show.legend = FALSE
    ) +
    geom_text_repel(
      data = diff.marker.xj,
      aes(x = xj, y = avg_log2FC_cap, label = topN),
      color = "grey50",
      max.overlaps = 100
    ) +
    geom_text(
      data = tile.data,
      aes(x = cluster, y = 0, label = str_wrap(cluster, width = 15)),
      color = "grey10",
      size = 3
    ) +
    scale_shape_manual(
      values = c(0, 1)
    ) +
    scale_color_manual(
      values = c(Down = "#a121f0", Up = "#ffad21")
    ) +
    scale_fill_manual(values = col) +
    scale_y_continuous(n.breaks = 10,guide = "prism_offset") +
    # scale_y_continuous(n.breaks = 10) +
    labs(x = "Clusters", y = "Average log2FoldChange",title = plotTitle) +
    guides(
      color = guide_legend(
        override.aes = list(size = 4)
      )
    )+
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5,face = "bold"),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.background = element_blank(),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_text(face = "bold")
    )
  return(p1)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalSampCellComm
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate Niche Cell Communication using CellChat
#'
#' Analyzes cell-cell communication within niches using the CellChat framework.
#' Supports both spatial and non-spatial modes with comprehensive visualization.
#'
#' @param STID_obj An STID object containing niche analysis results
#' @param niche_key Character, niche key to analyze (default: NULL)
#' @param meta_key Character, metadata key for cell annotations (default: NULL)
#' @param group_by Character, column name for cell type grouping
#' @param assay_id Character, assay name (default: "Spatial")
#' @param layer_id Character, layer name (default: "data")
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param col Named vector, colors for cell types
#' @param is_Spatial Logical, whether to use spatial mode (default: TRUE)
#' @param spatial.factors Data frame, spatial scaling factors (default: NULL)
#' @param interaction.range Numeric, interaction range for spatial mode (default: 250)
#' @param remove_genes Character vector, genes to exclude
#' @param return_data Logical, whether to return results list (default: TRUE)
#' @param grp_nm Character, group name for output organization (default: NULL)
#' @param dir_nm Character, directory name for output (default: "M3_CalSampCellComm")
#'
#' @return If return_data = TRUE, returns a list of CellChat objects per sample
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze cell communication in niches
#' results <- CalSampCellComm(
#'   STID_obj = STID_obj,
#'   niche_key = "niche_virulence",
#'   group_by = "cell_type",
#'   is_Spatial = TRUE,
#'   interaction.range = 200
#' )
#' }
CalSampCellComm <- function(STID_obj = NULL,
                            loop_id = "LoopAllSamp",
                            niche_key = NULL,
                            meta_key = NULL,
                            group_by = NULL,
                            assay_id = "Spatial",
                            layer_id = "data",
                            is_Spatial = TRUE, # 是否进行空间cellchat，默认TRUE
                            spatial.factors = NULL,
                            interaction.range = 250,
                            remove_genes = NULL,
                            col = NULL, # COLOR_LIST[["PALETTE_WHITE_BG"]]
                            return_data = TRUE,
                            grp_nm = NULL,dir_nm = "M3_CalSampCellComm"
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  clog_check()
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    clog_error("Package 'CellChat' is required but is not installed.")
  }

  # >>> Check input patameter
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  # .check_null_args(niche_key) # can be NULL
  if(is.null(group_by)){
    clog_warn("group_by is NULL, will use the default celltype_colnm from STID object.")
    group_by <- GetInfo(STID_obj, info_key = "data_info",sub_key = "celltype_colnm")[[1]]
  }
  if(!is.null(col)){
    if(is.null(names(col))){
      clog_error("The col vector must be a named vector, with names corresponding to cell types.")
    }
  }

  #>
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  host_org <- GetInfo(STID_obj, info_key = "samp_info", sub_key = "host_org")[[1]]
  data_format <- GetInfo(STID_obj, info_key = "data_info", sub_key = "data_format")[[1]]
  data_platform <- GetInfo(STID_obj, info_key = "data_info", sub_key = "data_platform")[[1]]
  binsize <- GetInfo(STID_obj, info_key = "data_info", sub_key = "binsize")[[1]]
  coord_interval <- GetInfo(STID_obj, info_key = "data_info", sub_key = "coord_interval")[[1]]
  base_unit <- GetInfo(STID_obj, info_key = "data_info", sub_key = "base_unit")[[1]]

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  clog_step("Execute Single-sample niche cell communication analysis using CellChat...")
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)
  valid_genes <- rownames(STID_obj)[!rownames(STID_obj) %in% remove_genes]

  #> main pipeline
  results_list <- list()
  for(i in seq_along(loop_single)){
    i_single <- loop_single[i]
    clog_loop(paste0("Processing samp_id: ", i_single, " (", i, "/", length(loop_single), ")"))

    #>
    if(!is.null(niche_key)){
      logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
      if(!logic_niche){
        clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
        next
      }
      clog_normal("Using niche_key for CalSampCellComm.")
      Niche_cells <- GetSSNicheCells(STID_obj = STID_obj, niche_key = niche_key, loop_id = i_single)[[1]]
      .check_column_exist(Niche_cells, "is_Niche")
      Niche_cells <- Niche_cells %>%
        filter(is_Niche)
    }else{
      if(is.null(meta_key)){
        clog_warn("Both niche_key and meta_key are NULL, will use meta_key: coord for CalSampCellComm.")
        meta_key <- "coord"
      }
      clog_normal("Using meta_key for CalSampCellComm.")
      Niche_cells <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]] %>%
        filter(!!sym(samp_colnm) == i_single) # actually all, not niche
    }

    #>
    i_coord_interval <- coord_interval[i_single]
    samp_seurat <- subset(STID_obj, cells = rownames(Niche_cells),features = valid_genes) # subset to Niche cells, and remove genes in remove_genes
    clog_normal("Start Cellchat in single sample niche...")
    data.input <- GetAssayData(
      samp_seurat,
      assay = assay_id,
      layer = layer_id,
    )
    meta_data <- GetMetaData(STID_obj = STID_obj,meta_key = "coord")[[1]][rownames(Niche_cells),]
    if(is.factor(meta_data[[group_by]])){
      meta_data[[group_by]] <- meta_data[[group_by]] %>% droplevels()
    }
    len_grp <- length(unique(meta_data[[group_by]]))

    #>
    spatial.locs <- meta_data[c("x","y")]
    if(is_Spatial){
      if(is.null(spatial.factors)){
        if(data_format == "square_grid"){
          spatial.factors <- data.frame(
            ratio = base_unit*binsize/i_coord_interval,
            tol = base_unit*binsize/2
          )
        }else if(data_format == "hex_grid" & data_platform == "Visium"){
          spatial.factors <- data.frame(
            ratio = 100/i_coord_interval,
            tol = base_unit/2
          )
        }
      }
      i_cellchat <- CellChat::createCellChat(
        object = data.input,
        meta = meta_data, # 创建即添加meta
        group.by = group_by,
        datatype = "spatial",
        coordinates = spatial.locs,
        spatial.factors = spatial.factors
      )
    }else{
      i_cellchat <- CellChat::createCellChat(
        object = data.input,
        meta = meta_data,
        group.by = group_by,
        datatype = "RNA"
      )
    }

    #> Set CellChatDB
    if(host_org == "human"){
      CellChatDB <- CellChat::CellChatDB.human
    }else if(host_org == "mouse"){
      CellChatDB <- CellChat::CellChatDB.mouse
    }else{
      clog_error("The host_org must be human or mouse")
    }
    CellChatDB.use <- CellChat::subsetDB(
      CellChatDB,
      search = "Secreted Signaling",
      key = "annotation"
    )
    i_cellchat@DB <- CellChatDB.use

    #> Preprocessing
    cat("\n")
    clog_normal("CellChat preprocessing...")
    i_cellchat <- CellChat::subsetData(i_cellchat)
    i_cellchat <- CellChat::identifyOverExpressedGenes(i_cellchat)
    i_cellchat <- CellChat::identifyOverExpressedInteractions(i_cellchat)

    #> replace this function
    computeCommunProbPathway <- function (object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05) {
      if (is.null(net)) {
        net <- object@net
      }
      if (is.null(pairLR.use)) {
        pairLR.use <- object@LR$LRsig
      }
      prob <- net$prob
      prob[net$pval > thresh] <- 0
      LR <- dimnames(prob)[[3]]
      LR.sig <- LR[apply(prob, 3, sum) != 0]
      pathways <- unique(pairLR.use$pathway_name)
      group <- factor(pairLR.use$pathway_name, levels = pathways)
      prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum), c(2, 3, 1))
      pathways.sig <- pathways[apply(prob.pathways, 3, sum) != 0]
      # >>> 修改
      pathways.sig <- pathways.sig[pathways.sig != ""]
      # <<<
      prob.pathways.sig <- prob.pathways[, , pathways.sig, drop = FALSE]
      idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing = TRUE,
                  index.return = TRUE)$ix
      pathways.sig <- pathways.sig[idx]
      prob.pathways.sig <- prob.pathways.sig[, , idx]
      if (is.null(object)) {
        netP = list(pathways = pathways.sig, prob = prob.pathways.sig)
        return(netP)
      }else {
        object@net$LRs <- LR.sig
        object@netP$pathways <- pathways.sig
        object@netP$prob <- prob.pathways.sig
        return(object)
      }
    }

    #> computeCommunProb
    cat("\n")
    clog_normal("CellChat computeCommunProb, it may take a while ...")
    scale.distance <- signif(1.5/i_coord_interval,digits = 1)
    i_cellchat <- CellChat::computeCommunProb(
      i_cellchat,
      type = "truncatedMean",
      trim = 0.1,
      distance.use = TRUE,
      interaction.range = interaction.range,
      scale.distance = scale.distance,
      contact.range = 100,
      raw.use = TRUE
    )
    results_list[[i_single]][["data"]][["cellchat_raw"]] <- i_cellchat

    #>
    i_cellchat <- CellChat::filterCommunication(
      i_cellchat,
      min.cells = 10
    )
    i_cellchat <- CellChat::computeCommunProbPathway(i_cellchat)
    i_cellchat <- CellChat::aggregateNet(i_cellchat)
    i_cellchat <- CellChat::netAnalysis_computeCentrality(i_cellchat, slot.name = "netP")
    results_list[[i_single]][["data"]][["cellchat_processed"]] <- i_cellchat
    clog_normal("Saving results to RDS file...")
    saveRDS(results_list, paste0(output_dir,"CalSampCellComm_results_list.rds"))

    #>
    celltypes <- i_cellchat@idents %>% levels()
    if(!is.null(col)){
      if(!all(celltypes %in% names(col))){
        clog_error("The names of col vector must include all cell types in the CellChat object.")
      }
      i_col <- col[celltypes]
    }else{
      i_col <- NULL
    }

    #> plot: netVisual_circle
    cat("\n")
    clog_normal("CellChat plotting ...")
    pdf(file = paste0(photo_dir,i_single,"_CellChat_circle_plot.pdf"),
        width = 10,height = 6)
    par(mfrow = c(1,2), xpd=TRUE)
    p1 <- CellChat::netVisual_circle(
      i_cellchat@net$count,
      weight.scale = T,
      label.edge= F,
      title.name = "Number of interactions",
      vertex.label.cex = 1,
      color.use = i_col
    )
    p2 <- CellChat::netVisual_circle(
      i_cellchat@net$weight, # 细胞间互作网络中每个边（互作关系）的权重或强度，与count呈正相关
      weight.scale = T, # 边的宽度会按照权重进行缩放
      label.edge= F, # 不显示边的标签
      title.name = "Weight of interactions",
      vertex.label.cex = 1,
      color.use = i_col
    )
    dev.off()

    #> netVisual_heatmap
    p1 <- CellChat::netVisual_heatmap(i_cellchat,measure = c("count"),color.heatmap = "Reds",color.use = i_col)
    p2 <- CellChat::netVisual_heatmap(i_cellchat,measure = c("weight"),color.heatmap = "Reds",color.use = i_col)
    pdf(file = paste0(photo_dir,i_single,"_CellChat_heatmap_plot.pdf"),
        width = len_grp,height = 1 + len_grp*0.5)
    print(p1+p2)
    dev.off()

    #> netAnalysis_signalingRole_heatmap
    tryCatch({
      p1 <- CellChat::netAnalysis_signalingRole_heatmap(
        i_cellchat,
        pattern = "outgoing", # outgoing/incoming/all
        color.heatmap = "Blues",
        width = len_grp*0.5, height = 10,
        font.size = 5
      )
      p2 <- CellChat::netAnalysis_signalingRole_heatmap(
        i_cellchat,
        pattern = "incoming", # outgoing/incoming/all
        color.heatmap = "Blues",
        width = len_grp*0.5, height = 10,
        font.size = 5
      )
      pdf(file = paste0(photo_dir,i_single,"_CellChat_signalingRole_heatmap_plot.pdf"),
          width = 1 + len_grp*0.5,height =  1 + len_grp*0.5)
      print(p1+p2)
      dev.off()
    },
    error = function(e){
      clog_warn(paste0("Error in netAnalysis_signalingRole_heatmap: ", e$message, " Skipping this plot."))
    })
  }

  # >>> Final
  .save_function_params("CalSampCellComm", envir = environment(), file = paste0(output_dir,"Log_function_params_(CalSampCellComm).log"))
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(CalSampCellComm).log")) %>% invisible()
  if(return_data){
    return(results_list)
  }
}

#' Plot Niche Cell Communication Results
#'
#' Visualizes CellChat analysis results with various plot types including
#' circle plots, heatmaps, bubble plots, and spatial signaling plots.
#'
#' @param STID_obj An STID object containing niche analysis results
#' @param CellComm_data List, CellChat results from CalSampCellComm
#' @param samp_mode  Character, sample type - "SS" or "MS" (default: "SS")
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param sources.use Character vector, source cell types to include
#' @param targets.use Character vector, target cell types to include
#' @param signaling Character vector, signaling pathways to plot
#' @param pairLR.use Data frame, specific ligand-receptor pairs to plot
#' @param col Named vector, colors for cell types
#'
#' @return NULL (invisible), generates plots
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot CellChat results
#' Plot_NicheCellComm(
#'   STID_obj = STID_obj,
#'   CellComm_data = cellcomm_results,
#'   signaling = c("TGFb", "WNT"),
#'   sources.use = "Epithelial",
#'   targets.use = "Immune"
#' )
#' }
Plot_NicheCellComm <- function(STID_obj = NULL,
                               samp_mode  = "SS",
                               loop_id = "LoopAllSamp",
                               CellComm_data = NULL,
                               sources.use = NULL,
                               targets.use = NULL,
                               signaling = NULL,
                               pairLR.use = NULL,
                               col = NULL # COLOR_LIST[["PALETTE_WHITE_BG"]]
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!requireNamespace("nichenetr", quietly = TRUE)) {
    clog_error("Package 'nichenetr' is required but is not installed.")
  }
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(CellComm_data)
  match.arg(samp_mode, .STID_globals$samp_modes)
  if(!is.null(col)){
    if(is.null(names(col))){
      clog_error("The col vector must be a named vector, with names corresponding to cell types.")
    }
  }
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  celltype_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "celltype_colnm")[[1]]

  # >>> Start main pipeline
  clog_step("Plotting single-sample cell communication analysi...")
  if(samp_mode  == "SS"){
    clog_normal("Execute single-sample analysis...")
    loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)
    results_list <- list()
    for(i in seq_along(loop_single)){
      i_single <- loop_single[i]
      clog_loop(paste0("Processing samp_id: ", i_single, " (", i, "/", length(loop_single), ")"))
      if(!i_single %in% names(CellComm_data)){
        clog_warn(paste0("samp_id ", i_single, " not found in CellComm_data, skipping..."))
        next
      }
      i_cellchat <- CellComm_data[[i_single]][["data"]][["cellchat_processed"]]
      datatype <- i_cellchat@options[["datatype"]]
      celltypes <- i_cellchat@idents %>% levels()
      if(!is.null(col)){
        if(!all(celltypes %in% names(col))){
          clog_error("The names of col vector must include all cell types in the CellChat object.")
        }
        i_col <- col[celltypes]
      }else{
        i_col <- NULL
      }

      #> plot
      clog_normal("Plotting circle plot ...")
      par(mfrow = c(1,2), xpd=TRUE) # xpd: 绘制的图形可以超出坐标系的边界
      p1 <- CellChat::netVisual_circle(
        i_cellchat@net$count, # 细胞间互作网络中每个节点（细胞）之间的互作关系数量
        color.use = i_col,
        sources.use = sources.use,
        targets.use = targets.use,
        remove.isolate = T,
        weight.scale = T, # 边的宽度会按照权重进行缩放
        label.edge= F, # 不显示边的标签
        title.name = "Number of interactions",
        vertex.label.cex = 1
      )
      p2 <- CellChat::netVisual_circle(
        i_cellchat@net$weight, # 细胞间互作网络中每个边（互作关系）的权重或强度，与count呈正相关
        color.use = i_col,
        sources.use = sources.use,
        targets.use = targets.use,
        remove.isolate = T,
        weight.scale = T, # 边的宽度会按照权重进行缩放
        label.edge= F, # 不显示边的标签
        title.name = "Weight of interactions",
        vertex.label.cex = 1
      )
      print(p1);print(p2)

      #>
      clog_normal("Plotting heatmap plot ...")
      p1 <- CellChat::netVisual_heatmap(i_cellchat,measure = c("count"),
                                        color.heatmap = "Reds",
                                        remove.isolate = T,
                                        color.use = i_col,
                                        sources.use = sources.use,
                                        targets.use = targets.use)
      p2 <- CellChat::netVisual_heatmap(i_cellchat,measure = c("weight"),
                                        color.heatmap = "Reds",
                                        remove.isolate = T,
                                        color.use = i_col,
                                        sources.use = sources.use,
                                        targets.use = targets.use)
      print(p1+p2)

      #>
      clog_normal("Plotting signalingRole heatmap plot ...")
      p1 <- CellChat::netAnalysis_signalingRole_heatmap(
        i_cellchat,
        signaling = signaling,
        color.heatmap = "Blues",
        pattern = "outgoing", # outgoing/incoming/all
        width = 4, height = 4,
        font.size = 5
      )
      p2 <- CellChat::netAnalysis_signalingRole_heatmap(
        i_cellchat,
        signaling = signaling,
        color.heatmap = "Blues",
        pattern = "incoming", # outgoing/incoming/all
        width = 4, height = 4,
        font.size = 5
      )
      print(p1+p2)

      #>
      if(datatype == "spatial"){
        clog_normal("Plotting spatial signaling plot ...")
        for(j in signaling){
          tryCatch({
            p1 <- CellChat::netVisual_aggregate(i_cellchat,
                                                sources.use = sources.use,
                                                targets.use = targets.use,
                                                signaling = j,
                                                color.use = i_col,
                                                layout = "spatial",
                                                edge.width.max = 2,
                                                vertex.size.max = 1,
                                                alpha.image = 0.4,
                                                vertex.label.cex = 4,
                                                point.size = 1)
            print(p1)
          },
          error = function(e){
            clog_warn(paste0("Error in plotting spatial signaling for pathway: ", j, ". Skipping this pathway."))
          })
        }
      }

      #>
      clog_normal("Plotting bubble plot ...")
      p1 <- CellChat::netVisual_bubble(i_cellchat,
                                       sources.use = sources.use,
                                       targets.use = targets.use,
                                       signaling = signaling,
                                       pairLR.use = pairLR.use, # 默认为所有，配体受体对
                                       # comparison = c(1,2,3,4), # 比较组，会都放在一起
                                       angle.x = 45,
                                       sort.by.target = T,
                                       remove.isolate = T)
      print(p1)
    }
  }else if(samp_mode  == "MS"){
    clog_normal("Execute multi-sample analysis...")
    loop_multi <- .check_loop_multi(STID_obj = STID_obj, loop_id = loop_id)
    for(i in seq_along(loop_multi)){
      i_multi <- loop_multi[i]
      clog_loop(paste0("Processing samp_id: ", i_multi, " (", i, "/", length(loop_multi), ")"))
      samp_id2grp <- data.frame(
        samp_id = STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@samp_info$samp_id,
        samp_grp = STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@samp_info$samp_grp
      )
      if(!all(samp_id2grp$samp_id %in% names(CellComm_data))){
        clog_warn(paste0("Not all samp_id in MultiSampNiche of loop ", i_multi, " are found in CellComm_data, skipping..."))
        next
      }
      compare_mode <- STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@samp_info$compare_mode

      #>
      clog_normal("liftCellChat for multiple samples ...")
      i_data <- CellComm_data[samp_id2grp$samp_id] %>% lapply(function(x){
        x$data[["cellchat_processed"]]
      })
      group.new <- lapply(i_data, function(x){
        x@idents %>% as.character() %>% unique()
      }) %>% unlist() %>% unique()
      i_data <- lapply(i_data, function(x){
        x <- CellChat::liftCellChat(x, group.new)
      })

      #>
      clog_normal("Merging CellChat objects for multiple samples ...")
      i_cellchat <- CellChat::mergeCellChat(i_data, add.names = names(i_data))
      celltypes <- i_cellchat@idents$joint %>% levels()
      if(!is.null(col)){
        if(!all(celltypes %in% names(col))){
          clog_error("The names of col vector must include all cell types in the CellChat object.")
        }
        i_col <- col[celltypes]
      }else{
        i_col <- NULL
      }
      comparison <- seq_along(samp_id2grp$samp_id)
      first_comp <- samp_id2grp$samp_id[1]
      last_comp <- samp_id2grp$samp_id[length(samp_id2grp$samp_id)]

      #>
      clog_normal("Plotting comparative circle plot ...")
      gg1 <- CellChat::compareInteractions(
        i_cellchat, group = comparison, show.legend = F, measure = "count",
        title.name = "Number of interactions")
      gg2 <- CellChat::compareInteractions(
        i_cellchat, group = comparison, show.legend = F, measure = "weight",
        title.name = "Weight of interactions")
      print(gg1 + gg2)

      #> important celltype
      if(compare_mode == "Comparative"){
        clog_normal("Plotting comparative heatmap plot ...")
        gg1 <- CellChat::netVisual_heatmap(
          i_cellchat,comparison = comparison,
          measure = "weight",
          title.name = paste0("Differential interaction strength: ", last_comp," vs ",first_comp),
          sources.use = sources.use, targets.use = targets.use,
          color.use = i_col,
          remove.isolate = T
        )
        print(gg1)
      }

      #> important signaling
      clog_normal("Plotting ranked signaling networks ...")
      gg1 <- CellChat::rankNet(
        i_cellchat,
        comparison = comparison,
        mode = "comparison",
        measure = "weight", # weight有影响
        sources.use = sources.use,
        targets.use = targets.use,
        stacked = F,
        do.stat = TRUE,
        title = paste0("Ranked signaling networks")
      )
      print(gg1)
      gg2 <- CellChat::rankNet(
        i_cellchat,
        comparison = comparison,
        mode = "comparison",
        measure = "weight", # weight有影响
        sources.use = sources.use,
        targets.use = targets.use,
        stacked = T,
        do.stat = TRUE,
        title = paste0("Ranked signaling networks")
      )
      print(gg2)

      #> important ligand-receptor
      clog_normal("Plotting ranked ligand-receptor pairs ...")
      tryCatch({
        gg1 <- CellChat::netVisual_bubble(
          i_cellchat,
          sources.use = sources.use,
          targets.use = targets.use,
          signaling = signaling,
          comparison = comparison,
          max.dataset = 1,
          title.name = paste0("max.dataset: ", first_comp),
          angle.x = 45,
          sort.by.target = T,
          remove.isolate = T)
        print(gg1)
        gg2 <- CellChat::netVisual_bubble(
          i_cellchat,
          sources.use = sources.use,
          targets.use = targets.use,
          signaling = signaling,
          comparison = comparison,
          max.dataset = length(comparison),
          title.name = paste0("max.dataset: ", last_comp),
          angle.x = 45,
          sort.by.target = T,
          remove.isolate = T
        )
        print(gg2)
      },
      error = function(e){
        clog_warn(paste0("Error in plotting ranked ligand-receptor pairs: ", e$message, " Skipping this plot."))
      })
    }

    #> expression of genes in important signaling pathways
    clog_normal("Plotting gene expression of signaling pathways ...")
    i_cellchat@meta$new_group <- paste0(i_cellchat@meta[[samp_colnm]],"_",i_cellchat@meta[[celltype_colnm]])
    tryCatch({
      gg1 <- CellChat::plotGeneExpression(
        i_cellchat,
        signaling = signaling,
        split.by = "new_group",
        group.by = "new_group",
        type = "dot",
        color.use = i_col) +
        labs(title = "Expression of genes in important signaling pathways")
      print(gg1)
    },
    error = function(e){
      clog_warn(paste0("Error in plotting gene expression of signaling pathways: ", e$message))
    })
  }

  #> Final
  clog_end()
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalSampGRN
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate Niche Gene Regulatory Networks using NicheNet
#'
#' Predicts ligand-target gene regulatory networks within niches using the
#' NicheNet framework. Identifies potential signaling from sender to receiver cells.
#'
#' @param STID_obj An STID object containing niche analysis results
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param niche_key Character, niche key to analyze
#' @param group_by Character, column name for cell type grouping
#' @param ref_data List, NicheNet reference data (default: NULL, downloads automatically)
#' @param sender_celltypes Character vector, sender cell types
#' @param receiver_celltypes Character vector, receiver cell types
#' @param target_features Character vector, target genes of interest
#' @param expression_pct Numeric, expression percentage threshold (default: 0.05)
#' @param top_ligand_num Integer, number of top ligands to prioritize (default: 10)
#' @param top_target_num Integer, number of top targets per ligand (default: 5)
#' @param remove_genes Character vector, genes to exclude
#' @param return_data Logical, whether to return results list (default: TRUE)
#' @param grp_nm Character, group name for output organization (default: NULL)
#' @param dir_nm Character, directory name for output (default: "M3_CalSampGRN")
#'
#' @return If return_data = TRUE, returns a list of NicheNet results per sample
#'
#' @importFrom igraph graph_from_data_frame as.undirected V E plot.igraph
#' @import ggraph
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Predict ligand-target networks
#' results <- CalSampGRN(
#'   STID_obj = STID_obj,
#'   niche_key = "niche_virulence",
#'   sender_celltypes = "Epithelial",
#'   receiver_celltypes = "Immune",
#'   target_features = c("gene1", "gene2", "gene3")
#' )
#' }
CalSampGRN <- function(STID_obj = NULL,
                       loop_id = "LoopAllSamp",
                       niche_key = NULL, # no meta_key
                       group_by = NULL,
                       ref_data = NULL,
                       sender_celltypes = NULL,
                       receiver_celltypes = NULL,
                       target_features = NULL,
                       expression_pct = 0.05,
                       top_ligand_num  = 10,
                       top_target_num = 5,
                       remove_genes = NULL,
                       return_data = TRUE, # no col
                       grp_nm = NULL,dir_nm = "M3_CalSampGRN"
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  tmp_dir <- tempdir()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!requireNamespace("nichenetr", quietly = TRUE)) {
    clog_error("Package 'nichenetr' is required but is not installed.")
  }
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(niche_key,sender_celltypes,target_features)
  if(is.null(receiver_celltypes)) {
    celltype_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "celltype_colnm")[[1]]
    receiver_celltypes <- STID_obj@meta.data[[celltype_colnm]]  %>% as.character() %>% unique() %>% sort()
  }
  if(is.null(group_by)){
    clog_warn("group_by is NULL, will use the default celltype_colnm from STID object.")
    group_by <- GetInfo(STID_obj, info_key = "data_info",sub_key = "celltype_colnm")[[1]]
  }
  if(is.null(ref_data)){
    clog_warn("ref_data is NULL,  will download NicheNet reference data, it may take a while...")
    host_org <- GetInfo(STID_obj, info_key = "samp_info", sub_key = "host_org")[[1]]
    if(host_org == "human"){
      lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
      ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
      weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
      gr_network <- readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))
      sig_network <- readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))
      ligand_tf_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final.rds"))
    } else if(host_org == "mouse"){
      lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
      ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
      weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
      gr_network <- readRDS(url("https://zenodo.org/record/7074291/files/gr_network_mouse_21122021.rds"))
      sig_network <- readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_mouse_21122021.rds"))
      ligand_tf_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final_mouse.rds"))
    }
    ref_data <- list(
      weighted_networks = nichenet_weighted_networks,
      lr_network = nichenet_lr_network,
      ligand_target_matrix = nichenet_ligand_target_matrix,
      gr_network = nichenet_gr_network,
      sig_networks = nichenet_signaling_network,
      ligand_tf_matrix = ligand_tf_matrix
    )
  }
  ref_data_nm <- c("weighted_networks","lr_network","ligand_target_matrix","gr_network","sig_networks","ligand_tf_matrix")
  if(!all(ref_data_nm %in% names(ref_data))){
    clog_error(paste0("ref_data must be a list contains: ", paste(ref_data_nm, collapse = ", ")))
  }

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)
  valid_genes <- rownames(STID_obj)[!rownames(STID_obj) %in% remove_genes]

  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  host_org <- GetInfo(STID_obj, info_key = "samp_info", sub_key = "host_org")[[1]]
  ligand_target_matrix <- ref_data$ligand_target_matrix
  lr_network <- ref_data$lr_network
  weighted_networks <- ref_data$weighted_networks
  results_list <- list()
  for(i in seq_along(loop_single)){
    i_single <- loop_single[i]
    clog_loop(paste0("Processing samp_id: ", i_single, " (", i, "/", length(loop_single), ")"))
    logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
    if(!logic_niche){
      clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
      next
    }

    Niche_cells <- GetSSNicheCells(STID_obj = STID_obj, niche_key = niche_key,loop_id = i_single)[[1]]
    .check_column_exist(Niche_cells, "is_Niche")
    Niche_cells <- Niche_cells %>%
      filter(is_Niche)
    clog_normal("Subset niche cells...")
    Idents(STID_obj) <- STID_obj@meta.data[[group_by]]
    samp_seurat <- subset(STID_obj, cells = rownames(Niche_cells),features = valid_genes)
    receiver_celltypes <- receiver_celltypes[receiver_celltypes %in% (Idents(samp_seurat) %>% as.character() %>% unique())]

    clog_normal("Indentify potential ligands and target genes...")
    #> Identify potential ligands from sender cells
    expressed_genes_receiver <- nichenetr::get_expressed_genes(receiver_celltypes, samp_seurat, pct = expression_pct)
    clog_normal(paste0("The number of expressed_genes_receiver is: ",length(expressed_genes_receiver)))
    all_receptors <- unique(lr_network$to)
    expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

    # 用于后续发送者不可知方法
    potential_ligands <- lr_network %>%
      filter(to %in% expressed_receptors) %>% # 根据受体细胞表达的受体，筛选预测的配体
      pull(from) %>% unique()

    # *用于后续以发送者为中心的方法，在不可知基础上，这才是符合我们要求的
    list_expressed_genes_sender <- sender_celltypes %>% unique() %>%
      lapply(nichenetr::get_expressed_genes, samp_seurat, expression_pct)
    expressed_genes_sender <- list_expressed_genes_sender %>%
      unlist() %>% unique()
    potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)
    clog_normal(paste0("The number of potential_ligands_focused is: ",length(potential_ligands_focused)))

    #>
    clog_normal(paste0("The number of target_features is: ",length(target_features)))
    background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    clog_normal(paste0("The number of background_expressed_genes is: ",length(background_expressed_genes)))
    clog_normal("Predict ligand activities, it may take a while...")
    ligand_activities <- nichenetr::predict_ligand_activities(
      geneset = target_features,
      background_expressed_genes = background_expressed_genes,
      ligand_target_matrix = ligand_target_matrix,
      potential_ligands = potential_ligands_focused
    )
    ligand_activities %>% head() # 5列：重要的就是aupr_corrected，其它指标可以考虑
    ligand_activities <- ligand_activities %>%
      arrange(-aupr_corrected) %>%
      mutate(rank = rank(-aupr_corrected))
    dim(ligand_activities)
    results_list[[i_single]][["data"]][["ligand_activities"]] <- ligand_activities

    #> 以发送者为中心的方法
    clog_normal("Visualize ligand activities and top ligands' targets and receptors...")
    ligand_activities <- ligand_activities %>%
      filter(test_ligand %in% potential_ligands_focused) # 根据指定发送者筛选的ligand_activities

    #> 绘制配体活性排名图
    p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) +
      geom_histogram(color="black", fill="darkorange")  +
      geom_vline(aes(xintercept=min(ligand_activities %>%
                                      top_n(30, aupr_corrected) %>%
                                      pull(aupr_corrected))),
                 color="red", linetype="dashed", size=1) +
      labs(x="ligand activity (PCC)", y = "# ligands") +
      theme_classic()
    p_hist_lig_activity
    ggsave(p_hist_lig_activity,
           filename = paste0(photo_dir,"/",i_single,"_Hist_ligand_activity_AUPR.pdf"),
           width = 5, height = 4)
    best_upstream_ligands <- ligand_activities %>%
      top_n(top_ligand_num, aupr_corrected) %>%
      arrange(-aupr_corrected) %>%
      pull(test_ligand) %>% unique()
    len_ligand <- length(best_upstream_ligands)

    # 顶级配体的配体活性测量 (AUPR)
    ligand_aupr_matrix <- ligand_activities %>%
      filter(test_ligand %in% best_upstream_ligands) %>%
      column_to_rownames("test_ligand") %>%
      dplyr::select(aupr_corrected) %>%
      arrange(aupr_corrected)
    head(ligand_aupr_matrix)
    vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1)
    p_ligand_aupr <- nichenetr::make_heatmap_ggplot(vis_ligand_aupr,
                                                    "Prioritized ligands", "Ligand activity",
                                                    legend_title = "AUPR", color = "darkorange") +
      theme(axis.text.x.top = element_blank())
    p_ligand_aupr
    ggsave(p_ligand_aupr,
           filename = paste0(photo_dir,"/",i_single,"_Heatmap_top_ligand_activity_AUPR.pdf"),
           width = 3, height = 2+0.25*len_ligand)

    #> 顶级配体的靶基因
    active_ligand_target_links_df <- best_upstream_ligands %>%
      lapply(nichenetr::get_weighted_ligand_target_links,
             geneset = target_features,
             ligand_target_matrix = ligand_target_matrix,
             n = top_target_num) %>%
      bind_rows() %>% drop_na()
    results_list[[i_single]][["data"]][["active_ligand_target_links_df"]] <- active_ligand_target_links_df
    active_ligand_target_links <- nichenetr::prepare_ligand_target_visualization(
      ligand_target_df = active_ligand_target_links_df,
      ligand_target_matrix = ligand_target_matrix,
      cutoff = 0.05)  # !!!
    dim(active_ligand_target_links)
    order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
    order_targets <- active_ligand_target_links_df$target %>%
      unique() %>%
      intersect(rownames(active_ligand_target_links))
    vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
    p_ligand_target <- nichenetr::make_heatmap_ggplot(
      vis_ligand_target, "Prioritized ligands", "Predicted target genes",
      color = "#1D91C0", legend_title = "Regulatory potential") +
      scale_fill_gradient2(low = "grey95",  high = "#1D91C0")
    p_ligand_target
    ggsave(p_ligand_target,
           filename = paste0(photo_dir,"/",i_single,"_Heatmap_top_ligand_target_regulatory_potential.pdf"),
           width = 2+0.25*ncol(vis_ligand_target), height = 2+0.25*len_ligand)

    #> 顶级配体的受体
    ligand_receptor_links_df <- nichenetr::get_weighted_ligand_receptor_links(
      best_upstream_ligands, expressed_receptors,
      lr_network, weighted_networks$lr_sig)
    results_list[[i_single]][["data"]][["ligand_receptor_links_df"]] <- ligand_receptor_links_df
    vis_ligand_receptor_network <- nichenetr::prepare_ligand_receptor_visualization(
      ligand_receptor_links_df,
      best_upstream_ligands,
      order_hclust = "both")
    p_ligand_receptor <- nichenetr::make_heatmap_ggplot(
      t(vis_ligand_receptor_network),
      y_name = "Ligands", x_name = "Receptors",
      color = "#EF3B2C",
      legend_title = "Prior interaction potential"
    )
    p_ligand_receptor
    ggsave(p_ligand_receptor,
           filename = paste0(photo_dir,"/",i_single,"_Heatmap_top_ligand_receptor_interaction_potential.pdf"),
           width = 2+0.25*(nrow(vis_ligand_receptor_network)), height = 2+0.25*len_ligand)

    #> 发送细胞中的顶级配体表达
    p_dotplot <- DotPlot(subset(samp_seurat, subset = !!sym(group_by) %in% sender_celltypes),
                         features = rev(best_upstream_ligands), cols = "RdYlBu") +
      coord_flip() +
      scale_y_discrete(position = "right")
    p_dotplot
    ggsave(p_dotplot,
           filename = paste0(photo_dir,"/",i_single,"_Dotplot_top_ligand_expression_sender_cells.pdf"),
           width = 3+1*length(sender_celltypes ), height = 2+0.25*len_ligand)

    #>
    clog_normal("Construct signaling network from top ligands to target genes, it may take a while...")
    ligands_all <- best_upstream_ligands # "Ccl3"
    targets_all <- active_ligand_target_links_df$target %>% unique() # "Ccr5"
    active_signaling_network <- nichenetr::get_ligand_signaling_path(
      ligands_all = ligands_all,
      targets_all = targets_all,
      weighted_networks = weighted_networks,
      ligand_tf_matrix = ligand_tf_matrix,
      top_n_regulators = 3,
      minmax_scaling = TRUE
    )

    # >
    tf_signaling <- active_signaling_network$sig
    tf_regulatory <- active_signaling_network$gr
    links <- bind_rows(tf_signaling %>% mutate(edge_type = "sig"), # sig red
                       tf_regulatory %>% mutate(edge_type = "gr")) # gr blue
    nodes <- data.frame(name = unique(c(links$from,links$to))) %>%
      mutate(node_type = if_else(name %in% ligands_all, "Ligands",
                                 if_else(name %in% targets_all, "Targets", "TFs")))
    clog_normal(paste0("The number of nodes is ",nrow(nodes)))
    clog_normal(paste0("The number of links is ",nrow(links)))
    net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = T)
    V(net)$deg <- igraph::degree(net)
    V(net)$size <- igraph::degree(net)/5
    results_list[[i_single]][["data"]][["ligand_target_signaling_paths"]] <- net

    #>
    p1 <- ggraph(net,layout = "kk")+
      geom_edge_fan(aes(edge_width=weight,color = edge_type),
                    arrow = arrow(length = unit(1, 'mm')),
                    end_cap = circle(5, 'mm'),
                    show.legend = F)+
      geom_node_point(aes(fill = node_type,colour = node_type,size = size),
                      # size = 10,
                      color="grey30",
                      stroke = 0.75,
                      shape = 21,
                      alpha=1)+
      geom_node_text(aes(label=name), size = 4, repel = T)+
      scale_edge_width(range = c(0.1,0.5))+
      scale_size_continuous(range = c(5,10) )+
      labs(fill = "Group",title = "ligand-to-target signaling paths")+
      scale_fill_manual(values = c("Ligands" = "#EF3B2C", "Targets" = "#1D91C0","TFs" = "#FFE0B3")) +
      scale_edge_color_manual(values = scales::alpha(c("sig" = "#EF3B2C", "gr" = "#1D91C0"),0.2)) +
      guides(size="none",  fill = guide_legend(override.aes = list(size = 5)))+
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    ggsave(p1,
           filename = paste0(photo_dir,"/",i_single,"_Ligand_to_target_signaling_path.pdf"),
           width = 5+0.05*nrow(nodes), height = 5+0.05*nrow(nodes))
  }
  clog_normal("Saving results to RDS file...")
  saveRDS(results_list, file = paste0(output_dir,"CalSampGRN_results_list.rds"))

  # >>> Final
  .save_function_params("CalSampGRN", envir = environment(), file = paste0(output_dir,"Log_function_params_(CalSampGRN).log"))
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(CalSampGRN).log")) %>% invisible()
  if(return_data){
    return(results_list)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalSampGeneCor: noMS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate Niche Gene Expression Correlations
#'
#' Computes correlation between pathogen and host gene expression within niches.
#' Supports both feature-based and gene-based correlation analysis.
#'
#' @param STID_obj An STID object containing niche analysis results
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param niche_key Character, niche key to analyze
#' @param meta_key Character, metadata key for analysis when niche_key is NULL
#' @param p.features Character vector, pathogen gene names
#' @param h.features Character vector, host gene names
#' @param p.feature_colnm Character vector, pathogen feature column names in metadata
#' @param h.feature_colnm Character vector, host feature column names in metadata
#' @param method Character, correlation method - "spearman" or "pearson"
#'        (default: "spearman")
#' @param assay_id Character, assay name (default: "Spatial")
#' @param layer_id Character, layer name (default: "data")
#' @param return_data Logical, whether to return results list (default: TRUE)
#' @param grp_nm Character, group name for output organization (default: NULL)
#' @param dir_nm Character, directory name for output (default: "M3_CalSampGeneCor")
#'
#' @return If return_data = TRUE, returns a list of correlation results per sample
#'
#' @import Matrix
#' @import ggrepel
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Correlate pathogen and host gene expression
#' results <- CalSampGeneCor(
#'   STID_obj = STID_obj,
#'   niche_key = "niche_virulence",
#'   p.features = c("pathogen_gene1", "pathogen_gene2"),
#'   h.features = c("host_gene1", "host_gene2"),
#'   method = "spearman"
#' )
#' }
CalSampGeneCor <- function(STID_obj = NULL,
                           loop_id = "LoopAllSamp",
                           meta_key = NULL,
                           niche_key = NULL,
                           p.features = NULL,
                           h.features = NULL,
                           p.feature_colnm = NULL,
                           h.feature_colnm = NULL,
                           method = "spearman",
                           assay_id = "Spatial",
                           layer_id = "data",
                           return_data = TRUE,
                           grp_nm = NULL, dir_nm = "M3_CalSampGeneCor"
){
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
  # .check_null_args(niche_key)
  .check_at_least_one_null(p.features,p.feature_colnm)
  .check_at_least_one_null(h.features,h.feature_colnm)
  match.arg(method, choices = c("spearman", "pearson"))
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)
  if(!is.null(p.features)){
    if(!all(p.features %in% rownames(STID_obj))){
      clog_warn(paste0("Some genes in p.features are not found in STID object: ",
                       paste(p.features[!p.features %in% rownames(STID_obj)], collapse = ", ")))
      p.features <- p.features[p.features %in% rownames(STID_obj)]
    }
  }
  if(!all(h.features %in% rownames(STID_obj))){
    clog_warn(paste0("Some genes in h.features are not found in STID object: ",
                     paste(h.features[!h.features %in% rownames(STID_obj)], collapse = ", ")))
    h.features <- h.features[h.features %in% rownames(STID_obj)]
  }
  if(!is.null(p.features)){
    if(length(p.features) == 0){
      clog_warn("No valid genes in p.features.")
    }
  }
  if(length(h.features) == 0){
    clog_error("No valid genes in h.features.")
  }
  clog_normal(paste0(length(p.features), " valid genes in p.features. ",
                     length(h.features), " valid genes in h.features."))
  valid_genes <- unique(c(p.features,h.features))
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  results_list <- list()
  for(i in seq_along(loop_single)){
    i_single <- loop_single[i]
    clog_loop(paste0("Processing samp_id: ", i_single, " (", i, "/", length(loop_single), ")"))

    #>
    if(!is.null(niche_key)){
      logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
      if(!logic_niche){
        clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
        next
      }
      clog_normal("Using niche_key...")
      Niche_cells <- GetSSNicheCells(STID_obj = STID_obj, niche_key = niche_key, loop_id = i_single)[[1]]
      .check_column_exist(Niche_cells, "is_Niche")
      Niche_cells <- Niche_cells %>%
        filter(is_Niche) # only niche, not all
    }else{
      if(is.null(meta_key)){
        clog_warn("Both niche_key and meta_key are NULL, will use meta_key: coord for analysis...")
        meta_key <- "coord"
      }
      clog_normal("Using meta_key...")
      Niche_cells <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]] %>%
        filter(!!sym(samp_colnm) == i_single) # actually all, not niche
    }
    if(!all(p.feature_colnm %in% colnames(Niche_cells))){
      clog_error(paste0("Some columns in p.feature_colnm are not found in Niche_cells: ",
                        paste(p.feature_colnm[!p.feature_colnm %in% colnames(Niche_cells)], collapse = ", ")))
    }else{
      feature_colnm.1_df <- Niche_cells[,p.feature_colnm,drop=F] %>%
        as.matrix() %>% as("dgCMatrix")
    }
    if(!all(h.feature_colnm %in% colnames(Niche_cells))){
      clog_error(paste0("Some columns in h.feature_colnm are not found in Niche_cells: ",
                        paste(h.feature_colnm[!h.feature_colnm %in% colnames(Niche_cells)], collapse = ", ")))
    }else{
      feature_colnm.2_df <- Niche_cells[,h.feature_colnm,drop=F] %>%
        as.matrix() %>% as("dgCMatrix")
    }
    valid_genes_df <- LayerData(STID_obj,assay_id = assay_id,layer = layer_id,
                                cells = rownames(Niche_cells),
                                features = valid_genes) %>%
      t()
    features.1_df <- valid_genes_df[,p.features,drop=F] %>%
      cbind(feature_colnm.1_df, .)
    features.2_df <- valid_genes_df[,h.features,drop=F] %>%
      cbind(feature_colnm.2_df, .)
    clog_normal(paste0("Calculating correlation..."))
    cor_results <- .Cor_pipeline(feat_exp.1 = features.1_df,
                                 feat_exp.2 = features.2_df,
                                 method = method,
                                 fast_cor = T, # must be TRUE
                                 grp_nm = dir_list$grp_nm,
                                 dir_nm = dir_list$dir_nm)
    results_list[[i_single]]$data <- cor_results$result

    #> plot
    plot_data <- cor_results$result
    sign_data_up <- plot_data %>%
      filter(`P-value` < 0.05 & Correlation >0) %>%
      group_by(Variable1) %>%
      arrange(desc(Correlation)) %>%
      slice_head(n = 3)
    sign_data_down <- plot_data %>%
      filter(`P-value` < 0.05 & Correlation <0) %>%
      group_by(Variable1) %>%
      arrange(Correlation) %>%
      slice_head(n = 3)
    sign_data <- rbind(sign_data_up, sign_data_down)
    p1 <- ggplot(plot_data, aes(x=Correlation)) +
      facet_wrap(~ Variable1,ncol = 4,scales = "free") +
      geom_density(fill=NA, color="#99CCFF",linewidth = 1) +
      # geom_histogram(binwidth = 0.05, fill=NA, color="#99CCFF") +
      labs(x="Spearman Correlation", y="Density",
           title=paste0("Niche Gene Correlation Distribution (", i_single,")")) +
      theme_common()
    ggsave(filename = paste0(photo_dir,"/",i_single,"_Correlation_Distribution.pdf"),
           plot = p1, width = 6, height = 5)
    print(p1)

    p2 <- ggplot(plot_data, aes(x=Correlation, y= -log10(`P-value`))) +
      facet_wrap(~ Variable1,ncol = 4,scales = "free") +
      geom_point(alpha = 0.5, size = 0.5,shape = 16,stroke = 0, color="#99CCFF") +
      geom_text_repel(
        data = sign_data,
        aes(x = Correlation, y = -log10(`P-value`), label = Variable2),
        color = "grey50",
        max.overlaps = 100
      ) +
      theme_bw() +
      labs(x="Spearman Correlation", y='-log10(P value)',
           title=paste0("Niche Gene Correlation Scatter Plot (", i_single,")")) +
      theme_common()
    ggsave(filename = paste0(photo_dir,"/",i_single,"_Correlation_Scatter.pdf"),
           plot = p2, width = 6, height = 5)
    print(p2)
  }
  clog_normal("Saving results to RDS file...")
  saveRDS(results_list, file = paste0(output_dir,"/CalSampGeneCor_Results_List.rds"))

  # >>> Final
  .save_function_params("CalSampGeneCor", envir = environment(), file = paste0(output_dir,"Log_function_params_(CalSampGeneCor).log"))
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(CalSampGeneCor).log")) %>% invisible()
  if(return_data){
    return(results_list)
  }
}


#' Correlation Analysis Pipeline
#'
#' Internal function to compute correlations between two feature matrices,
#' supporting both fast parallel and standard modes.
#'
#' @param feat_exp.1 Matrix, first feature matrix
#' @param feat_exp.2 Matrix, second feature matrix
#' @param method Character, correlation method
#' @param fast_cor Logical, whether to use fast parallel computation
#' @param grp_nm Character, group name for output
#' @param dir_nm Character, directory name for output
#'
#' @return List of correlation results
#'
#' @importFrom psych corr.test
#' @importFrom foreach foreach %dopar% %:%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#'
#' @keywords internal
#'
#' @noRd
.Cor_pipeline <- function(feat_exp.1 = NULL,
                          feat_exp.2 = NULL, # 启用feat_exp.2，则进行两个矩阵间的比较
                          method = NULL,
                          fast_cor = TRUE, # 只适用于feat_exp.1与feat_exp.2的比较
                          grp_nm = "stat_exp_cor1",
                          dir_nm = "stat_exp_cor"){
  output_dir <- paste0("./outputdata/",dir_nm,"/",grp_nm) # dir_nm下有子目录，grp_nm一般为gene或其他指定值

  # feat_exp.1 to feat_exp.2
  if(!is.null(feat_exp.1) & !is.null(feat_exp.2)){
    clog_normal("feat_exp.1 cor to feat_exp.2")
    # if((ncol(feat_exp.1) + ncol(feat_exp.2)) <1000){
    #   fast_cor <- F
    # }
    if(fast_cor){
      clog_normal("run fast mode")
      corr_p_value <- function(data1, data2, var1, var2) {
        # tryCatch({ # dense matrix
        #   data_var1 <- data1[[var1]][c(data1[[var1]] != 0 | data2[[var2]] != 0)]
        #   data_var2 <- data2[[var2]][c(data1[[var1]] != 0 | data2[[var2]] != 0)]
        #   corr <- cor(data_var1, data_var2, use = "pairwise.complete.obs",method ="spearman")
        #   res <- cor.test(data_var1, data_var2, use = "pairwise.complete.obs",method ="spearman")
        #   p_value <- res$p.value
        #   len <- min(c(length(na.omit(data_var1)), length(na.omit(data_var2))))
        #   return(c(var1, var2, len, corr, p_value))
        # }, error = function(e) {
        #   return(c(var1, var2, 0, 0, 1))
        # })
        tryCatch({ # dgcMatrix
          x <- data1[, var1]
          y <- data2[, var2]
          nz <- which(x != 0 | y != 0)
          # if (length(nz) == 0) {
          #   return(c(var1, var2, 0, 0, 1))
          # }
          x_clean <- as.numeric(x[nz])
          y_clean <- as.numeric(y[nz])
          corr <- cor(x_clean, y_clean, method = method)
          res <- cor.test(x_clean, y_clean, method = method)
          p_value <- res$p.value
          len <- length(x_clean)
          return(c(var1, var2, len, corr, p_value))
        }, error = function(e) {
          return(c(var1, var2, 0, 0, 1))
        })
      }

      #> parallel
      num_cores <- detectCores() - 1
      num_cores <- if_else(num_cores>12,12,num_cores)
      cl <- makeCluster(num_cores)
      clog_normal(paste0("num_cores:",num_cores))
      registerDoParallel(cl)
      clog_title("Start parallel processing, it may take a while...")
      res_cor <- foreach(var1 = colnames(feat_exp.1), .combine = rbind) %:%
        foreach(var2 = colnames(feat_exp.2), .combine = rbind) %dopar% {
          corr_p_value(feat_exp.1, feat_exp.2, var1, var2)
        }
      stopCluster(cl)

      #> process results
      res_cor <- res_cor %>%
        as.data.frame() %>%
        `colnames<-`(c("Variable1", "Variable2", "length","Correlation", "P-value")) %>%
        mutate(Correlation = as.numeric(Correlation),
               `P-value` = as.numeric(`P-value`)) %>%
        arrange(Variable2, "P-value")
      write.table(res_cor,
                  file = paste0(output_dir,'/cont2cont2_',method,'_fast.txt'),
                  sep = "\t",row.names = F,col.names = T,quote = F)
      clog_normal("End parallel processing...")

      res_list <- list(result = res_cor)
    }else{
      clog_normal("run common mode")
      df <- corr.test(feat_exp.1,feat_exp.2, #
                      use = "pairwise",method = method,
                      adjust = "BH")
      res_r <- df$r
      res_p <- df$p
      res_padj <- df$p.adj
      res_n <- df$n
      write.table(res_r,
                  file = paste0(output_dir,'/cont2cont2_r.txt'),
                  sep = "\t",row.names = T,col.names = NA,quote = F)
      write.table(res_p,
                  file = paste0(output_dir,'/cont2cont2_p.txt'),
                  sep = "\t",row.names = T,col.names = NA,quote = F)
      write.table(res_padj,
                  file = paste0(output_dir,'/cont2cont2_padj.txt'),
                  sep = "\t",row.names = T,col.names = NA,quote = F)
      write.table(res_n,
                  file = paste0(output_dir,'/cont2cont2_n.txt'),
                  sep = "\t",row.names = T,col.names = NA,quote = F)
      res_list <- list(r = res_r,
                       p = res_p,
                       padj = res_padj,
                       n = res_n)
    }
  }
  return(res_list)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CalSampPPI: 没有距离的限制
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Predict Niche Protein-Protein Interactions
#'
#' Predicts pathogen-host protein-protein interactions using BLAST homology
#' search and STRING database integration. Constructs and visualizes PPI networks.
#'
#' @param STID_obj An STID object containing niche analysis results
#' @param p.fasta_path Character, path to pathogen protein FASTA file
#' @param h.fasta_path Character, path to host protein FASTA file
#' @param p.symbol2protein_path Character, pathogen gene symbol to protein ID mapping file
#' @param h.symbol2protein_path Character, host gene symbol to protein ID mapping file
#' @param p.features Character vector, pathogen genes to include
#' @param h.features Character vector, host genes to include
#' @param BLAST_tool Character, BLAST tool to use (default: "blastp")
#' @param dbtype Character, database type (default: "prot")
#' @param BLAST_args Character, BLAST search arguments (default: "-evalue 1e-10 -max_target_seqs 5 -max_hsps 50 -num_threads 4")
#' @param BLAST_fil_args Character, BLAST result filtering criteria
#' @param version_STRING Character, STRING database version (default: "12.0")
#' @param score_thre_only Numeric, score threshold for pathogen-only network (default: 100)
#' @param score_thre_all Numeric, score threshold for full network (default: 700)
#' @param return_data Logical, whether to return results list (default: TRUE)
#' @param grp_nm Character, group name for output organization (default: NULL)
#' @param dir_nm Character, directory name for output (default: "M3_CalSampPPI")
#'
#' @return If return_data = TRUE, returns a list of PPI results
#'
#' @import rBLAST
#' @import rlang
#' @import STRINGdb
#' @importFrom igraph graph_from_data_frame V E degree vertex_attr edge_attr
#'           plot.igraph layout_with_fr
#' @import ggraph
#' @importFrom data.table fread fwrite setnames
#' @importFrom Biostrings readAAStringSet writeXStringSet
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Predict pathogen-host PPIs
#' results <- CalSampPPI(
#'   STID_obj = STID_obj,
#'   p.fasta_path = "pathogen_proteins.fasta",
#'   h.fasta_path = "host_proteins.fasta",
#'   p.symbol2protein_path = "pathogen_mapping.txt",
#'   h.symbol2protein_path = "host_mapping.txt"
#' )
#' }
CalSampPPI <- function(
    STID_obj = NULL,
    p.fasta_path = NULL, h.fasta_path = NULL,
    p.symbol2protein_path = NULL, h.symbol2protein_path = NULL,
    p.features = NULL,
    h.features = NULL, # 相当于限制在了Niche内部
    BLAST_tool = "blastp", dbtype = "prot",
    BLAST_args = "-evalue 1e-10 -max_target_seqs 5 -max_hsps 50 -num_threads 4",
    BLAST_fil_args = "pident > 40, evalue < 1e-10, bitscore > 100, length > 100",
    version_STRING = "12.0",
    score_thre_only = 100,score_thre_all = 700,
    return_data = TRUE,
    grp_nm = NULL,dir_nm = "M3_CalSampPPI"
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  tmp_dir <- tempdir()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(p.fasta_path, h.fasta_path, p.symbol2protein_path, h.symbol2protein_path)
  if(is.null(p.features)){
    p.features <- GetInfo(STID_obj, info_key = "data_info", sub_key = "pathogen_genes")[[1]]
  }
  p.features <- p.features[p.features %in% rownames(STID_obj)]
  if(is.null(h.features)){
    h.features <- rownames(STID_obj)
  }
  h.features <- h.features[h.features %in% rownames(STID_obj)]
  clog_normal(paste0(length(p.features), " valid genes in p.features. ",
                     length(h.features), " valid genes in h.features."))
  # >>> End check

  # >>> Check input
  if(has_blast()){
    clog_normal("You have successfully configured the BLAST tool.")
    clog_normal(paste0("The ",BLAST_tool," path is ",Sys.which(BLAST_tool)))
    clog_normal(paste0("The makeblastdb path is ",Sys.which("makeblastdb")))
  }else{
    clog_error("
      Can't detect the BLAST tools:
      1.Please install the BLAST tool from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/.
        For example, download ncbi-blast-2.16.0+-win64.exe.
      2.set the PATH using Sys.setenv(PATH = paste(Sys.getenv(\"PATH\"), <Your BLAST tool dir path>, sep=.Platform$path.sep)).
      3.check if the BLAST tool is in your PATH using Sys.which(\"blastn\").
      ")
  }
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start pipeline
  clog_step(paste0("Execute Niche PPI prediction using BLAST..."))
  results_list <- list()
  results_list$data <- list()
  clog_title(paste0("read host AAStringSet"))
  h.fasta_seq <- readAAStringSet(filepath = h.fasta_path)
  headers <- names(h.fasta_seq)
  h.protein <- strsplit(headers,split = " ") %>%
    sapply(., function(x){x[1]})
  h.symbol2protein <- read.table(file = h.symbol2protein_path,
                                 sep = "\t",header = T,stringsAsFactors = F)
  if(!all(c("symbol","protein") %in% colnames(h.symbol2protein))){
    clog_error("The h.symbol2protein must contain the columns: symbol, protein")
  }
  h.symbol2protein$symbol <- gsub("_","-",h.symbol2protein$symbol)
  h.valid_protein <- h.symbol2protein$protein[h.symbol2protein$symbol %in% h.features] %>% unique()
  h.fasta_seq_fil <- h.fasta_seq[h.protein %in% h.valid_protein]
  # writeXStringSet(h.fasta_seq_fil, filepath = paste0(tmp_dir,"/host_fil.fasta"))

  # >
  clog_title(paste0("read pathogen AAStringSet"))
  p.fasta_seq <- readAAStringSet(filepath = p.fasta_path)
  headers <- names(p.fasta_seq)
  p.protein <- strsplit(headers,split = " ") %>%
    sapply(., function(x){x[1]})
  p.symbol2protein <- read.table(file = p.symbol2protein_path,
                                 sep = "\t",header = T,stringsAsFactors = F)
  if(!all(c("symbol","protein") %in% colnames(p.symbol2protein))){
    clog_error("The p.symbol2protein must contain the columns: symbol, protein")
  }
  p.symbol2protein$symbol <- gsub("_","-",p.symbol2protein$symbol)
  p.valid_protein <- p.symbol2protein$protein[p.symbol2protein$symbol %in% p.features] %>% unique()
  p.fasta_seq_fil <- p.fasta_seq[p.protein %in% p.valid_protein]
  # gene <- strsplit(headers,split = " ") %>%
  #   sapply(., function(x){x[SYMBOL_col]}) %>%
  #   gsub("\\[locus_tag=","",.) %>%
  #   gsub("\\]","",.) %>%
  #   gsub("\\[gene=","",.)
  # ID <- strsplit(headers,split = " ") %>%
  #   sapply(., function(x){x[1]})
  # p.symbol2protein <- data.frame(ID = ID, gene = gene)

  #>
  clog_title(paste0("make host blastdb"))
  host_org <- GetInfo(STID_obj, info_key = "samp_info", sub_key = "host_org")[[1]]
  pathogen_org <- GetInfo(STID_obj, info_key = "samp_info", sub_key = "pathogen_org")[[1]]
  db_dir <- paste0(tmp_dir,"/db/",host_org)
  dir.create(db_dir,recursive = T, showWarnings = F)
  if(!is.null(h.fasta_path)){
    makeblastdb(
      # file = paste0(tmp_dir,"/host_fil.fasta"), # 不应该比对输入的宿主基因
      file = h.fasta_path,
      db_name = paste0(db_dir,"/",host_org),
      dbtype = dbtype #  ("nucl" or "prot")
    )
  }
  db <- blast(paste0(db_dir,"/",host_org),type = BLAST_tool)

  #>
  clog_step(paste0("Predict pathogen-host PPIs using BLAST"))
  clog_warn("Running a BLAST search with 5000 queries against 100000 sequences takes approximately 2 minutes using 16 threads.")
  clog_normal("Start BLAST prediction, it may take a while...")
  res_blast_predict <- predict(db, p.fasta_seq_fil,BLAST_args = BLAST_args)
  clog_normal("End BLAST prediction...")
  if(nrow(res_blast_predict)==0){
    clog_error("The res_blast_predict has zero rows, please modify the BLAST_args!")
  }else{
    clog_normal(paste0("The nrow of res_blast_predict is ",nrow(res_blast_predict)))
  }
  res_blast_predict <- res_blast_predict %>%
    mutate(query_gene = p.symbol2protein$gene[match(qseqid,p.symbol2protein$ID)],.before = 2) %>%
    mutate(ref_gene = h.symbol2protein$gene[match(sseqid,h.symbol2protein$ID)],.before = 4)
  write.table(res_blast_predict, file = paste0(output_dir,"/res_blast_predict.txt"),
              sep = "\t", quote = F,row.names = F,col.names = T)
  results_list$data[["res_blast_predict"]] <- res_blast_predict

  #>
  # BLAST_fil_args <- "pident > 60, evalue < 1e-10, bitscore > 100, length > 100"
  clog_title("Filter BLAST results using BLAST_fil_args")
  if(!is.null(BLAST_fil_args)){
    conditions <- strsplit(BLAST_fil_args, ",")[[1]]
    combined_condition <- paste(conditions, collapse = " & ")
    res_blast_predict_fil <- res_blast_predict %>%
      filter(eval(parse(text = combined_condition)))
  }else{
    res_blast_predict_fil <- res_blast_predict
  }
  res_blast_predict_fil <- res_blast_predict_fil %>%
    arrange(evalue) %>%
    dplyr::select(qseqid,sseqid) %>%
    unique() %>%
    mutate(p.symmbol = p.symbol2protein$symbol[match(qseqid,p.symbol2protein$protein)],
           p.symbol_mg = paste0("p.",p.symmbol),
           h.symbol = h.symbol2protein$symbol[match(sseqid,h.symbol2protein$protein)])

  if(nrow(res_blast_predict_fil)==0){
    clog_error("The res_blast_predict_fil has zero rows, please modify the BLAST_fil_args!")
  }else{
    clog_normal(paste0("The nrow of res_blast_predict_fil is ",nrow(res_blast_predict_fil)))
    print(head(res_blast_predict_fil))
  }
  write.table(res_blast_predict_fil, file = paste0(output_dir,"/res_blast_predict_fil.txt"),
              sep = "\t", quote = F,row.names = F,col.names = T)
  results_list$data[["res_blast_predict_fil"]] <- res_blast_predict_fil

  #
  clog_step("Construct PPI network using STRINGdb")
  if(host_org == "human"){
    species_id <- 9606
  }else if(host_org == "mouse"){
    species_id <- 10090
  }else{
    stop("The host_org must be human or mouse")
  }
  STRING_dir <- paste0(tmp_dir,"/db/STRING/")
  string_db <- STRINGdb$new( # 需要联网
    version=version_STRING,
    species=species_id,
    score_threshold = score_thre_only, # 可信度？
    network_type = "full",
    input_directory=STRING_dir # 这里必须要下面的三个文件，否则后面的许多分析报错
  )

  #
  dir.create(STRING_dir,recursive = T, showWarnings = F)
  link_aliases <- paste0("https://stringdb-downloads.org/download/protein.aliases.v",version_STRING,"/",species_id,".protein.aliases.v",version_STRING,".txt.gz")
  link_info <- paste0("https://stringdb-downloads.org/download/protein.info.v",version_STRING,"/",species_id,".protein.info.v",version_STRING,".txt.gz")
  link_links <- paste0("https://stringdb-downloads.org/download/protein.links.v",version_STRING,"/",species_id,".protein.links.v",version_STRING,".txt.gz")
  if(all(file.exists(paste0(STRING_dir,basename(link_aliases)),
                     paste0(STRING_dir,basename(link_info)),
                     paste0(STRING_dir,basename(link_links))))){
    clog_normal("STRINGdb files already exist, skip download.")
  }else{
    clog_normal(paste0("Download STRINGdb files to ",STRING_dir))
    download.file(link_aliases,destfile = paste0(STRING_dir,basename(link_aliases)),method = "wininet")
    download.file(link_info,destfile = paste0(STRING_dir,basename(link_info)),method = "wininet")
    download.file(link_links,destfile = paste0(STRING_dir,basename(link_links)),method = "wininet")
  }
  prot_id_map <- res_blast_predict_fil %>%
    mutate(sseqid2 = gsub("\\..*","",sseqid)) %>%
    string_db$map(  # 支持多种ID，自动转换HUGO，Entrez GeneID, ENSEMBL
      my_data_frame_id_col_names = "sseqid",  # gene symbol一般对应Ensembl_UniProt
      removeUnmappedRows = TRUE
    )

  #>>>> only
  cat("\n")
  clog_step(paste0("Create PPI network by only using pathogen proteins"))
  hit <- prot_id_map$STRING_id %>% unique()
  clog_normal(paste0("The nrow of hit pathogen proteins is ",length(hit)))
  if(length(hit) >20){
    clog_warn("The nrow of hit pathogen proteins is greater than 20, it may take a while...")
  }
  # from convert to pathogen genes, and to convert to host genes
  links <- string_db$get_interactions(hit) %>%  # !!! 不需要网络，但是需要下载的三个文件
    unique() # 其实是双向的，两个方向的权重略有不同。
  if(0){
    pdf(file = paste0(photo_dir,"PPI_STRING_network.pdf"),
        width = 6,height = 6)
    string_db$plot_network(hit)
    dev.off()
  }
  if(nrow(links)==0){
    clog_warn("No interactions found in STRINGdb for the given pathogen proteins.")
  }else{
    clog_normal(paste0("The nrow of links is ",nrow(links)))
    links <- links %>%
      mutate(from = prot_id_map$p.symbol_mg[match(from, prot_id_map$STRING_id)]) %>%
      mutate(to = prot_id_map$h.symbol[match(to, prot_id_map$STRING_id)]) %>%
      dplyr::select(from,to,combined_score) %>%
      transform(from = pmin(from, to), to = pmax(from, to)) %>%
      filter(to %in% h.features) %>% # !!!!
      unique() %>%
      # 由于from和to都可能是pathogen，所以拷贝一份反向的
      bind_rows(.,data.frame(from = .$to,to = .$from,combined_score = .$combined_score))
    nodes <- data.frame(name = unique(c(links$from,links$to))) %>%
      mutate(is_Mg = ifelse(grepl("p\\.",name),"pathogen","host")) %>%
      mutate(Mg_nm = ifelse(grepl("p\\.",name),name,NA))
    clog_normal(paste0("The number of nodes is ",nrow(nodes)))
    clog_normal(paste0("The number of links is ",nrow(links)))
    net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
    V(net)$deg <- igraph::degree(net)
    V(net)$size <- igraph::degree(net)/5

    # ggraph
    nodes <- igraph::as_data_frame(net, what = "vertices") # !!!!!! 重要，节点转变成df查看属性
    links <- igraph::as_data_frame(net, what = "edges") # !!!!!! 重要，边转变成df查看属性
    write.table(links,file = paste0(output_dir,"/links_only_pathogen.txt"),sep = "\t",
                row.names = F,col.names = T,quote = F)
    write.table(nodes,file = paste0(output_dir,"/nodes_only_pathogen.txt"),sep = "\t",
                row.names = F,col.names = T,quote = F)
    results_list$data[["PPI_network_only_pathogen"]] <- list(nodes = nodes, links = links)

    #>
    p1 <- ggraph(net,layout = "kk")+
      geom_edge_fan(aes(edge_width=combined_score), color = "grey95", show.legend = F)+
      geom_node_point(aes(size=size,color = is_Mg,fill = is_Mg),
                      color="grey95",
                      stroke = 0.15,
                      shape = 21,
                      alpha=0.8)+
      geom_node_text(aes(label=name), size = 2, repel = T)+
      scale_edge_width(range = c(0.1,0.5))+
      scale_size_continuous(range = c(3,8) )+
      labs(fill = "Group",title = paste0("combined_score > ",score_thre_only)) +
      scale_fill_manual(values = c("pathogen" = "#a121f0", "host" = "#ffad21")) +
      guides(size="none",  fill = guide_legend(override.aes = list(size = 5)))+
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
    pdf(file = paste0(photo_dir,"PPI_network_only_pathogen.pdf"),
        width = 6,height = 6)
    print(p1)
    dev.off()
  }

  #>>>> all
  # only其实也可以用这第二种方法，直接读取，非string_db$get_interactions
  clog_step(paste0("Create PPI network using all pathogen-host proteins"))
  clog_title("Read STRINGdb links file, it may take a while...")
  ID_gene_df <- fread(paste0(STRING_dir,basename(link_aliases)),
                      header = TRUE, sep = "\t", fill=TRUE) %>%
    as.data.frame() %>%
    filter(source == "UniProt_GN_Name")
  links <- fread(paste0(STRING_dir,basename(link_links)), header = TRUE, sep = " ") %>% # 数据是单向的，应该使用双向的
    as.data.frame() %>%
    `colnames<-`(c("from","to","combined_score")) %>%
    filter(combined_score>score_thre_all) %>%
    # transform(from = pmin(from, to), to = pmax(from, to)) %>%  # mutate会修改from了，但是transform不会
    filter(from %in% prot_id_map$STRING_id | to %in% prot_id_map$STRING_id ) %>%  # 两边要有病原的
    mutate(from = ID_gene_df$alias[match(from, ID_gene_df$`#string_protein_id`)]) %>%
    mutate(to = ID_gene_df$alias[match(to, ID_gene_df$`#string_protein_id`)]) %>%
    na.omit() %>%
    filter(from %in% h.features | to %in% h.features) %>%  # 两边必须要有宿主的
    filter(from %in% prot_id_map$h.symbol | to %in% prot_id_map$h.symbol ) %>% # 再次过滤两边也要有病原的
    arrange(from,to) %>%
    unique()
  links_all_pathogen <- links %>% # 和only一样，这种无法判断左右
    filter(from %in% prot_id_map$h.symbol & to %in% prot_id_map$h.symbol ) %>%
    mutate(from = prot_id_map$p.symbol_mg[match(from,prot_id_map$h.symbol)]) %>%
    unique() %>%
    bind_rows(.,data.frame(from = .$to,to = .$from,combined_score = .$combined_score))
  links_only_one <- links %>%
    filter(!(from %in% prot_id_map$h.symbol & to %in% prot_id_map$h.symbol) ) %>%
    mutate(from = ifelse(from %in% prot_id_map$h.symbol,
                         prot_id_map$p.symbol_mg[match(from,prot_id_map$h.symbol)],from)) %>%
    mutate(to = ifelse(to %in% prot_id_map$h.symbol,
                       prot_id_map$p.symbol_mg[match(to,prot_id_map$h.symbol)],to)) %>%
    unique()
  links_fil <- rbind(links_all_pathogen,links_only_one)
  nodes <- data.frame(name = unique(c(links_fil$from,links_fil$to))) %>%
    mutate(is_Mg = ifelse(grepl("p\\.",name),"pathogen","host")) %>%
    mutate(Mg_nm = ifelse(grepl("p\\.",name),name,NA))
  clog_normal(paste0("The number of nodes is ",nrow(nodes)))
  clog_normal(paste0("The number of links is ",nrow(links_fil)))
  net <- igraph::graph_from_data_frame(d=links_fil,vertices=nodes,directed = F)
  V(net)$deg <- igraph::degree(net)
  V(net)$size <- igraph::degree(net)/5
  nodes <- igraph::as_data_frame(net, what = "vertices") # !!!!!! 重要，节点转变成df查看属性
  links <- igraph::as_data_frame(net, what = "edges") # !!!!!! 重要，边转变成df查看属性
  write.table(nodes,file = paste0(output_dir,"/nodes_all.txt"),sep = "\t",
              row.names = F,col.names = T,quote = F)
  write.table(links_fil,file = paste0(output_dir,"/links_all.txt"),sep = "\t",
              row.names = F,col.names = T,quote = F)
  results_list$data[["PPI_network_all"]] <- list(nodes = nodes, links = links)

  #>
  clog_title("Visualize PPI network using ggraph")
  p1 <- ggraph(net,layout = "kk")+
    geom_edge_fan(aes(edge_width=combined_score), color = "grey95", show.legend = F)+
    geom_node_point(aes(size=size,color = is_Mg,fill = is_Mg),
                    color="grey95",
                    stroke = 0.15,
                    shape = 21,
                    alpha=0.8)+
    geom_node_text(aes(label=Mg_nm), size = 4, repel = T)+
    scale_edge_width(range = c(0.1,0.5))+
    scale_size_continuous(range = c(3,8) )+
    labs(fill = "Group",title = paste0("combined_score > ",score_thre_all)) +
    scale_fill_manual(values = c("pathogen" = "#a121f0", "host" = "#ffad21")) +
    guides(size="none",  fill = guide_legend(override.aes = list(size = 5)))+
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  pdf(file = paste0(photo_dir,"PPI_network_all.pdf"),
      width = 8,height = 7)
  print(p1)
  dev.off()

  p2 <- ggraph(net,layout = "kk")+
    geom_edge_fan(aes(edge_width=combined_score), color = "grey95", show.legend = F)+
    geom_node_point(aes(size=size,color = is_Mg,fill = is_Mg),
                    color="grey95",
                    stroke = 0.15,
                    shape = 21,
                    alpha=0.8)+
    geom_node_text(aes(label=name), size = 4, repel = T)+
    scale_edge_width(range = c(0.1,0.5))+
    scale_size_continuous(range = c(3,8) )+
    labs(fill = "Group",title = paste0("combined_score > ",score_thre_all)) +
    scale_fill_manual(values = c("pathogen" = "#a121f0", "host" = "#ffad21")) +
    guides(size="none",  fill = guide_legend(override.aes = list(size = 5)))+
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  pdf(file = paste0(photo_dir,"PPI_network_all_show_gene.pdf"),
      width = 8,height = 7)
  print(p2)
  dev.off()

  #> save results
  clog_normal("Saving results to RDS file...")
  saveRDS(results_list, file = paste0(output_dir,"CalSampPPI_results_list.rds"))

  # >>> Final
  .save_function_params("CalSampPPI", envir = environment(), file = paste0(output_dir,"Log_function_params_(CalSampPPI).log"))
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(CalSampPPI).log")) %>% invisible()
  if(return_data){
    return(results_list)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot_DistLine
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Plot Gene Expression Along Distance from Niche Center
#'
#' Visualizes gene expression patterns as a function of distance from niche
#' centers, with optional smoothing and scaling.
#'
#' @param STID_obj An STID object containing niche analysis results
#' @param features Character vector, gene names to plot
#' @param feature_colnm Character vector, metadata column names to plot
#' @param facet_grpnm Character, column name for faceting
#' @param meta_key Character, metadata key containing distance information
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param smooth_method Character, smoothing method for geom_smooth (default: "gam")
#' @param exp_scale Logical, whether to scale expression values (default: TRUE)
#' @param distance_scale Logical, whether to scale distance values (default: TRUE)
#' @param col Color palette (default: COLOR_LIST$PALETTE_WHITE_BG)
#' @param linewidth Numeric, line width (default: 1)
#' @param ncol Integer, number of facet columns (default: 4)
#' @param assay_id Character, assay name (default: "Spatial")
#' @param layer_id Character, layer name (default: "data")
#'
#' @return A ggplot object
#'
#' @importFrom ggh4x facet_wrap2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot gene expression along distance from niche center
#' Plot_DistLine_Exp(
#'   STID_obj = STID_obj,
#'   features = c("gene1", "gene2"),
#'   meta_key = "M2_NicheDetect_STS_20240101"
#' )
#' }
Plot_DistLine_Exp <- function(STID_obj = NULL,
                              loop_id = "LoopAllSamp",
                              meta_key = NULL,
                              facet_grpnm = NULL,
                              assay_id = "Spatial", layer_id = "data",
                              features = NULL, feature_colnm = NULL,
                              smooth_method = "gam",
                              exp_scale = TRUE,distance_scale = TRUE,
                              linewidth = 1,ncol = 4,
                              col = COLOR_LIST[["PALETTE_WHITE_BG"]]

){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(col)
  .check_at_least_one_null(features,feature_colnm,facet_grpnm,meta_key)
  if(length(STID_obj@STID_analysis@SingleSampNiche) == 0){
    loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id, mode = 1)
  }else{
    loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)
  }

  # >
  meta_data <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]]
  .check_column_exist(meta_data,facet_grpnm)
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  if(!all(feature_colnm %in% colnames(meta_data))){
    missing_cols <- feature_colnm[!feature_colnm %in% colnames(meta_data)]
    clog_error(paste0("The following feature_colnm are not found in meta_data: ", paste(missing_cols, collapse = ", ")))
  }else{
    meta_data[feature_colnm] <- lapply(meta_data[feature_colnm], as.numeric)
    meta_data <- meta_data[c(samp_colnm,"x","y","ROI_edge","All_Dist2ROIcenter",feature_colnm)]
  }
  if(!is.null(features)){
    valid_features <- .check_features_exist(STID_obj, features)
    feat_exp <- FetchData(STID_obj, vars = valid_features, assay = assay_id, layer = layer_id)
    meta_data <- bind_cols(meta_data[rownames(feat_exp),],feat_exp,.name_repair = "unique")
  }else{
    valid_features <- features
  }
  if(!is.null(feature_colnm)){
    if(all(feature_colnm %in% colnames(meta_data))){
      valid_features <- c(valid_features, feature_colnm)
    }
  }
  len_feat <- length(valid_features)

  # >>> Start main pipeline
  plot_data <- meta_data %>%
    na.omit() %>%
    filter(!!sym(samp_colnm) %in% loop_single) %>%
    group_by(!!sym(samp_colnm)) %>%
    {
      if(distance_scale){
        mutate(.,All_Dist2ROIcenter = All_Dist2ROIcenter/quantile(All_Dist2ROIcenter,0.95)) %>% # !!!!!!!
          mutate(.,All_Dist2ROIcenter = if_else(All_Dist2ROIcenter>1,1,All_Dist2ROIcenter))
      }else{
        mutate(.,All_Dist2ROIcenter = if_else(All_Dist2ROIcenter>quantile(All_Dist2ROIcenter,0.95),
                                              quantile(All_Dist2ROIcenter,0.95),All_Dist2ROIcenter))
      }
    } %>%
    mutate(All_Dist2ROIedge = median(All_Dist2ROIcenter[ROI_edge]),.after = "All_Dist2ROIcenter") %>%  # mean is better than median
    pivot_longer(cols = c(-!!sym(samp_colnm),-All_Dist2ROIcenter,-All_Dist2ROIedge,-ROI_edge,-x,-y),
                 names_to = "gene",values_to = "exp") %>%
    mutate(gene = factor(gene,levels = valid_features)) %>%
    mutate(exp = as.numeric(exp))
  if(exp_scale){
    plot_data <- plot_data %>%
      group_by(!!sym(samp_colnm),gene) %>%
      mutate(exp = as.numeric(scale(exp)))
  }
  clog_normal(paste0("Plotting gene expression along distance, this may take a while..."))
  rect_data <- plot_data %>%
    group_by(gene, !!sym(samp_colnm)) %>%
    summarise(All_Dist2ROIedge = unique(All_Dist2ROIedge),.groups = "drop") %>%
    mutate(
      xmin = -Inf,
      xmax = All_Dist2ROIedge,
      ymin = -Inf,
      ymax = Inf
    )
  p1 <- ggplot(plot_data,aes(x = All_Dist2ROIcenter,y = exp,color = gene))+
    facet_grid2(cols = vars(!!sym(samp_colnm)),rows = vars(gene),
                scales = "free", independent = "y")+
    # facet_wrap(~gene,ncol = ncol,scales  = "free") +
    geom_vline(
      data = rect_data,
      aes(xintercept = All_Dist2ROIedge),
      color = "#EB1E2C",
      lty = "dashed",
      lwd = 0.5
    ) +
    geom_rect(
      data = rect_data,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = scales::alpha("#E64B35FF", 0.1),
      color = NA, inherit.aes = FALSE
    ) +
    # annotate("rect",xmin = -Inf, xmax = unique(plot_data$All_Dist2ROIedge),ymin = -Inf ,ymax = Inf,fill = scales::alpha("#E64B35FF",0.1),color = NA)+
    geom_smooth(method = smooth_method,linewidth = linewidth,se = T,level = 0.95,fill = "grey90") +
    labs(y = "Average expression", x = paste0("Spot order along distance (Edge:",round(unique(plot_data$All_Dist2ROIedge),2),")"),
         color = "Features") +
    scale_color_manual(values = col) +
    theme_common()
  print(p1)

  clog_end()
}

#' Plot Cell Type Ratio Along Distance from Niche Center
#'
#' Visualizes the fraction of cell types as a function of distance from niche
#' centers, with distance binning based on spatial coord_intervals.
#'
#' @param STID_obj An STID object containing niche analysis results
#' @param celltypes Character vector, cell types to plot
#' @param group_by Character, column name for cell type grouping
#' @param facet_grpnm Character, column name for faceting
#' @param meta_key Character, metadata key containing distance information
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param coord_interval_ratio Numeric, multiplier for coord_interval to set bin width
#' @param col Color palette (default: COLOR_LIST$PALETTE_WHITE_BG)
#' @param linewidth Numeric, line width (default: 1)
#' @param ncol Integer, number of facet columns (default: 4)
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot cell type ratios along distance from niche center
#' Plot_DistLine_Ratio(
#'   STID_obj = STID_obj,
#'   celltypes = c("Tcell", "Bcell", "Macrophage"),
#'   group_by = "cell_type",
#'   meta_key = "M2_NicheDetect_STS_20240101",
#'   coord_interval_ratio = 5
#' )
#' }
Plot_DistLine_Ratio <- function(STID_obj = NULL,
                                loop_id = "LoopAllSamp",
                                meta_key = NULL,
                                group_by = NULL,
                                facet_grpnm = NULL,
                                celltypes = NULL,
                                coord_interval_ratio = NULL ,
                                linewidth = 1,ncol = 4,
                                col = COLOR_LIST[["PALETTE_WHITE_BG"]]

){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(col)
  .check_at_least_one_null(celltypes,group_by,facet_grpnm,meta_key,coord_interval_ratio)
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)

  # >
  meta_data <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]]
  .check_column_exist(meta_data, group_by,facet_grpnm)
  all_celltypes <- meta_data[[group_by]] %>% na.omit() %>% as.character() %>% unique()
  if(!all(celltypes %in% all_celltypes)){
    invalid_celltype <- setdiff(celltypes, all_celltypes)
    clog_error(paste0("Invalid celltypes: ", paste0(invalid_celltype,collapse = ",")))
  }

  # >>> Start main pipeline
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  coord_interval <- GetInfo(STID_obj, info_key = "data_info",sub_key = "coord_interval")[[1]]
  plot_data <-   meta_data  %>%
    filter(!!sym(samp_colnm) %in% loop_single) %>%
    mutate(i_coord_interval = coord_interval[!!sym(samp_colnm)]) %>%
    group_by(!!sym(samp_colnm)) %>%
    mutate(cut_value = cut_width(All_Dist2ROIcenter, width = unique(i_coord_interval)*coord_interval_ratio, boundary = 0, closed = "right"),
           cut_value = as.numeric(cut_value),
           cut_value = if_else(cut_value > 10, 10, cut_value),
           cut_label = paste0("L",cut_value),
           cut_label = factor(cut_label, levels = paste0("L",1:10))) %>%
    ungroup()
  clog_normal(paste0("The Results of spot counts in each distance bin:"))
  print(
    table(plot_data[c(samp_colnm,"cut_label")])
  )
  plot_data_line <- plot_data %>%
    group_by(!!sym(samp_colnm),cut_label) %>%
    mutate(cellnum_cut = n()) %>%
    group_by(!!sym(samp_colnm),cut_label,!!sym(group_by)) %>%
    summarise(celltype_ratio = n()/dplyr::first(cellnum_cut)) %>%
    as.data.frame() %>%
    complete(!!sym(samp_colnm),cut_label,!!sym(group_by),fill = list(celltype_ratio = 0))
  if(!is.null(celltypes)){
    plot_data_line <- plot_data_line %>% filter(!!sym(group_by) %in% celltypes)
  }
  if(nrow(plot_data_line) == 0){
    clog_error("No data to plot, please check the input parameters.")
  }
  plot_data_edge <- plot_data %>%
    filter(ROI_edge) %>%
    group_by(!!sym(samp_colnm)) %>%
    summarise(median_edge_value = median(cut_value)) %>%
    mutate(median_edge = paste0("L",median_edge_value))
  col_hex_add <- .hex_add(col,add_s = -0.2,add_v = 0.2)
  names(col_hex_add) <- names(col)
  p1 <- ggplot(data = plot_data_line,aes(x = cut_label, y = celltype_ratio,
                                         color = !!sym(group_by),fill = !!sym(group_by),
                                         group = !!sym(group_by)))+
    facet_grid(cols = vars(!!sym(samp_colnm)),rows = vars(!!sym(group_by)),scales = "free")+
    geom_rect(data = plot_data_edge,
              aes(x= NULL,y = NULL,xmin = -Inf, xmax = median_edge,ymin = -Inf ,ymax = Inf,),
              color = NA,fill = scales::alpha("#E64B35FF",0.1),inherit.aes = F)+
    geom_vline(data = plot_data_edge,aes(xintercept = median_edge),
               alpha = 0.8,color = scales::alpha("#E64B35FF",0.1),
               lty = "dashed", lwd = 0.5)+
    geom_line(linewidth=1,alpha = 0.8) +
    geom_point(size=2.5,shape=21,stroke = 1) +
    scale_color_manual(values = scales::alpha(col,0.9))+
    scale_fill_manual(values = scales::alpha(col_hex_add,0.9))+
    labs(x=paste0("Distance to ROI center (",coord_interval_ratio,"* coord_interval)"),
         y="Fraction of group",color = "Groups",fill = "Groups")+
    # scale_x_continuous(
    #   breaks = scales::breaks_width(1),  # coord_interval of 1
    #   labels = function(x) as.integer(x)
    # ) +
    scale_y_continuous(
      labels = scales::percent
    ) +
    theme_common()
  print(p1)

  clog_end()
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GeneEnrichment
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Perform Gene Enrichment Analysis (GO/KEGG or GSEA)
#'
#' Performs gene enrichment analysis using either traditional over-representation
#' analysis (GO/KEGG) or Gene Set Enrichment Analysis (GSEA). Supports both
#' human and mouse organisms.
#'
#' @param up_gene Character vector, upregulated genes for GO/KEGG analysis
#' @param down_gene Character vector, downregulated genes for GO/KEGG analysis
#' @param DEGs Data frame, differential expression results with columns 'SYMBOL'
#'        and 'LogFC' for GSEA analysis
#' @param STID_obj STID object (optional, used to extract host organism)
#' @param host_org Character, host organism - "human" or "mouse" (if NULL,
#'        extracted from STID_obj)
#' @param method Character, method to use - "GO_KEGG" or "GSEA_GO_KEGG"
#'        (default: "GO_KEGG")
#' @param plot_pathway_num Integer, number of top pathways to display in plots
#'        (default: 12)
#' @param col Color palette for plots (default: \code{viridis(100, option = "H")[15:85] \%>\% rev()})
#' @param internal Logical, whether to use internal KEGG database (default: FALSE)
#' @param return_data Logical, whether to return enrichment results (default: FALSE)
#' @param grp_nm Character, group name for output organization (default: NULL)
#' @param dir_nm Character, directory name for output (default: "M3_GeneEnrichment")
#' @param ... Additional arguments passed to internal functions
#'
#' @return If return_data = TRUE, returns list of enrichment results; otherwise NULL
#'
#' @importFrom viridis viridis
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # GO/KEGG enrichment with up and down regulated genes
#' GeneEnrichment(
#'   up_gene = upregulated_genes,
#'   down_gene = downregulated_genes,
#'   host_org = "human",
#'   plot_pathway_num = 10
#' )
#'
#' # GSEA enrichment with DEGs
#' GeneEnrichment(
#'   DEGs = deg_results,
#'   host_org = "mouse",
#'   method = "GSEA_GO_KEGG"
#' )
#' }
GeneEnrichment <- function(STID_obj = NULL,
                           up_gene = NULL,
                           down_gene = NULL,
                           DEGs = NULL,
                           host_org = NULL, # "human"/"mouse"
                           method = "GO_KEGG", # GO_KEGG, GSEA_GO_KEGG
                           plot_pathway_num = 12,
                           internal = FALSE,
                           col = viridis(100, option = "H")[15:85] %>% rev(),
                           return_data = FALSE, # !!!
                           grp_nm = NULL, dir_nm = "M3_GeneEnrichment",
                           ...) {
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check()
  if(is.null(host_org)){
    if (!inherits(STID_obj, "STID")) {
      clog_error("host_org is NULL and input object is not an STID object")
    }
    host_org <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "host_org") %>% unlist(use.names = F)
  }
  match.arg(method, choices = c("GO_KEGG", "GSEA_GO_KEGG"))
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  if(method == "GO_KEGG"){
    clog_step("Performing GO and KEGG enrichment analysis")
    if(is.null(up_gene) & is.null(down_gene)){
      clog_error("Please provide up or down gene lists for enrichment analysis.")
    }
    results_list <- .GO_KEGG(up = up_gene,
                             down = down_gene,
                             org = host_org,
                             grp_nm = dir_list$grp_nm,
                             dir_nm = dir_list$dir_nm,
                             col = col,
                             internal = internal,
                             GO_num = plot_pathway_num,
    )
  }else if(method == "GSEA_GO_KEGG"){
    clog_step("Performing GSEA GO and KEGG enrichment analysis")
    .check_null_args(DEGs)
    results_list <- .GSEA_GO_KEGG(DEGs = DEGs,
                                  org = host_org,
                                  grp_nm = dir_list$grp_nm,
                                  dir_nm = dir_list$dir_nm,
                                  col = col,
                                  internal = internal,
                                  GO_num = plot_pathway_num,
    )
  }
  clog_normal("Saving results to RDS file...")
  saveRDS(results_list, file = paste0(output_dir,"GeneEnrichment_results_list.rds"))

  # >>> Final
  .save_function_params("GeneEnrichment", envir = environment(), file = paste0(output_dir,"Log_function_params_(GeneEnrichment).log") )
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(GeneEnrichment).log")) %>% invisible()
  if (return_data) {
    return(results_list)
  }
}


#' Perform GO and KEGG Over-representation Analysis
#'
#' Internal function for traditional GO and KEGG enrichment analysis using
#' clusterProfiler. Performs gene ID conversion, enrichment testing, and
#' generates bar plots and bubble plots.
#'
#' @param updown Character vector, combined up and down regulated genes
#' @param up Character vector, upregulated genes
#' @param down Character vector, downregulated genes
#' @param simplify_index Logical, whether to simplify redundant GO terms
#' @param org Character, organism - "human" or "mouse" (default: "human")
#' @param go_ont Character vector, GO ontology types - "BP", "CC", "MF"
#'        (default: c("BP","CC","MF"))
#' @param internal Logical, whether to use internal KEGG database (default: TRUE)
#' @param GO_num Integer, number of top GO terms to display (default: 20)
#' @param col Color palette for plots (default: c('#ED0000FF','#00468BFF'))
#' @param maxGSSize Integer, maximum gene set size (default: 500)
#' @param grp_nm Character, group name for output
#' @param dir_nm Character, directory name for output
#'
#' @return List of enrichment results for each gene set
#'
#' @importFrom magrittr %>% add
#' @import clusterProfiler
#' @import DOSE
#' @import patchwork
#'
#' @keywords internal
#'
#' @noRd
.GO_KEGG <- function(updown = NULL, up = NULL, down = NULL,
                     simplify_index = F,org = "human", # human/mouse
                     go_ont = c("BP","CC","MF"), # c("BP","CC","MF")
                     internal = T, # T/F
                     GO_num = 20,
                     col = c('#ED0000FF','#00468BFF'), # 颜色
                     maxGSSize = 500,
                     grp_nm = "GO_KEGG",dir_nm = "GO_KEGG"){

  clog_normal(paste0("clusterProfiler version: ", packageVersion("clusterProfiler")))
  clog_normal(paste0("Your organism is: ", org))
  if(org == "human"){
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      clog_error("Package 'org.Hs.eg.db' is required for human annotation but is not installed.")
    }
    clog_normal(paste0("org.Hs.eg.db version: ", packageVersion("org.Hs.eg.db")))
    org_index <- org.Hs.eg.db
    kegg_org <- "hsa"
  }else if(org == "mouse"){
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      clog_error("Package 'org.Mm.eg.db' is required for mouse annotation but is not installed.")
    }
    clog_normal(paste0("org.Mm.eg.db version: ", packageVersion("org.Mm.eg.db")))
    org_index <- org.Mm.eg.db
    kegg_org <- "mmu"
  }
  output_dir <- paste0("./outputdata/",dir_nm,"/",grp_nm,"/")
  photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm,"/")
  patchwork_index <- F

  #> ID conversion
  tryCatch({
    ent_updown <- bitr(updown, # 基因名大写
                       fromType = 'SYMBOL',
                       toType = c('ENSEMBL','ENTREZID'), # 'ENTREZID'最适合富集
                       OrgDb = org_index)$ENTREZID
  },error = function(e){
    ent_updown <<- NULL
  })
  tryCatch({
    ent_up <- bitr(up,
                   fromType = 'SYMBOL',
                   toType = c('ENSEMBL','ENTREZID'),
                   OrgDb = org_index)$ENTREZID
  },error = function(e){
    ent_up <<- NULL
  })
  tryCatch({
    ent_down <- bitr(down,
                     fromType = 'SYMBOL',
                     toType = c('ENSEMBL','ENTREZID'),
                     OrgDb = org_index)$ENTREZID
  },error = function(e){
    ent_down <<- NULL
  })

  #> GO and KEGG
  vector_ent <- c()
  if(!is_empty(ent_updown)){
    vector_ent <- c(vector_ent,as.character(quote(ent_updown))) # "ent_updown"
  }
  if(!is_empty(ent_up)){
    vector_ent <- c(vector_ent,as.character(quote(ent_up)))
  }
  if(!is_empty(ent_down)){
    vector_ent <- c(vector_ent,as.character(quote(ent_down)))
  }
  # go_ont <- c("BP","CC","MF")
  results_list <- list()
  for (go_i in vector_ent) { # "ent_updown", "ent_up", "ent_down"
    results_list[[go_i]] <- list()

    tryCatch({
      ## 3.1 go
      for (go_j in go_ont) { # "BP","CC","MF
        clog_normal(paste0("Start GO: ",go_i," ",go_j))
        tryCatch({
          go <- enrichGO(gene = get(go_i),  # entrezgene_id
                         OrgDb = org_index , # 物种注释包，所有GO本身就是不联网的，取决于物种注释包的版本
                         keyType = "ENTREZID", # ID类型
                         ont = go_j, # ALL,MF:molecular function ,BP:biological process ,CC:cellular compotent
                         pAdjustMethod = "BH", # 多重假设检验矫正方法
                         pvalueCutoff  = ifelse(simplify_index,0.2,1), # 设置过滤标准，影响go的矩阵数据，影响后续许多分析
                         qvalueCutoff  = 1, # 同上
                         maxGSSize = maxGSSize,
                         pool = T, # ont=all,同时输出3个通路
                         readable = T ) # 将go结果中的ENTREZID转为gene symbol
          # assign(paste0('go_',go_i,'_',go_j),go)

          ## simplify
          if(simplify_index){
            print("Start simplify")
            go <- clusterProfiler::simplify(go,cutoff=0.7,by="p.adjust",select_fun=min,measure="Wang")
            print("End simplify")
          }

          ## go_res
          # go <- setReadable(go,OrgDb=org_index, keyType = "ENTREZID") # enrichGO中已经有readable = T
          go_res <- go@result # go@result是原始的所有结果，go保存的是按照指标过滤后的结果
          if(nrow(go_res)==0){
            patchwork_index=T
            next
          }
          go_res$Des_short <-  sapply(go_res$Description,.shorten_names)
          go_res$group <- grp_nm
          results_list[[go_i]][[go_j]] <- go_res
          # assign(paste0('go_',go_i,'_',go_j,'_res'),go_res)
          write.table(go_res,
                      file = paste0(output_dir,'/GO_',go_i,'_',go_j,'_res.txt'),
                      row.names =F,col.names = T,quote = F,sep = '\t')

          ## barplot
          if(nrow(go_res)>=GO_num){
            go_data <- go_res[1:GO_num,] # 自动挑选感兴趣的通路
          }else if(nrow(go_res)>=0){
            go_data <- go_res[1:nrow(go_res),]
          }else{
            # patchwork_index=T
            # next
          }
          go_data <- go_data %>%
            arrange(Count,-p.adjust) %>%
            mutate(Des_short_unique=paste0(Des_short,1:nrow(go_data))) %>%
            mutate(Des_short_unique=factor(Des_short_unique,levels=unique(Des_short_unique)))

          pdf(file = paste0(photo_dir,'/GO_',go_i,'_',go_j,'_bar.pdf'),
              width = 6.5, height = 4.5)
          p1 <- ggplot(go_data,aes(x=Count,y=Des_short_unique,fill=p.adjust))+
            geom_col(position="identity")+
            geom_text(aes(label=Count),hjust=-0.3)+
            theme_test()+
            labs(x='Counts',y='',title = paste0("GO of ",go_i,'_',go_j," of ",grp_nm))+
            # scale_fill_gradient(low = col[1],high = col[2])+
            scale_fill_gradientn(colors = col) +
            scale_x_continuous(expand = c(0,0,0.2,0))+
            scale_y_discrete(expand = c(0,1),labels=go_data$Des_short)+
            theme(plot.margin = margin(0.4,0.4,0.4,0.4,'cm'),
                  axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
                  axis.text  = element_text(size = 10,face = 'bold',hjust = 0.5),
                  legend.title = element_text(size = 11,face = 'bold',hjust = 0),
                  plot.title = element_text(size = 14,face = 'bold',hjust = 0.5))
          print(p1)
          assign(paste0('GO_',go_i,'_',go_j,'_bar'),p1)
          dev.off()

          ## bubble diagram
          GeneRatio_all <- go_res$GeneRatio[1] %>%
            gsub(".*/","",.) %>%
            as.numeric()
          go_data <- go_data %>%
            mutate(GeneRatio_all=GeneRatio_all,.before = 10)

          pdf(file = paste0(photo_dir,'/GO_',go_i,'_',go_j,'_bubble.pdf'),
              width = 6.5, height = 4.5)
          p2 <- ggplot(go_data,aes(x=Count/GeneRatio_all,y=Des_short_unique,size=Count,color=p.adjust))+
            geom_point()+
            theme_bw()+
            scale_size(range = c(3, 8))+
            scale_color_gradientn(colors = col) +
            scale_x_continuous(expand = c(0.1,0,0.1,0))+
            scale_y_discrete(expand = c(0,1),labels=go_data$Des_short)+
            labs(x='Gene Ratio',y='',color='P.adjust',title = paste0("GO of ",go_i,'_',go_j," of ",grp_nm))+
            guides(color=guide_colorbar(order=0),size=guide_legend(order=1))+
            theme(plot.margin = margin(0.4,0.4,0.4,0.4,'cm'),
                  axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
                  axis.text  = element_text(size = 10,face = 'bold',hjust = 0.5),
                  legend.title = element_text(size = 11,face = 'bold',hjust = 0),
                  plot.title = element_text(size = 14,face = 'bold',hjust = 0.5))
          print(p2)
          assign(paste0('GO_',go_i,'_',go_j,'_bubble'),p2)
          dev.off()

          rm(go,go_res,go_data,GeneRatio_all)
        },
        error=function(e) print(e))
      }
    },
    error=function(e) print(e))

    ## 3.2 KEGG
    tryCatch({
      clog_normal(paste0("Start KEGG: ",go_i))
      kegg <- enrichKEGG(get(go_i),
                         organism = kegg_org,
                         pvalueCutoff = ifelse(simplify_index,0.2,1),
                         keyType = "kegg",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 1,
                         maxGSSize = maxGSSize,
                         use_internal_data = internal)
      # assign(paste0('kegg_',go_i),kegg)

      ## kegg_res
      kegg <- setReadable(kegg, # 自动ent→sym
                          OrgDb = org_index,
                          keyType="ENTREZID")
      kegg_res <- kegg@result
      kegg_res$Des_short <-  sapply(kegg_res$Description,.shorten_names)
      kegg_res$group <- grp_nm
      results_list[[go_i]][["KEGG"]] <- kegg_res
      # assign(paste0('kegg_',go_i,'_res'),kegg_res)
      write.table(kegg_res,
                  file = paste0(output_dir,'/KEGG_',go_i,'_res.txt'),
                  row.names =F,col.names = T,quote = F,sep = '\t')

      ## barplot
      if(nrow(kegg_res)>=GO_num){
        kegg_data <- kegg_res[1:GO_num,] #挑选感兴趣的通路
      }else if(nrow(kegg_res)>=0){
        kegg_data <- kegg_res[1:nrow(kegg_res),]
      }else{
        # patchwork_index=T
        # next
      }
      kegg_data <- kegg_data %>%
        arrange(Count,-p.adjust) %>%
        mutate(Des_short_unique=paste0(Des_short,1:nrow(kegg_data))) %>%
        mutate(Des_short_unique=factor(Des_short_unique,levels=unique(Des_short_unique)))

      pdf(file = paste0(photo_dir,'/KEGG_',go_i,'_bar.pdf'),
          width = 6.5, height = 4.5)
      p3 <- ggplot(kegg_data,aes(x=Count,y=Des_short_unique,fill=p.adjust))+
        geom_col(position="identity")+
        geom_text(aes(label=Count),hjust=-0.3)+
        theme_test()+
        labs(x='Counts',y='',title = paste0("KEGG of ",go_i," of ",grp_nm))+
        # scale_fill_gradient(low = col[1],high = col[2])+
        scale_fill_gradientn(colors = col) +
        scale_x_continuous(expand = c(0,0,0.2,0))+
        scale_y_discrete(expand = c(0,1),labels=kegg_data$Des_short)+
        theme(plot.margin = margin(0.4,0.4,0.4,0.4,'cm'),
              axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
              axis.text  = element_text(size = 10,face = 'bold',hjust = 0.5),
              legend.title = element_text(size = 11,face = 'bold',hjust = 0),
              plot.title = element_text(size = 14,face = 'bold',hjust = 0.5))
      print(p3)
      assign(paste0('KEGG_',go_i,'_bar'),p3)
      dev.off()

      ## bubble diagram
      GeneRatio_all <- kegg_res$GeneRatio[1] %>%
        gsub(".*/","",.) %>%
        as.numeric()
      kegg_data <- kegg_data %>%
        mutate(GeneRatio_all=GeneRatio_all,.before = 10)

      pdf(file = paste0(photo_dir,'/KEGG_',go_i,'_bubble.pdf'),
          width = 6.5, height = 4.5)
      p4 <- ggplot(kegg_data,aes(x=Count/GeneRatio_all,y=Des_short_unique,size=Count,color=p.adjust))+
        geom_point()+
        theme_bw()+
        scale_size(range = c(3, 8))+
        scale_color_gradientn(colors = col) +
        scale_x_continuous(expand = c(0.1,0,0.1,0))+
        scale_y_discrete(expand = c(0,1),labels=kegg_data$Des_short)+
        labs(x='Gene Ratio',y='',color='P.adjust',title = paste0("KEGG of ",go_i," of ",grp_nm))+
        guides(color=guide_colorbar(order=0),size=guide_legend(order=1))+
        theme(plot.margin = margin(0.4,0.4,0.4,0.4,'cm'),
              axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
              axis.text  = element_text(size = 10,face = 'bold',hjust = 0.5),
              legend.title = element_text(size = 11,face = 'bold',hjust = 0),
              plot.title = element_text(size = 14,face = 'bold',hjust = 0.5))
      print(p4)
      assign(paste0('KEGG_',go_i,'_bubble'),p4)
      dev.off()

      rm(kegg,kegg_res,kegg_data)
    },
    error=function(e) print(e))
  }

  ### 4. combine pictures
  if(patchwork_index){
    return(NULL) %>% invisible()
  }

  ## up
  variables <- c("GO_ent_up_BP_bar", "GO_ent_up_MF_bar", "GO_ent_up_CC_bar",
                 "KEGG_ent_up_bar", "GO_ent_up_BP_bubble", "GO_ent_up_MF_bubble",
                 "GO_ent_up_CC_bubble", "KEGG_ent_up_bubble")
  missing_vars <- variables[!sapply(variables, function(x) exists(x))]
  if(length(missing_vars) == 0){
    tryCatch({
      if(!is_empty(ent_up)){
        pdf(file = paste0(photo_dir,'/GO_KEGG_up.pdf'),
            width = 26,height = 8,family = "serif")
        print(
          GO_ent_up_BP_bar+GO_ent_up_MF_bar+GO_ent_up_CC_bar+KEGG_ent_up_bar+
            GO_ent_up_BP_bubble+GO_ent_up_MF_bubble+GO_ent_up_CC_bubble+KEGG_ent_up_bubble+
            plot_layout(nrow = 2,guides = "auto")
        )
        dev.off()
      }
    },
    error=function(e) {
      print(e)
      dev.off()
    })
  }

  ## down
  variables <- c("GO_ent_down_BP_bar", "GO_ent_down_MF_bar", "GO_ent_down_CC_bar",
                 "KEGG_ent_down_bar", "GO_ent_down_BP_bubble", "GO_ent_down_MF_bubble",
                 "GO_ent_down_CC_bubble", "KEGG_ent_down_bubble")
  missing_vars <- variables[!sapply(variables, function(x) exists(x))]
  if(length(missing_vars) == 0){
    tryCatch({
      if(!is_empty(ent_down)){
        pdf(file = paste0(photo_dir,'/GO_KEGG_down.pdf'),
            width = 26,height = 8,family = "serif")
        print(
          GO_ent_down_BP_bar+GO_ent_down_MF_bar+GO_ent_down_CC_bar+KEGG_ent_down_bar+
            GO_ent_down_BP_bubble+GO_ent_down_MF_bubble+GO_ent_down_CC_bubble+KEGG_ent_down_bubble+
            plot_layout(nrow = 2,guides = "auto")
        )
        dev.off()
      }
    },
    error=function(e) {
      print(e)
      dev.off()
    })
  }


  ## updpwn
  variables <- c("GO_ent_updown_BP_bar", "GO_ent_updown_MF_bar", "GO_ent_updown_CC_bar",
                 "KEGG_ent_updown_bar", "GO_ent_updown_BP_bubble", "GO_ent_updown_MF_bubble",
                 "GO_ent_updown_CC_bubble", "KEGG_ent_updown_bubble")
  missing_vars <- variables[!sapply(variables, function(x) exists(x))]
  if(length(missing_vars) == 0){
    tryCatch({
      if(!is_empty(ent_updown)){
        pdf(file = paste0(photo_dir,'/GO_KEGG_updown.pdf'),
            width = 26,height = 8,family = "serif")
        print(
          GO_ent_updown_BP_bar+GO_ent_updown_MF_bar+GO_ent_updown_CC_bar+KEGG_ent_updown_bar+
            GO_ent_updown_BP_bubble+GO_ent_updown_MF_bubble+GO_ent_updown_CC_bubble+KEGG_ent_updown_bubble+
            plot_layout(nrow = 2,guides = "auto")
        )
        dev.off()
      }
    },
    error=function(e) {
      print(e)
      dev.off()
    })
  }
  return(results_list)
}

#' Perform GSEA GO and KEGG Analysis
#'
#' Internal function for Gene Set Enrichment Analysis using clusterProfiler.
#' Calculates enrichment scores for GO terms and KEGG pathways based on
#' ranked gene lists.
#'
#' @param DEGs Data frame, differential expression results with columns
#'        'gene' and 'avg_log2FC'
#' @param org Character, organism - "human" or "mouse" (default: "human")
#' @param simplify_index Logical, whether to simplify redundant GO terms
#' @param maxGSSize Integer, maximum gene set size (default: 500)
#' @param col Color palette for plots (default: c('#ED0000FF','#00468BFF'))
#' @param internal Logical, whether to use internal KEGG database
#' @param GO_num Integer, number of top pathways to display (default: 20)
#' @param grp_nm Character, group name for output
#' @param dir_nm Character, directory name for output
#'
#' @return List of GSEA results for GO and KEGG
#'
#' @import clusterProfiler
#' @import DOSE
#' @import patchwork
#'
#' @keywords internal
#'
#' @noRd
.GSEA_GO_KEGG <- function(DEGs = NULL,
                          org = "human", # human/mouse
                          simplify_index = F,
                          maxGSSize = 500,
                          col = c('#ED0000FF','#00468BFF'),
                          internal = T,
                          GO_num = 20,
                          grp_nm = "GSEA_GO_KEGG",
                          dir_nm = "GSEA_GO_KEGG") {


  #> check
  .check_null_args(DEGs)
  must_colnm <- c("gene","avg_log2FC")
  if(!all(c(must_colnm) %in% colnames(DEGs))){
    clog_error(paste0("DEGs must contain columns 'gene' and 'avg_log2FC', but the following columns are missing: ",
                      paste(setdiff(c(must_colnm), colnames(DEGs)), collapse = ", ")))
  }

  #>
  clog_normal(paste0("clusterProfiler version: ", packageVersion("clusterProfiler")))
  clog_normal(paste0("Your organism is: ", org))

  if(org == "human"){
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      clog_error("Package 'org.Hs.eg.db' is required for human annotation but is not installed.")
    }
    clog_normal(paste0("org.Hs.eg.db version: ", packageVersion("org.Hs.eg.db")))
    org_index <- org.Hs.eg.db
    kegg_org <- "hsa"
  }else if(org == "mouse"){
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      clog_error("Package 'org.Mm.eg.db' is required for mouse annotation but is not installed.")
    }
    clog_normal(paste0("org.Mm.eg.db version: ", packageVersion("org.Mm.eg.db")))
    org_index <- org.Mm.eg.db
    kegg_org <- "mmu"
  }else{
    clog_error("org must be either 'human' or 'mouse'")
  }
  output_dir <- paste0("./outputdata/", dir_nm, "/", grp_nm, "/")
  photo_dir <- paste0("./photo/", dir_nm, "/", grp_nm, "/")


  #>
  gsea_df <- DEGs[, c("gene", "avg_log2FC")]
  colnames(gsea_df) <- c('SYMBOL', 'LogFC')
  gsea_df <- gsea_df[order(gsea_df$LogFC, decreasing = T), ]

  #> ID转换
  gsea_SYMBOL <- gsea_df$SYMBOL
  gsea_ENTREZID <- bitr(gsea_SYMBOL,
                        fromType = 'SYMBOL',
                        toType = c('ENSEMBL', 'ENTREZID'),
                        OrgDb = org_index)
  gsea_df <- merge(gsea_df, gsea_ENTREZID, by = "SYMBOL")

  #> 命名LogFC向量
  gsea_LogFC <- gsea_df$LogFC
  names(gsea_LogFC) <- as.character(gsea_df$ENTREZID)
  gsea_LogFC <- gsea_LogFC[!is.na(names(gsea_LogFC)) | names(gsea_LogFC) == ""]
  gsea_LogFC[is.infinite(gsea_LogFC) & gsea_LogFC > 0] <- max(gsea_LogFC[!is.infinite(gsea_LogFC)]) + 1
  gsea_LogFC[is.infinite(gsea_LogFC) & gsea_LogFC < 0] <- min(gsea_LogFC[!is.infinite(gsea_LogFC)]) - 1
  gsea_LogFC <- sort(gsea_LogFC, decreasing = T)

  #> GSEA GO分析 (BP, CC, MF)
  go_ont <- c("BP", "CC", "MF")
  results_list <- list()

  for (go_i in go_ont) {
    clog_normal(paste0("Start GSEA GO: ", go_i))

    browser()

    tryCatch({
      gsea_GO <- gseGO(
        gsea_LogFC,
        ont = go_i,
        OrgDb = org_index,
        keyType = "ENTREZID",
        exponent = 1,
        minGSSize = 10,
        maxGSSize = maxGSSize,
        eps = 1e-16,
        pvalueCutoff = ifelse(simplify_index, 0.2, 1),
        pAdjustMethod = "BH",
        verbose = TRUE,
        seed = TRUE,
        by = "fgsea"
      )

      #> Simplify
      if(simplify_index && !is.null(gsea_GO) && nrow(gsea_GO@result) > 0){
        clog_normal(paste0("Start simplify for GO: ", go_i))
        gsea_GO <- clusterProfiler::simplify(gsea_GO, cutoff = 0.7, by = "p.adjust",
                                             select_fun = min, measure = "Wang")
      }

      #> 结果处理
      gsea_GO <- setReadable(gsea_GO, OrgDb = org_index, keyType = "ENTREZID")
      gsea_GO_res <- gsea_GO@result

      if(nrow(gsea_GO_res) == 0) {
        clog_warning(paste0("No significant results for GO: ", go_i))
        next
      }

      gsea_GO_res$Des_short <- sapply(gsea_GO_res$Description, .shorten_names)
      gsea_GO_res$group <- grp_nm
      results_list[[paste0("GO_", go_i)]] <- gsea_GO_res

      #> 保存结果
      write.table(gsea_GO_res,
                  file = paste0(output_dir, 'GSEA_GO_', go_i, '_res.txt'),
                  row.names = F, col.names = T, quote = F, sep = '\t')

      #> 绘制上调通路图 (NES > 0)
      go_data_up <- gsea_GO_res %>%
        filter(NES > 0) %>%
        mutate(new_NES = round(NES, digits = 2))

      if(nrow(go_data_up) > 0) {
        go_data <- go_data_up
        if(nrow(go_data) >= GO_num) {
          go_data <- go_data[1:GO_num, ]
        }

        go_data <- go_data %>%
          arrange(new_NES, -p.adjust) %>%
          mutate(Des_short_unique = paste0(Des_short, 1:nrow(go_data))) %>%
          mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique)))

        pdf(file = paste0(photo_dir, 'GSEA_GO_', go_i, '_up_bar.pdf'),
            width = 6.5, height = 4.5)
        p_up <- ggplot(go_data, aes(x = new_NES, y = Des_short_unique, fill = p.adjust)) +
          geom_col(position = "identity") +
          geom_text(aes(label = new_NES), hjust = -0.15) +
          theme_test() +
          labs(x = 'NES', y = '', title = paste0("GSEA GO (Up) of ", go_i, " in ", grp_nm)) +
          scale_fill_gradientn(colors = col) +
          scale_x_continuous(expand = c(0, 0, 0.2, 0)) +
          scale_y_discrete(expand = c(0, 1), labels = go_data$Des_short) +
          theme(plot.margin = margin(0.4, 0.4, 0.4, 0.4, 'cm'),
                axis.title = element_text(size = 12, face = 'bold', hjust = 0.5),
                axis.text = element_text(size = 10, face = 'bold', hjust = 0.5),
                legend.title = element_text(size = 11, face = 'bold', hjust = 0),
                plot.title = element_text(size = 14, face = 'bold', hjust = 0.5))
        print(p_up)
        dev.off()
      }

      #> 绘制下调通路图 (NES < 0)
      go_data_down <- gsea_GO_res %>%
        filter(NES < 0) %>%
        mutate(new_NES = -round(NES, digits = 2))

      if(nrow(go_data_down) > 0) {
        go_data <- go_data_down
        if(nrow(go_data) >= GO_num) {
          go_data <- go_data[1:GO_num, ]
        }

        go_data <- go_data %>%
          arrange(new_NES, -p.adjust) %>%
          mutate(Des_short_unique = paste0(Des_short, 1:nrow(go_data))) %>%
          mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique)))

        pdf(file = paste0(photo_dir, 'GSEA_GO_', go_i, '_down_bar.pdf'),
            width = 6.5, height = 4.5)
        p_down <- ggplot(go_data, aes(x = new_NES, y = Des_short_unique, fill = p.adjust)) +
          geom_col(position = "identity") +
          geom_text(aes(label = new_NES), hjust = -0.15) +
          theme_test() +
          labs(x = 'NES', y = '', title = paste0("GSEA GO (Down) of ", go_i, " in ", grp_nm)) +
          scale_fill_gradientn(colors = col) +
          scale_x_continuous(expand = c(0, 0, 0.2, 0)) +
          scale_y_discrete(expand = c(0, 1), labels = go_data$Des_short) +
          theme(plot.margin = margin(0.4, 0.4, 0.4, 0.4, 'cm'),
                axis.title = element_text(size = 12, face = 'bold', hjust = 0.5),
                axis.text = element_text(size = 10, face = 'bold', hjust = 0.5),
                legend.title = element_text(size = 11, face = 'bold', hjust = 0),
                plot.title = element_text(size = 14, face = 'bold', hjust = 0.5))
        print(p_down)
        dev.off()
      }

    }, error = function(e) {
      clog_error(paste0("Error in GSEA GO analysis for ", go_i, ": ", e$message))
    })
  }

  #> GSEA KEGG分析
  clog_normal("Start GSEA KEGG analysis")

  tryCatch({
    gsea_KEGG <- gseKEGG(
      gsea_LogFC,
      organism = kegg_org,
      keyType = "kegg",
      exponent = 1,
      minGSSize = 10,
      maxGSSize = maxGSSize,
      eps = 1e-16,
      pvalueCutoff = ifelse(simplify_index, 0.2, 1),
      pAdjustMethod = "BH",
      verbose = TRUE,
      use_internal_data = internal,
      seed = TRUE,
      by = "fgsea"
    )

    gsea_KEGG <- setReadable(gsea_KEGG, OrgDb = org_index, keyType = "ENTREZID")
    gsea_KEGG_res <- gsea_KEGG@result

    if(nrow(gsea_KEGG_res) > 0) {
      gsea_KEGG_res$Des_short <- sapply(gsea_KEGG_res$Description, .shorten_names)
      gsea_KEGG_res$group <- grp_nm
      results_list[["KEGG"]] <- gsea_KEGG_res

      #> 保存结果
      write.table(gsea_KEGG_res,
                  file = paste0(output_dir, 'GSEA_KEGG_res.txt'),
                  row.names = F, col.names = T, quote = F, sep = '\t')

      #> 绘制上调通路图 (NES > 0)
      kegg_data_up <- gsea_KEGG_res %>%
        filter(NES > 0) %>%
        mutate(new_NES = round(NES, digits = 2))

      if(nrow(kegg_data_up) > 0) {
        kegg_data <- kegg_data_up
        if(nrow(kegg_data) >= GO_num) {
          kegg_data <- kegg_data[1:GO_num, ]
        }

        kegg_data <- kegg_data %>%
          arrange(new_NES, -p.adjust) %>%
          mutate(Des_short_unique = paste0(Des_short, 1:nrow(kegg_data))) %>%
          mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique)))

        pdf(file = paste0(photo_dir, 'GSEA_KEGG_up_bar.pdf'),
            width = 6.5, height = 4.5)
        p_kegg_up <- ggplot(kegg_data, aes(x = new_NES, y = Des_short_unique, fill = p.adjust)) +
          geom_col(position = "identity") +
          geom_text(aes(label = new_NES), hjust = -0.15) +
          theme_test() +
          labs(x = 'NES', y = '', title = paste0("GSEA KEGG (Up) in ", grp_nm)) +
          scale_fill_gradientn(colors = col) +
          scale_x_continuous(expand = c(0, 0, 0.2, 0)) +
          scale_y_discrete(expand = c(0, 1), labels = kegg_data$Des_short) +
          theme(plot.margin = margin(0.4, 0.4, 0.4, 0.4, 'cm'),
                axis.title = element_text(size = 12, face = 'bold', hjust = 0.5),
                axis.text = element_text(size = 10, face = 'bold', hjust = 0.5),
                legend.title = element_text(size = 11, face = 'bold', hjust = 0),
                plot.title = element_text(size = 14, face = 'bold', hjust = 0.5))
        print(p_kegg_up)
        dev.off()
      }

      #> 绘制下调通路图 (NES < 0)
      kegg_data_down <- gsea_KEGG_res %>%
        filter(NES < 0) %>%
        mutate(new_NES = -round(NES, digits = 2))

      if(nrow(kegg_data_down) > 0) {
        kegg_data <- kegg_data_down
        if(nrow(kegg_data) >= GO_num) {
          kegg_data <- kegg_data[1:GO_num, ]
        }

        kegg_data <- kegg_data %>%
          arrange(new_NES, -p.adjust) %>%
          mutate(Des_short_unique = paste0(Des_short, 1:nrow(kegg_data))) %>%
          mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique)))

        pdf(file = paste0(photo_dir, 'GSEA_KEGG_down_bar.pdf'),
            width = 6.5, height = 4.5)
        p_kegg_down <- ggplot(kegg_data, aes(x = new_NES, y = Des_short_unique, fill = p.adjust)) +
          geom_col(position = "identity") +
          geom_text(aes(label = new_NES), hjust = -0.15) +
          theme_test() +
          labs(x = 'NES', y = '', title = paste0("GSEA KEGG (Down) in ", grp_nm)) +
          scale_fill_gradientn(colors = col) +
          scale_x_continuous(expand = c(0, 0, 0.2, 0)) +
          scale_y_discrete(expand = c(0, 1), labels = kegg_data$Des_short) +
          theme(plot.margin = margin(0.4, 0.4, 0.4, 0.4, 'cm'),
                axis.title = element_text(size = 12, face = 'bold', hjust = 0.5),
                axis.text = element_text(size = 10, face = 'bold', hjust = 0.5),
                legend.title = element_text(size = 11, face = 'bold', hjust = 0),
                plot.title = element_text(size = 14, face = 'bold', hjust = 0.5))
        print(p_kegg_down)
        dev.off()
      }
    } else {
      clog_warning("No significant results for KEGG analysis")
    }

  }, error = function(e) {
    clog_error(paste0("Error in GSEA KEGG analysis: ", e$message))
  })

  return(results_list)
}


#' Shorten Long Pathway Names for Visualization
#'
#' Internal function to truncate long pathway names for better display in plots.
#' Shortens names that exceed word count or character length thresholds.
#'
#' @param x Character, pathway name to shorten
#' @param n_word Integer, maximum number of words (default: 6)
#' @param n_char Integer, maximum number of characters (default: 40)
#'
#' @return Shortened pathway name
#'
#' @keywords internal
#'
#' @noRd
.shorten_names <- function(x, n_word=6, n_char=40){ # 单词数>n_word,单词长度>n_char
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > n_char))
  {
    if (nchar(x) > n_char) x <- substr(x, 1, n_char)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  }else{
    return(x)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GOID2Genelist
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Convert GO IDs to Gene Lists
#'
#' Retrieves gene lists associated with specified Gene Ontology (GO) terms.
#' This function extracts genes annotated to given GO IDs from the appropriate
#' organism database.
#'
#' @param GOID Character vector, GO term identifiers to query
#' @param STID_obj STID object (optional, used to extract host organism)
#' @param host_org Character, host organism - "human" or "mouse" (if NULL,
#'        extracted from STID_obj)
#' @param ont Character, GO ontology type - "BP", "CC", "MF", or "ALL"
#'        (default: "BP")
#' @param keyType Character, gene identifier type - "SYMBOL" or "ENTREZID"
#'        (default: "SYMBOL")
#'
#' @return A data frame where columns are GO terms (with names as descriptions)
#'         and rows are genes. Contains gene lists for each requested GO term.
#'
#' @importFrom GOSemSim load_OrgDb
#' @importFrom AnnotationDbi keys mapIds keytypes
#' @import GO.db
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get gene list for specific GO terms
#' go_genes <- GOID2Genelist(
#'   GOID = c("GO:0006955", "GO:0002376"),
#'   host_org = "human",
#'   ont = "BP"
#' )
#'
#' # View the gene list for the first GO term
#' head(go_genes[,1])
#' }
GOID2Genelist <- function(GOID = NULL,
                          STID_obj = NULL, host_org = NULL,
                          ont = "BP", keyType = "SYMBOL"){
  clog_start()

  #> function
  get_GO_data <- function(OrgDb, ont, keytype) {
    GO_Env <- get_GO_Env()
    use_cached <- FALSE

    ont2 <- NULL
    if (exists("ont", envir = GO_Env, inherits = FALSE))
      ont2 <- get("ont", envir = GO_Env)

    if (exists("organism", envir=GO_Env, inherits=FALSE) &&
        exists("keytype", envir=GO_Env, inherits=FALSE) &&
        !is.null(ont2)) {

      org <- get("organism", envir=GO_Env)
      kt <- get("keytype", envir=GO_Env)

      if (org == clusterProfiler:::get_organism(OrgDb) &&
          keytype == kt &&
          (ont == ont2 || ont2 == "ALL") &&
          exists("goAnno", envir=GO_Env, inherits=FALSE)) {

        use_cached <- TRUE
      }
    }

    if (use_cached) {
      goAnno <- get("goAnno", envir=GO_Env)
      if (!is.null(ont2) && ont2 != ont) { ## ont2 == "ALL"
        goAnno <- goAnno[goAnno$ONTOLOGYALL == ont,]
      }
    } else {
      OrgDb <- GOSemSim::load_OrgDb(OrgDb)
      kt <- AnnotationDbi::keytypes(OrgDb)
      if (! keytype %in% kt) {
        stop("keytype is not supported...")
      }

      kk <- AnnotationDbi::keys(OrgDb, keytype=keytype) # 后面没用？？？用于从注释数据库中获取键（keys）

      goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
      if (ont != "ALL") {
        goterms <- goterms[goterms == ont]
      }
      go2gene <- suppressMessages(
        AnnotationDbi::mapIds(OrgDb, keys=names(goterms), column=keytype,
                              keytype="GOALL", multiVals='list')
      )
      goAnno <- stack(go2gene)
      colnames(goAnno) <- c(keytype, "GOALL")
      goAnno <- unique(goAnno[!is.na(goAnno[,1]), ])
      goAnno$ONTOLOGYALL <- goterms[goAnno$GOALL]

      assign("goAnno", goAnno, envir=GO_Env)
      assign("keytype", keytype, envir=GO_Env)
      assign("ont", ont, envir = GO_Env)
      assign("organism", clusterProfiler:::get_organism(OrgDb), envir=GO_Env) #     OrgDb <- load_OrgDb(OrgDb); AnnotationDbi::species(OrgDb)
    }

    GO2GENE <- unique(goAnno[, c(2,1)])
    GO_DATA <- clusterProfiler:::build_Anno(GO2GENE, clusterProfiler:::get_GO2TERM_table())
    goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique

    if (!is.null(ont2) && ont2 == "ALL") {
      return(GO_DATA)
    }

    goOnt <- goOnt.df[,2]
    names(goOnt) <- goOnt.df[,1]
    assign("GO2ONT", goOnt, envir=GO_DATA)

    return(GO_DATA)
  }
  get_GO_Env <- function () {
    if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
      pos <- 1
      envir <- as.environment(pos)
      assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
    }
    get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
  }

  #>
  if(is.null(host_org)){
    if (!inherits(STID_obj, "STID")) {
      clog_error("host_org is NULL and input object is not an STID object")
    }
    host_org <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "host_org") %>% unlist(use.names = F)
  }
  if(host_org == "human"){
    clog_normal("Your host organism is human, loading annotation database for human...")
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      clog_error("Package 'org.Hs.eg.db' is required for human annotation but is not installed.")
    }
    clog_normal(paste0("org.Hs.eg.db version: ", packageVersion("org.Hs.eg.db")))
    OrgDb <- "org.Hs.eg.db"
  }else if(host_org == "mouse"){
    clog_normal("Your host organism is mouse, loading annotation database for mouse...")
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      clog_error("Package 'org.Mm.eg.db' is required for mouse annotation but is not installed.")
    }
    clog_normal(paste0("org.Mm.eg.db version: ", packageVersion("org.Mm.eg.db")))
    OrgDb <- "org.Mm.eg.db"
  }

  #>
  clog_normal("get_GO_data, it may take a while...")
  GO_DATA <- get_GO_data(OrgDb, ont, keyType)
  PATHID2EXTID <- GO_DATA$PATHID2EXTID
  PATHID2NAME <- GO_DATA$PATHID2NAME
  select_GO <- GOID[GOID %in% names(PATHID2EXTID)]
  PATHID2EXTID_sub <- PATHID2EXTID[select_GO]
  names(PATHID2EXTID_sub) <- PATHID2NAME[names(PATHID2EXTID_sub)]
  geneset_df <- do.call(cbind, lapply(lapply(PATHID2EXTID_sub, unlist), `length<-`, max(lengths(PATHID2EXTID_sub))))
  geneset_df <- geneset_df %>% as.data.frame()

  clog_end()
  return(geneset_df)
}



