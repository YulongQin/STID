

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CorrectBackground
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Correct background expression in spatial transcriptomics data
#'
#' This function corrects background expression in spatial transcriptomics data
#' using specified background samples and features. It calculates background gene
#' expression statistics and subtracts them from all samples.
#'
#' @param STID_obj An STID object containing spatial transcriptomics data
#' @param loop_id Must be "LoopAllSamp" (default: "LoopAllSamp")
#' @param bg_samp_id Character vector specifying background sample IDs
#' @param bg_features Character vector specifying background features (genes) for correction
#' @param PosThres_prob Numeric, probability threshold (0-1) for determining positive expression
#'        (default: 0.95)
#' @param adjust_UMI Logical, whether to adjust correction values by mean UMI (default: TRUE)
#' @param assay_id Character, name of the assay to use (default: "Spatial")
#' @param layer_id Character, name of the layer/data slot to use (default: "counts")
#' @param grp_nm Character, group name for output organization (default: NULL, uses timestamp)
#' @param dir_nm Character, directory name for output (default: "M1_CorrectBackground")
#'
#' @return Returns the modified STID object with corrected counts in the Spatial assay
#'
#' @import Seurat
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @import stringr
#' @importFrom tidyr pivot_longer replace_na pivot_wider separate_longer_delim separate_wider_delim fill unite complete drop_na
#' @import parallel
#' @importFrom tibble rownames_to_column column_to_rownames remove_rownames
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Correct background using specified background samples and features
#' STID_obj <- CorrectBackground(
#'   STID_obj = STID_object,
#'   bg_samp_id = c("sample1", "sample2"),
#'   bg_features = c("gene1", "gene2", "gene3"),
#'   PosThres_prob = 0.95,
#'   adjust_UMI = TRUE
#' )
#' }
CorrectBackground <- function(STID_obj = NULL,
                             loop_id = "LoopAllSamp", # must be LoopAllSamp
                             bg_samp_id = NULL, bg_features = NULL,
                             PosThres_prob = 0.95, adjust_UMI = TRUE,
                             assay_id = "Spatial", layer_id = "counts",
                             grp_nm = NULL, dir_nm = "M1_CorrectBackground"
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check( )
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_at_least_one_null(bg_samp_id,bg_features)
  loop_single <- .check_loop_single(STID_obj = STID_obj,loop_id = loop_id, mode = 1)

  #>
  if(!all(bg_samp_id %in% c(loop_single))){
    clog_error(paste0("bg_samp_id must be in: ", paste(loop_single, collapse = ", ")))
  }
  clog_normal(paste0("Your bg_samp_id: ", paste(bg_samp_id, collapse = ", ")))
  if(PosThres_prob <= 0 | PosThres_prob >1){
    clog_error("PosThres_prob must be between 0 and 1")
  }
  clog_normal(paste0("Your PosThres_prob: ", PosThres_prob, " (greater than the value will be retained)"))
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  output_dir_bg <- paste0(output_dir,"bg_samp/")
  output_dir_correct <- paste0(output_dir,"correct_samp/")
  photo_dir_correct <- paste0(photo_dir,"correct_samp/")
  dir.create(output_dir_bg,recursive = T,showWarnings = F)
  dir.create(output_dir_correct,recursive = T,showWarnings = F)
  dir.create(photo_dir_correct,recursive = T,showWarnings = F)
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  clog_step("Start correct background expression")
  if(!is.null(bg_features)){
    valid_features <- .check_features_exist(STID_obj, bg_features)
  }else{
    valid_features <- bg_features
  }
  len_feat <- length(valid_features)
  clog_normal(paste0("Number of background features for correction: ", len_feat))

  # >>> calculate mean UMI
  clog_step("Calculate mean UMI of background samples and all samples")
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  meta_data <- STID_obj@meta.data
  count_assay <- GetAssayData(STID_obj, assay = assay_id, layer = layer_id )
  meta_data$nCount_RNA <- Matrix::colSums(count_assay)
  # meta_data$nFeature_RNA <- Matrix::colSums(GetAssayData(object = STID_obj, layer = "counts") > 0)
  mean_UMI_samp <- meta_data %>%
    group_by(!!sym(samp_colnm)) %>%
    summarise(mean_UMI = mean(nCount_RNA, na.rm = T)) %>%
    as.data.frame() %>%
    `row.names<-`(.[[samp_colnm]])
  mean_UMI_C <- mean_UMI_samp %>%
    filter(!!sym(samp_colnm) %in% bg_samp_id) %>%
    summarise(mean_UMI_C = mean(mean_UMI, na.rm = T)) %>%
    pull(mean_UMI_C)
  mean_UMI_all <- mean(mean_UMI_samp$mean_UMI, na.rm = T)
  clog_normal(paste0("Mean UMI of background samples: ", round(mean_UMI_C,2)))
  clog_normal(paste0("Mean UMI of all samples: ", round(mean_UMI_all,2)))

  # >>>
  clog_step("Calculate background gene expression statistics")
  res_df_C <- as.data.frame(matrix(NA,nrow = length(valid_features),ncol = length(bg_samp_id)))
  colnames(res_df_C) <- bg_samp_id
  rownames(res_df_C) <- valid_features
  for(i in 1:length(bg_samp_id)){
    # i = 1
    i_single <- bg_samp_id[i]
    clog_loop(paste0("Processing background sample: ", i_single," (",i,"/",length(bg_samp_id),")"))
    i_samp_meta_data <- meta_data[meta_data[[samp_colnm]] == i_single, ]
    i_samp_cells <- rownames(i_samp_meta_data)
    i_count_assay <- count_assay[valid_features,i_samp_cells, drop = FALSE]
    # table((rownames(i_count_assay) == valid_features))
    clog_normal(paste0("Dimension of count matrix: ", paste(dim(i_count_assay), collapse = " x ")))

    # >
    bg_gene_stat_df <- as.data.frame(matrix(NA,nrow = length(valid_features),ncol = 5))
    colnames(bg_gene_stat_df) <- c("nonzero_95","nonzero_95_adjust","nonzero_num","nonzero_ratio","nonzero_count")
    rownames(bg_gene_stat_df) <- valid_features
    bg_gene_stat_df$nonzero_95 <- apply(i_count_assay, 1, function(x) {
      quantile(x[x>0], probs = PosThres_prob , na.rm = T)
    })
    bg_gene_stat_df$nonzero_95_adjust <- bg_gene_stat_df$nonzero_95 * (mean_UMI_all / mean_UMI_C)
    bg_gene_stat_df$nonzero_num <- rowSums(i_count_assay>0)
    bg_gene_stat_df$nonzero_ratio <- bg_gene_stat_df$nonzero_num / ncol(i_count_assay)
    bg_gene_stat_df$nonzero_count <- rowSums(i_count_assay)
    bg_gene_stat_df <- bg_gene_stat_df %>%
      mutate(across(everything(), ~ replace_na(., 0)))
    fwrite(bg_gene_stat_df, file = paste0(output_dir_bg,i_single,"_bg_gene_stat_(CorrectBackground).txt"),
           sep = "\t", quote = F,col.names = T,row.names = T,na = "NA")
    if(adjust_UMI){
      res_df_C[,i] <- bg_gene_stat_df$nonzero_95_adjust
    }else{
      res_df_C[,i] <- bg_gene_stat_df$nonzero_95
    }
  }
  res_df_C  <- rowMeans(res_df_C, na.rm = T)
  names(res_df_C) <- valid_features

  # >>> Start correct background
  clog_step("Start correct background expression for all samples")
  corrected_objects <- list()
  plot_histgram_list <- list()
  plot_barplot_list <- list()
  for(i in 1:length(loop_single)){
    i_single <- loop_single[i]
    clog_loop(paste0("Processing sample: ", i_single," (",i,"/",length(loop_single),")"))
    i_count_assay <- count_assay[,meta_data[[samp_colnm]] == i_single, drop = FALSE]
    assay_pathogen_before <- i_count_assay[valid_features, , drop = FALSE]
    clog_normal(paste0("Dimension of count matrix: ", paste(dim(i_count_assay), collapse = " x ")))
    clog_normal(paste0("Number of common genes for background correction: ", length(valid_features)))
    i_mean_UMI <- mean_UMI_samp[i_single, "mean_UMI"]
    clog_normal(paste0("Mean UMI for ", i_single, ": ", round(i_mean_UMI, 2)))
    clog_normal(paste0("Mean UMI for all samples: ", round(mean_UMI_all, 2)))
    correction_vec <- rep(0, nrow(i_count_assay))
    names(correction_vec) <- rownames(i_count_assay)
    correction_vec[names(res_df_C)] <- res_df_C
    clog_normal(paste0("Max correction value before scaling: ", round(max(correction_vec, na.rm = TRUE), 2)))
    if(adjust_UMI){
      clog_normal("Adjust correction values by mean UMI")
      scaled_correction <- (correction_vec / mean_UMI_all) * i_mean_UMI
    }else{
      clog_normal("No adjust correction values by mean UMI")
      scaled_correction <- correction_vec
    }
    clog_normal(paste0("Max correction value after scaling: ", round(max(scaled_correction, na.rm = TRUE), 2)))

    # >
    mat <- as(i_count_assay, "dgCMatrix")
    non_zero_values <- mat@x
    row_indices_1based <- mat@i + 1
    corrected_values <- non_zero_values - scaled_correction[row_indices_1based]
    corrected_values[corrected_values < 0] <- 0
    mat@x <- corrected_values %>% round(0)
    assay_pathogen_after <- mat[valid_features, , drop = FALSE]
    corrected_objects[[i_single]] <- mat

    #> plot distribution
    stat_pathogen <- data.frame(
      TotalCount_before = rowSums(assay_pathogen_before),
      TotalCount_after = rowSums(assay_pathogen_after)
    )
    write.table(stat_pathogen,
                file = paste0(output_dir_correct,i_single,"_BgCorrect_stat_(CorrectBackground).txt"),
                sep = "\t", quote = F,col.names = NA,row.names = T,na = "NA")
    plot_data <- stat_pathogen %>%
      filter(rowSums(.)>0) %>%
      pivot_longer(
        cols = everything(),
        names_to = "grp",
        values_to = "count"
      ) %>%
      mutate( grp = factor(grp, levels = c("TotalCount_before","TotalCount_after")))
    p1 <- ggplot(plot_data, aes( x = count, fill = grp,color = grp)) +
      # geom_density( alpha = 0.5,bw = 10,trim = T)
      geom_histogram(aes(y = after_stat(count)), bins = 20, alpha = 0.8,color = "grey95",linewidth = 0.25) +
      geom_text(aes(label=ifelse(after_stat(count) == 0, "", after_stat(count)),
                    y=after_stat(count)), stat='bin', bins = 20, vjust=-0.5, size=3, color='black') +
      facet_wrap(~grp, ncol = 1) +
      scale_fill_manual(values = c("#99CCFF", "#FF9999")) +
      labs(
        title = paste0("Background correction for sample: ", i_single),
        x = "Pathogen gene counts",
        y = "Number of genes",
        fill = NULL,
        color = NULL
      ) +
      theme_test() +
      scale_y_continuous(trans = "log10",expand = expansion(mult = c(0.05, 0.15))) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_rect(fill="grey95"))
    ggsave(p1, filename = paste0(photo_dir_correct,i_single,"_BgCorrect_histogram_(CorrectBackground).pdf"),
           width = 6, height = 7,limitsize = FALSE)
    plot_histgram_list[[i]] <- p1

    #> plot count ratio
    stat_pathogen_ratio <- data.frame(
      grp = c("Before","After"),
      count = c(sum(stat_pathogen$TotalCount_before), sum(stat_pathogen$TotalCount_after))
    )
    write.table(stat_pathogen_ratio,
                file = paste0(output_dir_correct,i_single,"_BgCorrect_totalcount_(CorrectBackground).txt"),
                sep = "\t", quote = F,col.names = NA,row.names = T,na = "NA")
    plot_data <- stat_pathogen_ratio %>%
      mutate( grp = factor(grp, levels = c("Before","After")))
    p2 <- ggplot(plot_data, aes(x = grp, y = count, fill = grp)) +
      geom_bar(stat = "identity", color = "grey95", alpha = 0.8,linewidth = 0.25) +
      geom_text(aes(label = count), vjust = -0.5, size = 5, color = 'black') +
      scale_fill_manual(values = c("#99CCFF", "#FF9999")) +
      labs(
        title = i_single,
        x = NULL,
        y = "Total counts",
        fill = NULL
      ) +
      theme_test() +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_rect(fill="grey95"))
    ggsave(p2, filename = paste0(photo_dir_correct,i_single,"_BgCorrect_totalcount_barplot_(CorrectBackground).pdf"),
           width = 2, height = 4,limitsize = FALSE)
    plot_barplot_list[[(2*i-1)]] <- p2

    #> plot spot ratio
    pathogen_spot_ratio <- data.frame(
      grp = c("Before","After"),
      spot_num = c(sum(colSums(assay_pathogen_before)>0), sum(colSums(assay_pathogen_after)>0))
    )
    write.table(pathogen_spot_ratio,
                file = paste0(output_dir_correct,i_single,"_BgCorrect_PosSpotNum_(CorrectBackground).txt"),
                sep = "\t", quote = F,col.names = NA,row.names = T,na = "NA")
    plot_data <- pathogen_spot_ratio %>%
      mutate( grp = factor(grp, levels = c("Before","After")))
    p3 <- ggplot(plot_data, aes(x = grp, y = spot_num, fill = grp)) +
      geom_bar(stat = "identity", color = "grey95", alpha = 0.8,linewidth = 0.25) +
      geom_text(aes(label = spot_num), vjust = -0.5, size = 5, color = 'black') +
      scale_fill_manual(values = c("#99CCFF", "#FF9999")) +
      labs(
        title = i_single,
        x = NULL,
        y = "Total positive spots",
        fill = NULL
      ) +
      theme_test() +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_rect(fill="grey95"))
    ggsave(p3, filename = paste0(photo_dir_correct,i_single,"_BgCorrect_PosSpotNum_barplot_(CorrectBackground).pdf"),
           width = 2, height = 4,limitsize = FALSE)
    plot_barplot_list[[2*i]] <- p3
  }
  p_histgram_all <- wrap_plots(
    plot_histgram_list,
    ncol = 3,
  )
  p_barplot_all <- wrap_plots(
    plot_barplot_list,
    ncol = 6,
  )
  ggsave(p_histgram_all, filename = paste0(photo_dir,"All_samp_BgCorrect_histogram_(CorrectBackground).pdf"),
         width = 18, height = 7*ceiling(length(loop_single)/3),limitsize = FALSE)
  ggsave(p_barplot_all, filename = paste0(photo_dir,"All_samp_BgCorrect_barplot_(CorrectBackground).pdf"),
         width = 12, height = 4*ceiling(length(loop_single)/3),limitsize = FALSE)

  # >>> Save corrected data
  clog_step("Save background corrected data")
  corrected_count_assay <- do.call(cbind, corrected_objects)
  STID_obj@assays[["Spatial"]]["counts"] <- corrected_count_assay

  # >>>
  .save_function_params("CorrectBackground", envir = environment(),
                        file = paste0(output_dir,"Log_function_params_(CorrectBackground).log") )
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(SpotDetect_Geneset).log")) %>% invisible()
  return(STID_obj)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SpotDetect_Gene
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Detect infection spots based on individual gene expression
#'
#' This function identifies positive spots based on expression thresholds for
#' individual genes. It can apply spatial smoothing and generates visualization
#' of gene expression patterns.
#'
#' @param STID_obj An STID object containing spatial transcriptomics data
#' @param features Character vector of gene names to analyze
#' @param feature_colnm Character, column name in metadata containing pre-computed features
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param assay_id Character, name of the assay to use (default: "Spatial")
#' @param layer_id Character, name of the layer/data slot to use (default: "counts")
#' @param PosThres_prob Numeric, probability threshold (0-1) for positive detection
#'        (default: 0)
#' @param PosThres_count Numeric, absolute count threshold for positive detection
#'        (default: 0)
#' @param PosThres_gene Numeric, gene number threshold (default: 1)
#' @param col Color palette for visualization (default: COLOR_DIS_CON)
#' @param pt_size Numeric, point size in spatial plots (default: 0.5)
#' @param vmin Numeric or character, minimum value for color scale (default: NULL)
#' @param vmax Numeric or character, maximum value for color scale (default: "p99")
#' @param black_bg Logical, whether to use black background in plots (default: FALSE)
#' @param plot_method Character, plotting mode - "merge" or "single" (default: "merge")
#' @param blur_method Character, spatial smoothing method - "isoblur" or "medianblur"
#'        (default: NULL)
#' @param blur_n Numeric, number of iterations for median blur (default: 1)
#' @param blur_sigma Numeric, sigma parameter for isoblur (default: 0.25)
#' @param description Character, description of the analysis (default: NULL)
#' @param grp_nm Character, group name for output organization (default: NULL)
#' @param dir_nm Character, directory name for output (default: "M1_SpotDetect_Gene")
#'
#' @return Returns the modified STID object with added metadata containing
#'         gene expression labels and thresholds
#'
#' @import Seurat
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @import stringr
#' @importFrom furrr future_map furrr_options
#' @importFrom future nbrOfWorkers plan sequential multisession
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Detect spots based on pathogen gene expression
#' STID_obj <- SpotDetect_Gene(
#'   STID_obj = STID_object,
#'   features = c("geneA", "geneB", "geneC"),
#'   PosThres_prob = 0.95,
#'   PosThres_count = 1,
#'   plot_method = "merge"
#' )
#' }
SpotDetect_Gene <- function(
    STID_obj = NULL, features = NULL, feature_colnm = NULL,
    loop_id= "LoopAllSamp", assay_id = "Spatial", layer_id = "counts",
    PosThres_prob = 0, PosThres_count = 0, PosThres_gene = 1,
    col = COLOR_DIS_CON, pt_size = 0.5, vmin = NULL, vmax = "p99",
    black_bg = FALSE, plot_method = "merge",
    blur_method = NULL, blur_n = 1, blur_sigma = 0.25,
    description = NULL,
    grp_nm = NULL,dir_nm = "M1_SpotDetect_Gene"
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check( )
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_at_least_one_null(feature_colnm,features)
  match.arg(plot_method, .STID_globals$plot_methods)

  #>
  loop_single <- .check_loop_single(STID_obj = STID_obj,loop_id = loop_id, mode = 1)
  if(PosThres_prob <0 | PosThres_prob >1){
    clog_error("PosThres_prob must be between 0 and 1")
  }
  if(PosThres_count <0){
    clog_error("PosThres_count must be >=0")
  }
  basic_thres_value <- 0
  clog_normal(paste0("Your PosThres_prob: ", PosThres_prob, " (greater than the value will be retained)"))
  clog_normal(paste0("Your PosThres_count: ", PosThres_count, " (greater than the value will be retained)"))
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  clog_step("Start detect infection spot by specific genes")
  meta_data <- GetMetaData(STID_obj = STID_obj, meta_key = "raw")[[1]]
  fwrite(meta_data, file = paste0(output_dir,"Spot_meta.data_(SpotDetect_Gene).txt"),
         sep = "\t", quote = F,col.names = T,row.names = T)
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  x_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "x_colnm")[[1]]
  y_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "y_colnm")[[1]]
  if(!is.null(features)){
    valid_features <- .check_features_exist(STID_obj, features)
    feat_exp <- FetchData(STID_obj, vars = valid_features, assay = assay_id, layer = layer_id)
    meta_data <- cbind(meta_data[rownames(feat_exp),],feat_exp)
  }else{
    valid_features <- features
  }
  if(!is.null(feature_colnm)){
    if(all(feature_colnm %in% colnames(meta_data))){
      valid_features <- c(valid_features, feature_colnm)
    }
  }
  len_feat <- length(valid_features)
  fwrite(meta_data[,c(samp_colnm,x_colnm,y_colnm,valid_features)],
         file = paste0(output_dir,"Spot_GeneEXP_(SpotDetect_Gene).txt"),
         sep = "\t", quote = F,col.names = T,row.names = T)

  #> blur
  if(!is.null(blur_method)){
    match.arg(blur_method, .STID_globals$blur_methods)
    clog_step(paste0("Start spatial smooth for features by mode: ", blur_method))
    for(i in 1:length(loop_single)){
      i_single <- loop_single[i]
      clog_normal(paste0("Processing blur for sample: ", i_single," (",i,"/",length(loop_single),")"))
      i_samp_meta_data <- meta_data[meta_data[[samp_colnm]] == i_single, c(x_colnm,y_colnm,valid_features)]

      for(j in 1:len_feat){
        j_feat <- valid_features[j]
        # clog_normal(paste0("Processing blur for feature: ", j_feat," (",j,"/",len_feat,")"))
        j_feat_exp <- i_samp_meta_data[c("x","y",j_feat)]
        j_feat_exp_smooth <- .SpatialSmooth(meta_data = j_feat_exp, blur_method = blur_method, blur_n = blur_n, blur_sigma = blur_sigma)
        meta_data[meta_data[[samp_colnm]] == i_single, j_feat] <- j_feat_exp_smooth[,j_feat]
      }
    }
    fwrite(meta_data[,c(samp_colnm,x_colnm,y_colnm,valid_features)],
           file = paste0(output_dir,"Spot_GeneEXP_blur_(SpotDetect_Gene).txt"),
           sep = "\t", quote = F,col.names = T,row.names = T)
  }

  #> meta2STID
  new_meta_key <- paste0(dir_nm,"_",grp_nm)
  meta2STID <- meta_data[,c(samp_colnm,x_colnm,y_colnm,valid_features),drop = F]
  meta2STID_Label <- meta2STID %>%
    dplyr::select(all_of(valid_features)) %>%
    mutate( across(everything(), ~ NA_real_)) %>%
    rename_with( ~ paste0("Label_", .x), .cols = all_of(valid_features))
  meta2STID_Label <- meta2STID %>%
    bind_cols(meta2STID_Label) %>%
    as.matrix()

  # >>> meta_data_longer
  meta_data_longer <- meta2STID %>%
    mutate(Spot_id = rownames(.),.before = 1) %>%
    pivot_longer(
      cols = all_of(valid_features),
      names_to = "feature",
      values_to = "Expression"
    ) %>%
    group_by(!!sym(samp_colnm), feature) %>%
    mutate(basic_thres_value = basic_thres_value, # !!!!
           PosThres_count = PosThres_count,
           PosThres_prob = PosThres_prob) %>%
    mutate(PosThres_prob_value = if(unique(PosThres_prob) >0) quantile(Expression[Expression > basic_thres_value],PosThres_prob) else NA, # !!!
           max_thres_value = if(is.na(unique(PosThres_prob_value))) PosThres_count else ifelse(PosThres_prob_value >= PosThres_count, PosThres_prob_value, PosThres_count),
           max_thres_value = ifelse(max_thres_value>basic_thres_value, max_thres_value, basic_thres_value),
           exp_pos = ifelse(Expression > max_thres_value, Expression, basic_thres_value)) %>%
    ungroup()
  Sap_GeneExp <- meta_data_longer %>%
    group_by(!!sym(samp_colnm), feature) %>%
    summarise(
      all_num = n(),
      all_mean = mean(Expression,na.rm = T),
      all_sd = sd(Expression,na.rm = T),
      all_Q50 = quantile(Expression, 0.50,na.rm = T),
      all_Q75 = quantile(Expression, 0.75,na.rm = T),
      all_Q95 = quantile(Expression, 0.95,na.rm = T),
      all_Q99 = quantile(Expression, 0.99,na.rm = T),
      all_max = max(Expression,na.rm = T),
      basic_thres = unique(basic_thres_value),
      pos_thres = unique(max_thres_value),
      pos_num = sum(Expression > max_thres_value,na.rm = T),
      pos_ratio = pos_num / all_num ,
      pos_mean = mean(Expression[Expression > max_thres_value],na.rm = T),
      pos_sd = sd(Expression[Expression > max_thres_value],na.rm = T),
      pos_Q5 = quantile(Expression[Expression > max_thres_value], 0.05,na.rm = T),
      pos_Q25 = quantile(Expression[Expression > max_thres_value], 0.25,na.rm = T),
      pos_Q50 = quantile(Expression[Expression > max_thres_value], 0.50,na.rm = T),
      pos_Q75 = quantile(Expression[Expression > max_thres_value], 0.75,na.rm = T),
      pos_Q95 = quantile(Expression[Expression > max_thres_value], 0.95,na.rm = T),
    )
  fwrite(Sap_GeneExp, file = paste0(output_dir,"Sap_GeneExp_stat_(SpotDetect_Gene).txt"),
         sep = "\t", quote = F,col.names = T,row.names = F,na = "NA")

  Sap_thres <- meta_data_longer %>%
    dplyr::select(!!sym(samp_colnm), feature, basic_thres_value,PosThres_prob, PosThres_prob_value, PosThres_count, max_thres_value) %>%
    distinct()
  fwrite(Sap_thres, file = paste0(output_dir,"Sap_GeneExp_threshold_(SpotDetect_Gene).txt"),
         sep = "\t", quote = F,col.names = T,row.names = F,na = "NA")

  # >>> Plot_Spatial
  pathogen_org <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "pathogen_org")[[1]]
  host_org <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "host_org")[[1]]
  pathogen_grp <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "pathogen_grp")[[1]]
  subtitle_lable <- paste0(pathogen_grp,"_",pathogen_org," (",host_org,")")

  #> process_samp function
  process_samp <- function(i){
    i_single <- loop_single[i]
    clog_loop(paste0("Processing sample: ", i_single," (",i,"/",length(loop_single),")"))
    i_output_dir <- paste0(output_dir,i_single,"/")
    i_photo_dir <- paste0(photo_dir,i_single,"/")
    dir.create(i_output_dir,recursive = T, showWarnings =F)
    dir.create(i_photo_dir,recursive = T, showWarnings = F)
    i_meta_data_longer <- meta_data_longer[meta_data_longer[[samp_colnm]] == i_single,]
    i_meta2STID_Label <- meta2STID_Label[meta_data[[samp_colnm]] == i_single, , drop = FALSE]

    #> .Plot_histgram
    p_hist <- .Plot_histgram(plot_data = i_meta_data_longer, samp_colnm = samp_colnm,col = "#99C5E3",
                             title = paste0("Gene expression distribution in ", i_single),
                             subtitle = subtitle_lable)
    if(len_feat <= 4){
      width_value <- 4*len_feat
    }else{
      width_value <- 12
    }
    ggsave(p_hist, filename = paste0(i_photo_dir,i_single,"_histgram_(SpotDetect_Gene).pdf"),
           width = width_value, height = 3.5*ceiling(len_feat/4) + 0.5,limitsize = FALSE)

    #> Plot_Spatial
    if(plot_method == "single"){
      dir.create(paste0(i_photo_dir,"/SpatialPlot_samp/"),recursive = T, showWarnings = F)
      plot_dis_list <- list()
      plot_con_list <- list()
      plot_SpaialPlot_list <- list()

      for(j in 1:len_feat){
        j_feat <- valid_features[j]
        clog_normal(paste0("Processing feature: ", j_feat," (",j,"/",len_feat,")"))
        j_meta_data_longer <- i_meta_data_longer[i_meta_data_longer$feature == j_feat,]
        j_PosThres_prob_value <- j_meta_data_longer$PosThres_prob_value %>% unique()
        j_PosThre_count <- j_meta_data_longer$PosThres_count %>% unique()
        j_basic_thres_value <- j_meta_data_longer$basic_thres_value %>% unique()
        j_max_thres_value <- j_meta_data_longer$max_thres_value %>% unique()
        clog_normal(paste0("basic_thres_value: ", round(j_basic_thres_value,3),
                           "; PosThres_count: ", round(j_PosThre_count,3),
                           "; PosThres_prob_value: ", round(j_PosThres_prob_value,3),
                           "; max_thres_value: ", j_max_thres_value %>% signif(3)))
        subtitle_lable_pos <- paste0(subtitle_lable,"\n","Pos_Thres: ", round(j_max_thres_value,3))

        # >
        j_heat_gene <- j_meta_data_longer$Expression
        title_lable <- paste0(j_feat," (",i_single,")")
        plot_SpaialPlot <- j_meta_data_longer %>%
          mutate(
            Expression = as.numeric(Expression),
            Label = ifelse(Expression > max_thres_value, "pos","neg"),
            Label = factor(Label, levels = c("neg","pos"))
          )
        plot_SpaialPlot_list[[j]] <- plot_SpaialPlot
        fwrite(plot_SpaialPlot, file = paste0(i_output_dir,i_single,"_",j_feat,"_PosInfo_(SpotDetect_Gene).txt"), # !!!
               sep = "\t", quote = F,col.names = T,row.names = F)

        #> coord_Label
        coord_Label <- plot_SpaialPlot[,c("Spot_id","feature")] %>%
          mutate(feature = paste0("Label_", feature)) %>%
          as.matrix()
        i_meta2STID_Label[coord_Label] <- plot_SpaialPlot$Label %>%
          as.character()

        # >
        P_spatial_dis <- Plot_Spatial(plot_data = plot_SpaialPlot, x_colnm = x_colnm, y_colnm = y_colnm, group_by = "Label",
                                      facet_grpnm = NULL, datatype = "discrete",
                                      col = col,pt_size = pt_size, vmin = vmin, vmax = vmax,
                                      title = title_lable, subtitle = subtitle_lable_pos,
                                      black_bg = black_bg)
        P_spatial_con <- Plot_Spatial(plot_data = plot_SpaialPlot, x_colnm = x_colnm, y_colnm = y_colnm, group_by = "Expression",
                                      facet_grpnm = NULL, datatype = "continuous",
                                      col = col,pt_size = pt_size, vmin = vmin, vmax = vmax,
                                      title = title_lable, subtitle = subtitle_lable,
                                      black_bg = black_bg)
        ggsave(P_spatial_dis, filename = paste0(i_photo_dir,"/SpatialPlot_samp/",i_single,"_",j_feat,"_dis_(SpotDetect_Gene).pdf"),
               width = 5, height = 5,limitsize = FALSE)
        ggsave(P_spatial_con, filename = paste0(i_photo_dir,"/SpatialPlot_samp/",i_single,"_",j_feat,"_con_(SpotDetect_Gene).pdf"),
               width = 5, height = 5,limitsize = FALSE)
        plot_dis_list[[j]] <- P_spatial_dis
        plot_con_list[[j]] <- P_spatial_con
      }

      clog_normal(paste0("Save all features plot to a merge pdf"))
      p_dis_all <- wrap_plots(
        plot_dis_list,
        ncol = 4,
      ) +
        plot_annotation() +
        theme(
          plot.background = element_rect(
            fill = ifelse(black_bg, "black", "white"),
            color = ifelse(black_bg, "black", "white")
          )
        )
      p_con_all <- wrap_plots(
        plot_con_list,
        ncol = 4,
      ) +
        plot_annotation() +
        theme(
          plot.background = element_rect(
            fill = ifelse(black_bg, "black", "white"),
            color = ifelse(black_bg, "black", "white")
          )
        )
      fwrite(do.call(rbind, plot_SpaialPlot_list),
             file = paste0(i_output_dir,i_single,"_allgene_PosInfo_(SpotDetect_Gene).txt"),
             sep = "\t", quote = F,col.names = T,row.names = F)

    }else if(plot_method == "merge"){
      plot_SpaialPlot <- i_meta_data_longer %>%
        mutate(
          Expression = as.numeric(Expression),
          Label = ifelse(Expression >max_thres_value, "pos","neg"),
          Label = factor(Label, levels = c("neg","pos"))
        )
      fwrite(plot_SpaialPlot, file = paste0(i_output_dir,i_single,"_allgene_PosInfo_(SpotDetect_Gene).txt"),
             sep = "\t", quote = F,col.names = T,row.names = F)

      #> coord_Label
      coord_Label <- plot_SpaialPlot[,c("Spot_id","feature")] %>%
        mutate(feature = paste0("Label_", feature)) %>%
        as.matrix()
      i_meta2STID_Label[coord_Label] <- plot_SpaialPlot$Label %>%
        as.character()

      # >
      p_dis_all <- Plot_Spatial(plot_data = plot_SpaialPlot, x_colnm = x_colnm, y_colnm = y_colnm, group_by = "Label",
                                facet_grpnm = "feature", datatype = "discrete",
                                col = col,pt_size = pt_size, vmin = vmin, vmax = vmax,
                                title = paste0("Gene expression in ", i_single), subtitle = subtitle_lable,
                                black_bg = black_bg)
      p_con_all <- Plot_Spatial(plot_data = plot_SpaialPlot, x_colnm = x_colnm, y_colnm = y_colnm, group_by = "Expression",
                                facet_grpnm = "feature", datatype = "continuous",
                                col = col,pt_size = pt_size, vmin = vmin, vmax = vmax,
                                title = paste0("Gene expression in ", i_single), subtitle = subtitle_lable,
                                black_bg = black_bg)
    }
    ggsave(p_dis_all, filename = paste0(i_photo_dir,i_single,"_all_dis_(SpotDetect_Gene).pdf"),
           width = 20, height = 5*ceiling(len_feat/4),limitsize = FALSE)
    ggsave(p_con_all, filename = paste0(i_photo_dir,i_single,"_all_con_(SpotDetect_Gene).pdf"),
           width = 20, height = 5*ceiling(len_feat/4),limitsize = FALSE)

    return(i_meta2STID_Label)
  }

  #> parallel or not
  meta2STID_Label_list <- list()
  parallel_workers <- getOption("parallel_workers")
  default_plan <- getOption("default_plan")
  if(parallel_workers > 1){
    plan(default_plan, workers = parallel_workers)
    clog_step(paste0("Start parallel processing sample with ", nbrOfWorkers(), " workers"))
    meta2STID_Label_list <- future_map(seq_along(loop_single), process_samp,
                                      .progress=T, .options = furrr_options(seed = 123))
    names(meta2STID_Label_list) <- loop_single
    plan(sequential)
  }else{
    clog_step("Start processing sample without parallelization")
    plan(sequential)
    for(i in seq_along(loop_single)){
      meta2STID_Label_list[[i]] <- process_samp(i)
    }
  }

  #> AddMetaData
  meta2STID <- lapply(meta2STID_Label_list, as.data.frame) %>%
    bind_rows() %>%
    mutate(across(c(x,y,all_of(valid_features)), as.numeric))
  meta2STID <- meta2STID[rownames(meta_data), , drop = FALSE] # keep order
  STID_obj <- AddMetaData(STID_obj = STID_obj,
                         meta_key = new_meta_key,
                         add_data = meta2STID,
                         dir_nm = dir_nm,
                         grp_nm = grp_nm,
                         asso_key = NULL,
                         description = description)

  # >>> Final
  .save_function_params("SpotDetect_Gene", envir = environment(), file = paste0(output_dir,"Log_function_params_(SpotDetect_Gene).log") )
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(SpotDetect_Gene).log")) %>% invisible()
  return(STID_obj)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SpotDetect_Geneset
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Detect infection spots based on gene set scores
#'
#' This function calculates scores for gene sets using various methods and
#' identifies positive spots based on score thresholds. Supports multiple scoring
#' algorithms including AddModuleScore, AUCell, UCell, and basic statistics.
#'
#' @param STID_obj An STID object containing spatial transcriptomics data
#' @param geneset_list A named list of gene sets for scoring
#' @param score_method Character, scoring method to use. Options: "AddModuleScore",
#'        "AUCell", "UCell", "MeanExp", "SumExp" (default: "AddModuleScore")
#' @param n_iter Numeric, number of iterations for AddModuleScore (default: 5)
#' @param nbin Numeric, number of bins for AddModuleScore (default: 10)
#' @param seed Numeric, random seed for reproducibility (default: 10)
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param assay_id Character, name of the assay to use (default: "Spatial")
#' @param layer_id Character, name of the layer/data slot to use (default: "data")
#' @param PosThres_prob Numeric, probability threshold (0-1) for positive detection
#'        (default: 0)
#' @param PosThres_score Numeric, absolute score threshold for positive detection
#'        (default: 0)
#' @param col Color palette for visualization (default: COLOR_DIS_CON)
#' @param pt_size Numeric, point size in spatial plots (default: 0.5)
#' @param vmin Numeric or character, minimum value for color scale (default: NULL)
#' @param vmax Numeric or character, maximum value for color scale (default: "p99")
#' @param black_bg Logical, whether to use black background in plots (default: FALSE)
#' @param plot_method Character, plotting mode - "merge" or "single"
#'        (default: "merge")
#' @param blur_method Character, spatial smoothing method - "isoblur" or "medianblur"
#'        (default: NULL)
#' @param blur_n Numeric, number of iterations for median blur (default: 1)
#' @param blur_sigma Numeric, sigma parameter for isoblur (default: 0.25)
#' @param description Character, description of the analysis (default: NULL)
#' @param grp_nm Character, group name for output organization (default: NULL)
#' @param dir_nm Character, directory name for output (default: "M1_SpotDetect_Geneset")
#'
#' @return Returns the modified STID object with added metadata containing
#'         gene set scores and labels
#'
#' @import Seurat
#' @import AUCell
#' @import UCell
#' @import Matrix
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @import stringr
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define gene sets
#' pathogen_genes <- list(
#'   "Virulence" = c("gene1", "gene2", "gene3"),
#'   "Toxin" = c("gene4", "gene5", "gene6")
#' )
#'
#' # Detect spots based on gene set scores
#' STID_obj <- SpotDetect_Geneset(
#'   STID_obj = STID_object,
#'   geneset_list = pathogen_genes,
#'   score_method = "AddModuleScore",
#'   PosThres_prob = 0.95
#' )
#' }
SpotDetect_Geneset <- function(
    STID_obj = NULL, geneset_list = NULL,
    score_method = "AddModuleScore", n_iter = 5, nbin = 10, seed = 10,
    loop_id= "LoopAllSamp", assay_id = "Spatial", layer_id = "data",
    PosThres_prob = 0, PosThres_score = 0,
    col = COLOR_DIS_CON,
    pt_size = 0.5, vmin = NULL, vmax = "p99",
    black_bg = FALSE, plot_method = "merge",
    blur_method = NULL, blur_n = 1, blur_sigma = 0.25,
    description = NULL,
    grp_nm = NULL,dir_nm = "M1_SpotDetect_Geneset"
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check( )
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(geneset_list)
  match.arg(score_method, .STID_globals$score_methods)
  match.arg(plot_method, .STID_globals$plot_methods)
  if(!is.list(geneset_list)){
    clog_error("geneset_list must be a list")
  }

  # >
  loop_single <- .check_loop_single(STID_obj = STID_obj,loop_id = loop_id, mode = 1)
  if(PosThres_prob <0 | PosThres_prob >1){
    clog_error("PosThres_prob must be between 0 and 1")
  }
  if(PosThres_score <0){
    clog_error("PosThres_score must be >=0")
  }
  clog_normal(paste0("Your PosThres_prob: ", PosThres_prob, " (greater than the value will be retained)"))
  clog_normal(paste0("Your PosThres_score: ", PosThres_score, " (greater than the value will be retained)"))
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start score
  clog_step("Start calculate geneset score")
  geneset_list <- lapply(geneset_list, function(x){
    x[nzchar(x) & !is.na(x)]
  })
  len_geneset <- length(geneset_list)
  meta_data_raw <- STID_obj@meta.data

  if(score_method %in% c("AddModuleScore","UCell")){
    if(score_method == "Ucell"){
      clog_warn("UCell method is don't need n_iter, set n_iter to 1")
      n_iter <- 1
    }
    seed_vec <- seed + seq_len(n_iter)*10
    for (i in 1:n_iter) {
      clog_normal(paste0("Processing n_iter: ",i,"/",n_iter))
      if(score_method == "UCell"){
        clog_normal("Calculate UCell score..., it may take a while")
        score_meta_data <- AddModuleScore_UCell(
          obj = STID_obj,
          features = geneset_list,
          name = "_Infect_geneset_123",
          assay = assay_id,
          slot = layer_id)@meta.data
      }else if(score_method == "AddModuleScore"){
        score_meta_data <- AddModuleScore(object = STID_obj, features = geneset_list, name = "Infect_geneset_123_",
                                          nbin = 10, seed = seed_vec[i], assay = assay_id,layer = layer_id)@meta.data
      }
      if(i ==1){
        score_df <- score_meta_data %>%
          dplyr::select(contains("Infect_geneset_123")) %>%
          dplyr::rename_with(.fn = ~names(geneset_list)) %>%
          mutate(n_iter = paste0("Iter_", i),
                 row_nm = rownames(.),
                 id = row_number())
      }else{
        score_df_i <- score_meta_data %>%
          dplyr::select(contains("Infect_geneset_123")) %>%
          dplyr::rename_with(.fn = ~names(geneset_list)) %>%
          mutate(n_iter = paste0("Iter_", i),
                 row_nm = rownames(.),
                 id = row_number())
        score_df <- rbind(score_df, score_df_i)
      }
    }
    score_df <- score_df %>%
      dplyr::select(-n_iter,-id) %>%
      group_by(row_nm) %>%
      summarise(across(everything(), mean)) %>%
      ungroup() %>%
      dplyr::slice(match(rownames(meta_data_raw), row_nm)) %>%
      column_to_rownames(var = "row_nm")
    STID_obj@meta.data <- cbind(STID_obj@meta.data, score_df)
    basic_thres_value <- rep(0, len_geneset)
  }else if(score_method == "AUCell"){
    clog_warn("AUCell method is computationally intensive, set n_iter to 1")
    n_iter <- 1
    clog_normal("Build cells rankings..., it may take a while")
    cells_rankings <- AUCell_buildRankings(
      GetAssayData(STID_obj, assay = assay_id, layer = layer_id ),
      splitByBlocks=TRUE,
      plotStats = FALSE,
    )
    clog_normal("Calculate AUC..., it may take a while")
    cells_AUC <- AUCell_calcAUC(
      geneset_list,
      cells_rankings,
      aucMaxRank = ceiling(0.1 * nrow(cells_rankings)),
    )
    score_df <- getAUC(cells_AUC) %>% t() %>% as.data.frame()
    STID_obj@meta.data <- cbind(STID_obj@meta.data, score_df)
    set.seed(seed)
    clog_normal("Explore thresholds..., it may take a while")
    cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=FALSE, assignCells=TRUE)
    basic_thres_value <- sapply(cells_assignment, function(x) x$aucThr$selected) # !!!!

  }else if(score_method %in% c("MeanExp",'SumExp')){
    expr_mat <- GetAssayData(STID_obj, assay = assay_id, layer = layer_id )
    score_df <- lapply(geneset_list, function(gene_set) {
      common_genes <- intersect(gene_set, rownames(expr_mat))
      if (length(common_genes) == 0) {
        clog_warn(paste0("No genes matched the gene set: ", deparse(substitute(gene_set))))
        return(rep(NA, ncol(expr_mat)))
      } else {
        sub_expr_mat <- expr_mat[common_genes, , drop = FALSE]
        if (score_method == "SumExp") {
          return(colSums(sub_expr_mat))
        } else  {
          return(colMeans(sub_expr_mat))
        }
      }
    })
    score_df <- do.call(cbind, score_df) %>%
      as.data.frame() %>%
      `rownames<-`(colnames(expr_mat))
    STID_obj@meta.data <- cbind(STID_obj@meta.data, score_df)
    basic_thres_value <- rep(0, len_geneset)
  }
  names(basic_thres_value) <- names(geneset_list)
  valid_features <- colnames(score_df)

  # >>> Start detect infection spot
  clog_step("Start detect infection spot by specific genesets")
  meta_data <- STID_obj@meta.data
  fwrite(meta_data, file = paste0(output_dir,"Spot_meta.data_(SpotDetect_Geneset).txt"),
         sep = "\t", quote = F,col.names = T,row.names = T)
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  x_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "x_colnm")[[1]]
  y_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "y_colnm")[[1]]
  len_feat <- length(valid_features)
  fwrite(meta_data[,c(samp_colnm,x_colnm,y_colnm,valid_features)],
         file = paste0(output_dir,"Spot_GeneSetScore_(SpotDetect_Geneset).txt"),
         sep = "\t", quote = F,col.names = T,row.names = T)

  #> blur
  if(!is.null(blur_method)){
    match.arg(blur_method, .STID_globals$blur_methods)
    clog_step(paste0("Start spatial smooth for features by mode: ", blur_method))
    for(i in 1:length(loop_single)){
      i_single <- loop_single[i]
      clog_normal(paste0("Processing blur for sample: ", i_single," (",i,"/",length(loop_single),")"))
      i_samp_meta_data <- meta_data[meta_data[[samp_colnm]] == i_single, c(x_colnm,y_colnm,valid_features)]

      for(j in 1:len_feat){
        j_feat <- valid_features[j]
        # clog_normal(paste0("Processing blur for feature: ", j_feat," (",j,"/",len_feat,")"))
        j_feat_exp <- i_samp_meta_data[c("x","y",j_feat)]
        j_feat_exp_smooth <- .SpatialSmooth(meta_data = j_feat_exp, blur_method = blur_method, blur_n = blur_n, blur_sigma = blur_sigma)
        meta_data[meta_data[[samp_colnm]] == i_single, j_feat] <- j_feat_exp_smooth[,j_feat]
      }
    }
    fwrite(meta_data[,c(samp_colnm,x_colnm,y_colnm,valid_features)],
           file = paste0(output_dir,"Spot_GeneSetScore_blur_(SpotDetect_Geneset).txt"),
           sep = "\t", quote = F,col.names = T,row.names = T)
  }

  #> meta2STID
  new_meta_key <- paste0(dir_nm,"_",grp_nm)
  meta2STID <- meta_data[,c(samp_colnm,x_colnm,y_colnm,valid_features),drop = F]
  meta2STID_Label <- meta2STID %>%
    dplyr::select(all_of(valid_features)) %>%
    mutate( across(everything(), ~ NA_real_)) %>%
    rename_with( ~ paste0("Label_", .x), .cols = all_of(valid_features))
  meta2STID_Label <- meta2STID %>%
    bind_cols(meta2STID_Label) %>%
    as.matrix()

  # >>> meta_data_longer
  meta_data_longer <- meta2STID %>%
    mutate(Spot_id = rownames(.),.before = 1) %>%
    pivot_longer(
      cols = all_of(valid_features),
      names_to = "feature",
      values_to = "Expression"
    ) %>%
    group_by(!!sym(samp_colnm), feature) %>%
    mutate(basic_thres_value = basic_thres_value[feature], # !!!! different from the gene
           PosThres_score = PosThres_score,
           PosThres_prob = PosThres_prob) %>%
    mutate(PosThres_prob_value = if(unique(PosThres_prob) >0) quantile(Expression[Expression > basic_thres_value],PosThres_prob) else NA, # !!!
           max_thres_value = if(is.na(unique(PosThres_prob_value))) PosThres_score else ifelse(PosThres_prob_value >= PosThres_score, PosThres_prob_value, PosThres_score),
           max_thres_value = ifelse(max_thres_value>basic_thres_value, max_thres_value, basic_thres_value),
           exp_pos = ifelse(Expression > max_thres_value, Expression, basic_thres_value)) %>%
    ungroup()
  Sap_GeneExp <- meta_data_longer %>%
    group_by(!!sym(samp_colnm), feature) %>%
    summarise(
      all_num = n(),
      all_mean = mean(Expression,na.rm = T),
      all_sd = sd(Expression,na.rm = T),
      all_Q50 = quantile(Expression, 0.50,na.rm = T),
      all_Q75 = quantile(Expression, 0.75,na.rm = T),
      all_Q95 = quantile(Expression, 0.95,na.rm = T),
      all_Q99 = quantile(Expression, 0.99,na.rm = T),
      all_max = max(Expression,na.rm = T),
      basic_thres = unique(basic_thres_value),
      pos_thres = unique(max_thres_value),
      pos_num = sum(Expression > max_thres_value,na.rm = T),
      pos_ratio = pos_num / all_num ,
      pos_mean = mean(Expression[Expression > max_thres_value],na.rm = T),
      pos_sd = sd(Expression[Expression > max_thres_value],na.rm = T),
      pos_Q5 = quantile(Expression[Expression > max_thres_value], 0.05,na.rm = T),
      pos_Q25 = quantile(Expression[Expression > max_thres_value], 0.25,na.rm = T),
      pos_Q50 = quantile(Expression[Expression > max_thres_value], 0.50,na.rm = T),
      pos_Q75 = quantile(Expression[Expression > max_thres_value], 0.75,na.rm = T),
      pos_Q95 = quantile(Expression[Expression > max_thres_value], 0.95,na.rm = T),
    )
  fwrite(Sap_GeneExp, file = paste0(output_dir,"Sap_GeneSetScore_stat_(SpotDetect_Geneset).txt"),
         sep = "\t", quote = F,col.names = T,row.names = F,na = "NA")

  Sap_thres <- meta_data_longer %>%
    dplyr::select(!!sym(samp_colnm), feature, basic_thres_value,PosThres_prob, PosThres_prob_value, PosThres_score, max_thres_value) %>%
    distinct()
  fwrite(Sap_thres, file = paste0(output_dir,"Sap_GeneSetScore_threshold_(SpotDetect_Geneset).txt"),
         sep = "\t", quote = F,col.names = T,row.names = F,na = "NA")

  # >>> Plot_Spatial
  pathogen_org <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "pathogen_org")[[1]]
  host_org <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "host_org")[[1]]
  pathogen_grp <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "pathogen_grp")[[1]]
  subtitle_lable <- paste0(pathogen_grp,"_",pathogen_org," (",host_org,")")

  #>
  process_samp <- function(i){
    i_single <- loop_single[i]
    clog_loop(paste0("Processing sample: ", i_single," (",i,"/",length(loop_single),")"))
    i_output_dir <- paste0(output_dir,i_single,"/")
    i_photo_dir <- paste0(photo_dir,i_single,"/")
    dir.create(i_output_dir,recursive = T, showWarnings =F)
    dir.create(i_photo_dir,recursive = T, showWarnings = F)
    i_meta_data_longer <- meta_data_longer[meta_data_longer[[samp_colnm]] == i_single,]
    i_meta2STID_Label <- meta2STID_Label[meta_data[[samp_colnm]] == i_single, , drop = FALSE]

    #> .Plot_histgram
    p_hist <- .Plot_histgram(plot_data = i_meta_data_longer, samp_colnm = samp_colnm,col = "#99C5E3",
                             title = paste0("GeneSet Score distribution in ", i_single),
                             subtitle = subtitle_lable)
    if(len_feat <= 4){
      width_value <- 4*len_feat
    }else{
      width_value <- 12
    }
    ggsave(p_hist, filename = paste0(i_photo_dir,i_single,"_histgram_(SpotDetect_Geneset).pdf"),
           width = width_value, height = 3.5*ceiling(len_feat/4) + 0.5,limitsize = FALSE)

    #> Plot_Spatial
    if(plot_method == "single"){
      dir.create(paste0(i_photo_dir,"/SpatialPlot_samp/"),recursive = T, showWarnings = F)
      plot_dis_list <- list()
      plot_con_list <- list()
      plot_SpaialPlot_list <- list()

      for(j in 1:len_feat){
        j_feat <- valid_features[j]
        clog_normal(paste0("Processing feature: ", j_feat," (",j,"/",len_feat,")"))
        if(nchar(j_feat)>20){
          j_feat_short <- paste0(str_sub(j_feat,start = 1,end = 20),"(",paste0(sample(letters, 3, replace = TRUE), collapse = ""),")")
          clog_warn(paste0("The length of geneset is greater than 20, use short name '",j_feat_short,"'"))
        }else{
          j_feat_short <- j_feat
        }

        # >
        j_meta_data_longer <- i_meta_data_longer[i_meta_data_longer$feature == j_feat,]
        j_PosThres_prob_value <- j_meta_data_longer$PosThres_prob_value %>% unique()
        j_PosThre_score <- j_meta_data_longer$PosThres_score %>% unique()
        j_basic_thres_value <- j_meta_data_longer$basic_thres_value %>% unique()
        j_max_thres_value <- j_meta_data_longer$max_thres_value %>% unique()
        clog_normal(paste0("basic_thres_value: ", round(j_basic_thres_value,3),
                           "; PosThre_score: ", round(j_PosThre_score,3),
                           "; PosThres_prob_value: ", round(j_PosThres_prob_value,3),
                           "; max_thres_value: ", j_max_thres_value %>% signif(3)))
        subtitle_lable_pos <- paste0(subtitle_lable,"\n","Pos_Thres: ", round(j_max_thres_value,3))

        # >
        j_heat_gene <- j_meta_data_longer$Expression
        title_lable <- paste0(j_feat," (",i_single,")")
        plot_SpaialPlot <- j_meta_data_longer %>%
          mutate(
            Expression = as.numeric(Expression),
            Label = ifelse(Expression > max_thres_value, "pos","neg"),
            Label = factor(Label, levels = c("neg","pos"))
          )
        plot_SpaialPlot_list[[j]] <- plot_SpaialPlot
        fwrite(plot_SpaialPlot, file = paste0(i_output_dir,i_single,"_",j_feat_short,"_PosInfo_(SpotDetect_Geneset).txt"), # !!!
               sep = "\t", quote = F,col.names = T,row.names = F)

        #> coord_Label
        coord_Label <- plot_SpaialPlot[,c("Spot_id","feature")] %>%
          mutate(feature = paste0("Label_", feature)) %>%
          as.matrix()
        i_meta2STID_Label[coord_Label] <- plot_SpaialPlot$Label %>%
          as.character()

        # >
        P_spatial_dis <- Plot_Spatial(plot_data = plot_SpaialPlot, x_colnm = x_colnm, y_colnm = y_colnm, group_by = "Label",
                                      facet_grpnm = NULL, datatype = "discrete",
                                      col = col,pt_size = pt_size, vmin = vmin, vmax = vmax,
                                      title = title_lable, subtitle = subtitle_lable_pos,
                                      black_bg = black_bg)
        P_spatial_con <- Plot_Spatial(plot_data = plot_SpaialPlot, x_colnm = x_colnm, y_colnm = y_colnm, group_by = "Expression",
                                      facet_grpnm = NULL, datatype = "continuous",
                                      col = col,pt_size = pt_size, vmin = vmin, vmax = vmax,
                                      title = title_lable, subtitle = subtitle_lable,
                                      black_bg = black_bg)
        ggsave(P_spatial_dis, filename = paste0(i_photo_dir,"/SpatialPlot_samp/",i_single,"_",j_feat_short,"_dis_(SpotDetect_Geneset).pdf"),
               width = 5, height = 5,limitsize = FALSE)
        ggsave(P_spatial_con, filename = paste0(i_photo_dir,"/SpatialPlot_samp/",i_single,"_",j_feat_short,"_con_(SpotDetect_Geneset).pdf"),
               width = 5, height = 5,limitsize = FALSE)
        plot_dis_list[[j]] <- P_spatial_dis
        plot_con_list[[j]] <- P_spatial_con
      }

      clog_normal(paste0("Save all features plot to a merge pdf"))
      p_dis_all <- wrap_plots(
        plot_dis_list,
        ncol = 4,
      )
      p_con_all <- wrap_plots(
        plot_con_list,
        ncol = 4,
      )
      fwrite(do.call(rbind, plot_SpaialPlot_list),
             file = paste0(i_output_dir,i_single,"_allgene_PosInfo_(SpotDetect_Geneset).txt"),
             sep = "\t", quote = F,col.names = T,row.names = F)

    }else if(plot_method == "merge"){
      plot_SpaialPlot <- i_meta_data_longer %>%
        mutate(
          Expression = as.numeric(Expression),
          Label = ifelse(Expression >max_thres_value, "pos","neg"),
          Label = factor(Label, levels = c("neg","pos"))
        )
      fwrite(plot_SpaialPlot, file = paste0(i_output_dir,i_single,"_allgene_PosInfo_(SpotDetect_Geneset).txt"),
             sep = "\t", quote = F,col.names = T,row.names = F)

      #> coord_Label
      coord_Label <- plot_SpaialPlot[,c("Spot_id","feature")] %>%
        mutate(feature = paste0("Label_", feature)) %>%
        as.matrix()
      i_meta2STID_Label[coord_Label] <- plot_SpaialPlot$Label %>%
        as.character()

      # >
      p_dis_all <- Plot_Spatial(plot_data = plot_SpaialPlot, x_colnm = x_colnm, y_colnm = y_colnm, group_by = "Label",
                                facet_grpnm = "feature", datatype = "discrete",
                                col = col,pt_size = pt_size, vmin = vmin, vmax = vmax,
                                title = paste0("Gene expression in ", i_single), subtitle = subtitle_lable,
                                black_bg = black_bg)
      p_con_all <- Plot_Spatial(plot_data = plot_SpaialPlot, x_colnm = x_colnm, y_colnm = y_colnm, group_by = "Expression",
                                facet_grpnm = "feature", datatype = "continuous",
                                col = col,pt_size = pt_size, vmin = vmin, vmax = vmax,
                                title = paste0("Gene expression in ", i_single), subtitle = subtitle_lable,
                                black_bg = black_bg)
    }
    ggsave(p_dis_all, filename = paste0(i_photo_dir,i_single,"_all_dis_(SpotDetect_Geneset).pdf"),
           width = 20, height = 5*ceiling(len_feat/4),limitsize = FALSE)
    ggsave(p_con_all, filename = paste0(i_photo_dir,i_single,"_all_con_(SpotDetect_Geneset).pdf"),
           width = 20, height = 5*ceiling(len_feat/4),limitsize = FALSE)

    return(i_meta2STID_Label)
  }

  #> parallel or not
  meta2STID_Label_list <- list()
  parallel_workers <- getOption("parallel_workers")
  default_plan <- getOption("default_plan")
  if(parallel_workers > 1){
    plan(default_plan, workers = parallel_workers)
    clog_step(paste0("Start parallel processing sample with ", nbrOfWorkers(), " workers"))
    meta2STID_Label_list <- future_map(seq_along(loop_single), process_samp,
                                      .progress=T, .options = furrr_options(seed = 123))
    names(meta2STID_Label_list) <- loop_single
    plan(sequential)
  }else{
    clog_step("Start processing sample without parallelization")
    plan(sequential)
    for(i in seq_along(loop_single)){
      meta2STID_Label_list[[i]] <- process_samp(i)
    }
  }

  #> AddMetaData
  meta2STID <- lapply(meta2STID_Label_list, as.data.frame) %>%
    bind_rows() %>%
    mutate(across(c(x,y,all_of(valid_features)), as.numeric))
  meta2STID <- meta2STID[rownames(meta_data), , drop = FALSE] # keep order
  STID_obj <- AddMetaData(STID_obj = STID_obj,
                         meta_key = new_meta_key,
                         add_data = meta2STID,
                         dir_nm = dir_nm,
                         grp_nm = grp_nm,
                         asso_key = NULL,
                         description = description)
  STID_obj@meta.data <- meta_data_raw

  # >>> Final
  .save_function_params("SpotDetect_Geneset", envir = environment(), file = paste0(output_dir,"Log_function_params_(SpotDetect_Geneset).log") )
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(SpotDetect_Geneset).log")) %>% invisible()
  return(STID_obj)
}


#' Create spatial plots for STID data
#'
#' This function generates spatial visualizations with discrete or continuous
#' coloring. Can be used directly with an STID object or with a custom data frame.
#'
#' @param STID_obj An STID object (optional if plot_data provided)
#' @param meta_key Character, metadata key to use (default: NULL)
#' @param plot_data Data frame containing plotting data (optional if STID_obj provided)
#' @param x_colnm Character, column name for x coordinates
#' @param y_colnm Character, column name for y coordinates
#' @param group_by Character, column name for grouping/coloring variable
#' @param facet_grpnm Character, column name for faceting (default: NULL)
#' @param datatype Character, data type - "discrete" or "continuous" (default: "discrete")
#' @param col Color palette for visualization (default: COLOR_DIS_CON)
#' @param pt_size Numeric, point size (default: 1.1)
#' @param vmin Numeric or character, minimum value for color scale (default: NULL)
#' @param vmax Numeric or character, maximum value for color scale (default: "p99")
#' @param title Character, plot title (default: NULL)
#' @param subtitle Character, plot subtitle (default: NULL)
#' @param black_bg Logical, whether to use black background (default: FALSE)
#'
#' @return A ggplot object
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create spatial plot from STID object
#' p <- Plot_Spatial(
#'   STID_obj = STID_object,
#'   group_by = "Pathogen_Score",
#'   datatype = "continuous",
#'   black_bg = TRUE
#' )
#'
#' # Create spatial plot from custom data frame
#' p <- Plot_Spatial(
#'   plot_data = my_data,
#'   x_colnm = "x_coord",
#'   y_colnm = "y_coord",
#'   group_by = "cell_type",
#'   datatype = "discrete"
#' )
#' }
Plot_Spatial <- function(STID_obj = NULL,meta_key = NULL,plot_data = NULL,
                         x_colnm = NULL,y_colnm = NULL,group_by = NULL,
                         facet_grpnm = NULL, datatype = "discrete",
                         col = COLOR_DIS_CON, pt_size = 1.1,vmin = NULL, vmax = "p99",
                         title = NULL, subtitle = NULL,black_bg = FALSE){

  # >>> Check input patameter
  # .check_at_least_one_null(STID_obj,plot_data)
  if(!is.null(STID_obj)){
    if (!inherits(STID_obj, "STID")) {
      clog_error("Input object is not an STID object.")
    }
    if(is.null(meta_key)){
      meta_key <- "coord"
    }
    plot_data <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]]
    data_info <- GetInfo(STID_obj, info_key = "data_info")
    if(is.null(x_colnm)){
      x_colnm <- data_info$x_colnm
    }
    if(is.null(y_colnm)){
      y_colnm <- data_info$y_colnm
    }
    if(is.null(group_by)){
      group_by <- data_info$celltype_colnm
    }
    if(is.null(facet_grpnm)){
      facet_grpnm <- data_info$samp_colnm
    }
  }else{
    .check_null_args(plot_data, x_colnm, y_colnm, group_by) # facet_grpnm can be NULL
    plot_data <- as.data.frame(plot_data)
  }

  # >>> Start main pipeline
  if(datatype == "continuous"){
    plot_data[[group_by]] <- as.numeric(plot_data[[group_by]])
    if(!is.null(vmin)){
      if(grepl("^p",vmin)){
        probs <- as.numeric(sub("^p","",vmin))/100
        if(is.na(probs) | probs <0 | probs >1){
          clog_error("vmin with 'pXX' format must be between p0 to p100")
        }
        vmin <- quantile(plot_data[[group_by]], probs = probs, na.rm = TRUE)
      }
      plot_data[[group_by]][plot_data[[group_by]] < vmin] <- vmin
    }
    if(!is.null(vmax)){
      if(grepl("^p",vmax)){
        probs <- as.numeric(sub("^p","",vmax))/100
        if(is.na(probs) | probs <0 | probs >1){
          clog_error("vmax with 'pXX' format must be between p0 to p100")
        }
        vmax <- quantile(plot_data[[group_by]], probs = probs, na.rm = TRUE)
      }
      plot_data[[group_by]][plot_data[[group_by]] > vmax] <- vmax
    }
  }else if(datatype == "discrete"){
    if(is.numeric(plot_data[[group_by]])){
      plot_data[[group_by]] <- plot_data[[group_by]] %>% as.factor()
    }
  }else{
    clog_error("datatype must be 'discrete' or 'continuous'")
  }

  # >
  if(is.list(col)){
    if(!all(c("dis","con") %in% names(col))){
      clog_error("When 'col' is a list, it must contain 'dis' and 'con' elements")
    } else {
      col_dis <- col$dis; col_con <- col$con
    }
  } else if(is.vector(col)){
    col_dis <- col; col_con <- col
  } else {
    clog_error("'col' must be a vector or a list")
  }

  # >
  p1 <- ggplot(plot_data, aes(x = !!sym(x_colnm), y = !!sym(y_colnm), color = !!sym(group_by))) + # x = col, y = row
    theme_void() +
    # facet_wrap(~sample, ncol = 3) +
    geom_point(size = pt_size, shape = 16,stroke = 0) +
    coord_fixed(ratio = 1) +
    scale_y_reverse() +
    labs(x = "Column", y = "Row", color = group_by,
         title = title,subtitle = subtitle) +
    theme(plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))

  if(datatype == "discrete"){
    p1 <- p1 + scale_color_manual(values = c(col_dis), na.value = "grey90") +
      guides(color = guide_legend(override.aes = list(size = 4))) # legend中的点一般都比较小
  }else if(datatype == "continuous"){
    # p1 <- p1 + scale_color_gradient(low = col[1], high = col[2])
    p1 <- p1 + scale_color_gradientn(colours = grDevices::colorRampPalette(col_con)(100),
                                     guide = guide_colorbar(),
                                     na.value = "grey90")

  }

  if(!is.null(facet_grpnm)){
    p1 <- p1 + facet_wrap(as.formula(paste0("~",facet_grpnm)), ncol = 4)
  }

  if(black_bg){
    p1 <- p1 +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "black"), # 注释掉就是白底黑字
            plot.background = element_rect(fill = "black"),
            legend.background = element_rect(fill = "black"),
            strip.background = element_rect(fill=NA,color = NA),
            legend.text = element_text(color = "white"),
            legend.title = element_text(color = "white"),
            plot.title = element_text(color = "white", hjust = 0.5),
            plot.subtitle = element_text(color = "white", hjust = 0.5),
            strip.text = element_text(color = "white", hjust = 0.5)
      )
  }else{
    p1 <- p1 +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_rect(fill=NA,color = NA)
      )
  }
  return(p1)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# internal: SpotDetect_Gene/SpotDetect_Geneset
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Plot histogram of gene expression or gene set scores
#'
#' Internal function to create histogram plots with percentile markers
#'
#' @param plot_data Data frame containing plotting data
#' @param samp_colnm Character, sample column name
#' @param col Character, color for histogram (default: "#99C5E3")
#' @param title Character, plot title (default: NULL)
#' @param subtitle Character, plot subtitle (default: NULL)
#'
#' @return A ggplot object
#'
#' @keywords internal
#'
#' @noRd
.Plot_histgram <- function(plot_data = NULL, samp_colnm = NULL,col = "#99C5E3",
                           title = NULL, subtitle = NULL){

  if(max(plot_data$exp_pos, na.rm = T) <=15){
    if(all(is_integer(plot_data$exp_pos))){
      binwidth = 1; bins = NULL
    }else{
      binwidth = NULL; bins = 15
    }
  }else{
    binwidth = NULL; bins = 15
  }
  vline_df <- plot_data %>%
    group_by(feature) %>%
    summarise(
      x_len = max(exp_pos, na.rm = T) - min(exp_pos, na.rm = T),
      min = min(exp_pos, na.rm = T),
      all_p50 = quantile(exp_pos,0.50,na.rm = T),
      all_p75 = quantile(exp_pos,0.75,na.rm = T),
      all_p95 = quantile(exp_pos,0.95,na.rm = T),
      all_p99 = quantile(exp_pos,0.99,na.rm = T),
      pos_p5 = quantile(exp_pos[exp_pos > basic_thres_value],0.05,na.rm = T) %>% replace_na(0),
      pos_p25 = quantile(exp_pos[exp_pos > basic_thres_value],0.25,na.rm = T) %>% replace_na(0),
      pos_p50 = quantile(exp_pos[exp_pos > basic_thres_value],0.50,na.rm = T) %>% replace_na(0),
      pos_p75 = quantile(exp_pos[exp_pos > basic_thres_value],0.75,na.rm = T) %>% replace_na(0),
      pos_p95 = quantile(exp_pos[exp_pos > basic_thres_value],0.95,na.rm = T) %>% replace_na(0))
  p1 <- ggplot(plot_data, aes(x = exp_pos)) +
    facet_wrap(~feature,ncol = 4,scales = "free") +
    # geom_density(aes(y = ..density..), fill = col, alpha = 0.5) +
    geom_histogram(aes(y = after_stat(count)), binwidth  = binwidth, bins = bins, fill = col, alpha = 0.8) +
    geom_vline(data = vline_df,mapping = aes(xintercept = pos_p5), color = "#FAD510", linetype = "dashed", linewidth = 0.25) +
    geom_vline(data = vline_df,mapping = aes(xintercept = pos_p25), color = "#00A08A", linetype = "dashed", linewidth = 0.25) +
    geom_vline(data = vline_df,mapping = aes(xintercept = pos_p50), color = "#F98400", linetype = "dashed", linewidth = 0.25) +
    geom_vline(data = vline_df,mapping = aes(xintercept = pos_p75), color = "#5BBCD6", linetype = "dashed", linewidth = 0.25) +
    geom_vline(data = vline_df,mapping = aes(xintercept = pos_p95), color = "#FF0000", linetype = "dashed", linewidth = 0.25) +
    geom_text(aes(label = after_stat(count)), stat="bin", binwidth = binwidth, bins = bins,vjust=-0.5, size=2) +
    geom_text(data = vline_df,mapping = aes(x = min+x_len*0.1,label = paste0("pos_p5: ",round(pos_p5,2))), y = 0, color = "#FAD510", vjust = -17*1.2, hjust = 0, size=3) +
    geom_text(data = vline_df,mapping = aes(x = min+x_len*0.1,label = paste0("pos_p25: ",round(pos_p25,2))), y = 0, color = "#00A08A", vjust = -18*1.2, hjust = 0, size=3) +
    geom_text(data = vline_df,mapping = aes(x = min+x_len*0.1,label = paste0("pos_p50: ",round(pos_p50,2))), y = 0, color = "#F98400", vjust = -19*1.2, hjust = 0, size=3) +
    geom_text(data = vline_df,mapping = aes(x = min+x_len*0.1,label = paste0("pos_p75: ",round(pos_p75,2))), y = 0, color = "#5BBCD6", vjust = -20*1.2, hjust = 0, size=3) +
    geom_text(data = vline_df,mapping = aes(x = min+x_len*0.1,label = paste0("pos_p95: ",round(pos_p95,2))), y = 0, color = "#FF0000", vjust = -21*1.2, hjust = 0, size=3) +
    geom_text(data = vline_df,mapping = aes(x = min+x_len*0.1,label = paste0("all_p50: ",round(all_p50,2))), y = 0, color = "black", vjust = -12*1.2, hjust = 0, size=3) +
    geom_text(data = vline_df,mapping = aes(x = min+x_len*0.1,label = paste0("all_p75: ",round(all_p75,2))), y = 0, color = "black", vjust = -13*1.2, hjust = 0, size=3) +
    geom_text(data = vline_df,mapping = aes(x = min+x_len*0.1,label = paste0("all_p95: ",round(all_p95,2))), y = 0, color = "black", vjust = -14*1.2, hjust = 0, size=3) +
    geom_text(data = vline_df,mapping = aes(x = min+x_len*0.1,label = paste0("all_p99: ",round(all_p99,2))), y = 0, color = "black", vjust = -15*1.2, hjust = 0, size=3) +
    theme_test() +
    labs(x = "Expression", y = "Spot number", title = title, subtitle = subtitle) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_rect(fill="grey95"))

  return(p1)
}


#' Apply spatial smoothing to expression data
#'
#' Internal function to smooth spatial expression patterns using image processing
#' techniques
#'
#' @param meta_data Data frame containing spatial coordinates and expression values
#' @param blur_method Character, smoothing method - "isoblur" or "medianblur"
#'        (default: "isoblur")
#' @param blur_n Numeric, number of iterations for median blur (default: 1)
#' @param blur_sigma Numeric, sigma parameter for isoblur (default: 0.8)
#'
#' @return Data frame with smoothed expression values
#'
#' @importFrom imager as.cimg isoblur medianblur
#'
#' @keywords internal
#'
#' @noRd
.SpatialSmooth <- function(meta_data = NULL, blur_method = "isoblur", blur_n = 1, blur_sigma = 0.8){

  col_nm_raw <- colnames(meta_data)
  row_nm_raw <- rownames(meta_data)
  colnames(meta_data) <- c("x","y","raw_score")
  names(row_nm_raw) <- paste0(meta_data$x, "_", meta_data$y)
  meta_data_wider <- meta_data %>%
    pivot_wider(names_from = x, values_from = raw_score) %>%
    column_to_rownames("y") %>%
    as.matrix()
  col_nm <- colnames(meta_data_wider)
  row_nm <- rownames(meta_data_wider)
  if(blur_method == "isoblur"){
    matrix_img <- as.cimg(meta_data_wider)
    matrix_smooth <- isoblur(matrix_img, sigma = blur_sigma,
                             neumann = TRUE, gaussian = TRUE, na.rm = TRUE)
    matrix_smooth <- as.matrix(matrix_smooth)
    arg <- blur_sigma
  }else if(blur_method == "medianblur"){
    matrix_img <- as.cimg(meta_data_wider)
    matrix_smooth <- medianblur(matrix_img, n = blur_n)
    matrix_smooth <- as.matrix(matrix_smooth)
    arg <- blur_n
  }else{
    clog_error("blur_method must be 'isoblur' or 'medianblur'")
  }
  res_df <- matrix_smooth %>%
    as.matrix() %>%
    as.data.frame() %>%
    `colnames<-`(col_nm) %>%
    `rownames<-`(row_nm) %>%
    rownames_to_column("y") %>%
    pivot_longer(-y, names_to = "x", values_to = col_nm_raw[3]) %>%
    mutate(x = as.numeric(x), y = as.numeric(y)) %>%
    mutate(x_y = paste0(x, "_", y), .before = 3) %>%
    filter(x_y %in% names(row_nm_raw)) %>%
    as.data.frame()
  rownames(res_df) <- row_nm_raw[res_df$x_y]
  return(res_df)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GetSapThreshold
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Retrieve sample-specific thresholds from previous analyses
#'
#' This function retrieves stored threshold information from previous
#' SpotDetect analyses (both gene and gene set based)
#'
#' @param STID_obj An STID object
#' @param meta_key Character vector of metadata keys to retrieve thresholds for
#' @param ... Additional arguments passed to methods
#'
#' @return A list containing threshold information for each requested metadata key,
#'         with components 'stat' and 'threshold'
#'
#' @import Seurat
#' @import rlang
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Retrieve thresholds from previous analyses
#' thresholds <- GetSapThreshold(
#'   STID_obj = STID_object,
#'   meta_key = c("M1_SpotDetect_Gene_20240101", "M1_SpotDetect_Geneset_20240101")
#' )
#' }
GetSapThreshold <- function(
    STID_obj,
    meta_key,
    ...
) {
  UseMethod("GetSapThreshold", STID_obj)
}

#' @rdname GetSapThreshold
#' @export
GetSapThreshold.STID <- function(
    STID_obj = NULL,
    meta_key = NULL,
    ...
) {

  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_check( )
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_null_args(meta_key)
  clog_normal(paste0("Your meta_key: ", paste0(meta_key, collapse = ", ")))
  # >>> End check

  # >>> Start main pipeline
  clog_step("Getting GetSapThreshold...")
  meta_info <- GetMetaInfo(STID_obj_before)
  res_list <- list()
  for (i_key in meta_key){
    if(!i_key %in% rownames(meta_info)){
      clog_warn(paste0("meta_key ", i_key, " not found in STID meta_data, skipping..."))
      next
    }
    clog_normal(paste0("Getting meta data for meta_key: ", i_key))
    res_list[[i_key]] <- list()
    i_dir_path <- paste0("./outputdata/",meta_info[i_key,"dir_nm"],"/",meta_info[i_key,"grp_nm"])
    clog_normal(paste0("Getting data from directory: ", i_dir_path))
    i_file_path <- grep(paste0(i_dir_path,"/Sap_"),list.files(i_dir_path, full.names = T), value = T)
    i_stat_path <- i_file_path[grep("_stat_", i_file_path)]
    i_thres_path <- i_file_path[grep("_threshold_", i_file_path)]
    clog_normal(paste0("Stat file name: ", basename(i_stat_path)))
    clog_normal(paste0("Threshold file name: ", basename(i_thres_path)))
    i_stat <- read.table(i_stat_path, header = T, sep = "\t", stringsAsFactors = F)
    i_thres <- read.table(i_thres_path, header = T, sep = "\t", stringsAsFactors = F)
    res_list[[i_key]]$stat <- i_stat
    res_list[[i_key]]$threshold <- i_thres
  }
  clog_end()
  return(res_list)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GetTopGenes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Identify top expressed genes in the dataset
#'
#' This function identifies the top expressed genes either globally or grouped
#' by sample and/or cell type. Useful for finding marker genes or highly
#' expressed pathogen genes.
#'
#' @param STID_obj A Seurat or STID object
#' @param top_n Numeric, number of top genes to return (default: 10)
#' @param pattern Character, regular expression pattern to filter genes (default: NULL)
#' @param grp_by_samp Logical, whether to group by sample (default: FALSE)
#' @param grp_by_celltype Logical, whether to group by cell type (default: FALSE)
#' @param assay_id Character, name of the assay to use (default: "Spatial")
#' @param layer_id Character, name of the layer/data slot to use (default: "counts")
#'
#' @return If no grouping, returns a character vector of top genes.
#'         If grouping, returns a named list of top genes per group.
#'
#' @import Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get top 20 genes globally
#' top_genes <- GetTopGenes(STID_obj = STID_object, top_n = 20)
#'
#' # Get top 10 pathogen genes per sample
#' pathogen_genes <- GetTopGenes(
#'   STID_obj = STID_object,
#'   top_n = 10,
#'   pattern = "^VP",
#'   grp_by_samp = TRUE
#' )
#' }
GetTopGenes <- function(STID_obj, top_n = 10, pattern = NULL,
                        grp_by_samp = FALSE, grp_by_celltype = FALSE,
                        assay_id = "Spatial", layer_id = "counts") {
  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_check( )
  if (!inherits(STID_obj, "Seurat")) {
    clog_error("Input STID_obj is not a Seurat STID_obj.")
  }
  expr_matrix <- GetAssayData(STID_obj, assay = assay_id, layer = layer_id)
  if (ncol(expr_matrix) == 0) {
    clog_error("No expression data found in the specified assay and layer.")
  }
  metadata <- STID_obj@meta.data
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info", sub_key = "samp_colnm")[[1]]
  celltype_colnm <- GetInfo(STID_obj, info_key = "data_info", sub_key = "celltype_colnm")[[1]]
  if (grp_by_samp && !samp_colnm %in% colnames(metadata)) {
    clog_error(paste("Grouping by sample: column '", samp_colnm, "' not found in meta.data."))
  }
  if (grp_by_celltype && !celltype_colnm %in% colnames(metadata)) {
    clog_error(paste("Grouping by celltype: column '", celltype_colnm, "' not found in meta.data."))
  }
  clog_normal(paste0("Your grp_by_samp: ", grp_by_samp))
  clog_normal(paste0("Your grp_by_celltype: ", grp_by_celltype))
  # >>> End check

  # >>> Start main pipeline
  groups <- NULL
  if (grp_by_samp && grp_by_celltype) {
    groups <- interaction(metadata[[samp_colnm]], metadata[[celltype_colnm]],
                          sep = " | ", drop = TRUE) %>%
      as.factor() %>% droplevels()

  } else if (grp_by_samp) {
    groups <- metadata[[samp_colnm]] %>%
      as.factor() %>% droplevels()
  } else if (grp_by_celltype) {
    groups <- metadata[[celltype_colnm]] %>%
      as.factor() %>% droplevels()
  } else {
    avg_expr <- rowMeans(expr_matrix)
    if (!is.null(pattern)) {
      selected_genes <- grep(pattern, names(avg_expr), value = TRUE)
      avg_expr <- avg_expr[selected_genes]
    }
    if (length(avg_expr) == 0) {
      clog_warn("No genes match the pattern.")
      return(character(0))
    }
    top_n <- min(top_n, length(avg_expr))
    top_genes <- names(sort(avg_expr, decreasing = TRUE)[1:top_n])
    clog_end()
    return(top_genes)
  }
  names(groups) <- rownames(metadata)

  # >>> grp_by result
  group_names <- levels(groups)
  result <- vector("list", length(group_names))
  names(result) <- group_names
  for (grp in group_names) {
    clog_normal(paste("Processing group:", grp))
    cells_in_group <- names(groups)[groups == grp]
    if (length(cells_in_group) == 0) {
      result[[grp]] <- character(0)
      next
    }
    expr_sub <- expr_matrix[, cells_in_group, drop = FALSE]
    avg_expr <- rowMeans(expr_sub)
    if (!is.null(pattern)) {
      selected_genes <- grep(pattern, names(avg_expr), value = TRUE)
      avg_expr <- avg_expr[selected_genes]
    }

    if (length(avg_expr) == 0) {
      result[[grp]] <- character(0)
      next
    }
    top_n_adj <- min(top_n, length(avg_expr))
    if (top_n_adj < top_n) {
      clog_warn(paste0("In group '", grp, "': fewer than 'top_n' genes available after filtering."))
    }
    top_genes <- names(sort(avg_expr, decreasing = TRUE)[1:top_n_adj])
    result[[grp]] <- top_genes
  }
  clog_end()
  return(result)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GetGeneStat
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate summary statistics for selected genes
#'
#' This function computes summary statistics (nCount and nFeature) for a set of
#' selected genes, useful for creating pathogen-specific metrics.
#'
#' @param STID_obj A Seurat or STID object
#' @param pattern Character, regular expression pattern to select genes (default: NULL)
#' @param features Character vector of specific gene names to include (default: NULL)
#' @param prefix Character, prefix for output column names
#' @param func Character, function to apply for nCount calculation. Options:
#'        "sum", "mean", "median", "max", "min" (default: "sum")
#' @param assay_id Character, name of the assay to use (default: "Spatial")
#' @param layer_id Character, name of the layer/data slot to use (default: "counts")
#'
#' @return A data frame with two columns:
#'         \item{prefix_nCount(func)}{Summary statistic for total counts}
#'         \item{prefix_nFeature(sum)}{Number of detected features}
#'
#' @import Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate pathogen gene statistics
#' pathogen_stats <- GetGeneStat(
#'   STID_obj = STID_object,
#'   pattern = "^MTB",
#'   prefix = "Pathogen",
#'   func = "sum"
#' )
#'
#' # Add to metadata
#' STID_object <- AddMetaData(STID_object, pathogen_stats)
#' }
GetGeneStat <- function(STID_obj = NULL, pattern = NULL, features = NULL,
                        prefix = NULL, func = "sum",
                        assay_id = "Spatial",layer_id = "counts"
) {

  # >>> Check input patameter
  clog_step("GetGeneStat")
  if (!is(STID_obj, "Seurat")) clog_error("STID_obj must be a Seurat or STID object")
  .check_null_args(STID_obj,prefix)
  if(is.null(pattern) & is.null(features)) clog_error("pattern or features must be provided")

  if(func %in% c("sum","mean","median","max","min")){
    clog_normal(paste0("You choose ",func," to summarize the gene expression, it only available for nCount"))
  }else{
    clog_error("The function must be one of sum,mean,median,max,min")
  }
  # >>> End check

  # >>> Start main pipeline
  selected_genes <- character(0)

  # >>> pattern
  if (!is.null(pattern)) {
    matchingGenes <- grep(pattern = pattern, x = rownames(STID_obj), value = TRUE)
    if (length(matchingGenes) == 0) {
      clog_warn(paste("No genes matched the pattern:", pattern))
    } else {
      selected_genes <- c(selected_genes, matchingGenes)
      clog_normal(paste0("Pattern match: ", length(matchingGenes), " genes (e.g.: ",
                         paste0(head(matchingGenes, 3), collapse = ", "), "...)"))
    }
  }

  # >>> features
  if (!is.null(features)){
    valid_features <- .check_features_exist(STID_obj, features)
    selected_genes <- c(selected_genes, valid_features)
    clog_normal(paste0("Manual features added: ",length(valid_features), " genes (e.g.: ",
                       paste0(head(valid_features, 3), collapse = ", "), " ...)"))
  }

  # >>> merge
  selected_genes <- unique(selected_genes)
  if (length(selected_genes) == 0) {
    clog_error("No genes selected from pattern and feature. Cannot compute statistics.")
  }

  clog_normal(paste0("Total ", length(selected_genes), " genes will be used for ", func, " (e.g.: ",
                     paste0(head(selected_genes, 3), collapse = ", "), " ...)"))
  # >>> calculate statistics
  exp_data <- GetAssayData(STID_obj, assay = assay_id, layer = layer_id)[selected_genes, , drop = FALSE]
  nCount_stat <- switch(func,
                        sum = colSums(exp_data),
                        mean = colMeans(exp_data),
                        median = apply(exp_data, 2, median),
                        max = apply(exp_data, 2, max),
                        min = apply(exp_data, 2, min))
  # STID_obj[[paste0(prefix,"_nCount(",func,")")]] <- nCount_stat
  # STID_obj[[paste0(prefix,"_nFeature(sum)")]] <- colSums(exp_data>0)
  # clog_normal(paste0("Added new column to STID_obj@meta.data: ", prefix))
  res_df <- data.frame(A = nCount_stat,
                       B = colSums(exp_data>0),
                       row.names = names(nCount_stat))
  colnames(res_df) <- c(paste0(prefix,"_nCount(",func,")"),
                        paste0(prefix,"_nFeature(sum)"))
  # >>> End main pipeline

  # >>> return
  return(res_df)
}






