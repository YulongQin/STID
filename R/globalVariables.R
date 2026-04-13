
utils::globalVariables(c(
  # 常用操作符和占位符
  ".", "%:%", "%dopar%", "V<-",

  # 坐标与空间变量
  "x", "y", "x1", "x2", "y1", "y2", "x_y", "x_len",
  "xmin", "xmax", "ymin", "ymax", "xj",
  "All_centernm", "All_ROI_label", "All_Dist2ROIcenter", "All_Dist2ROIedge",
  "ROI_label", "ROI_edge", "ROI_region", "ROI_center", "ROI_count",
  "posROI", "mindist", "minspot",

  # 细胞类型与聚类相关
  "celltype", "cellnum", "cellnum_cut", "celltype_ratio", "label", "label_num",
  "grp", "group", "group_id", "samp_id", "samp_colnm",
  "is_Niche", "Niche_label", "hdr_membership", "probs",
  "color_stRNA1", "scale.factors_spot",

  # 基因表达与分析指标
  "feature", "Expression", "gene", "Count", "Des_short", "Des_short_unique",
  "avg_log2FC", "p_val_adj", "pct.1", "pct.2", "pct_mean", "change",
  "nCount_RNA", "mean_UMI", "spot_num", "all_num", "pos_num",
  "PosThres_prob_value", "max_thres_value", "raw_score",
  "exp_pos", "basic_thres_value", "pos_p5", "pos_p25", "pos_p50", "pos_p75", "pos_p95",
  "all_p50", "all_p75", "all_p95", "all_p99",
  "Correlation", "P-value", "Variable1", "Variable2", "var1", "var2",
  "NES", "new_NES", "aupr_corrected", "test_ligand",
  "f_ij", "N_bin50", "N_bin50_raw", "neighbors", "expected_nicheTye", "entropy_adj",
  "d_s", "d_e", "cellid", "type_counts",

  # 富集分析 (GO/KEGG)
  "GO_ent_up_BP_bar", "GO_ent_up_MF_bar", "GO_ent_up_CC_bar", "KEGG_ent_up_bar",
  "GO_ent_up_BP_bubble", "GO_ent_up_MF_bubble", "GO_ent_up_CC_bubble", "KEGG_ent_up_bubble",
  "GO_ent_down_BP_bar", "GO_ent_down_MF_bar", "GO_ent_down_CC_bar", "KEGG_ent_down_bar",
  "GO_ent_down_BP_bubble", "GO_ent_down_MF_bubble", "GO_ent_down_CC_bubble", "KEGG_ent_down_bubble",
  "GO_ent_updown_BP_bar", "GO_ent_updown_MF_bar", "GO_ent_updown_CC_bar", "KEGG_ent_updown_bar",
  "GO_ent_updown_BP_bubble", "GO_ent_updown_MF_bubble", "GO_ent_updown_CC_bubble", "KEGG_ent_updown_bubble",

  # 网络与互作 (PPI/GRN/CellChat)
  "from", "to", "name", "weight", "size", "node_type", "edge_type",
  "combined_score", "evalue", "qseqid", "sseqid", "p.symmbol",
  "CellChatDB.human", "CellChatDB.mouse",
  "nichenet_weighted_networks", "nichenet_lr_network", "nichenet_ligand_target_matrix",
  "nichenet_gr_network", "nichenet_signaling_network",

  # 聚合与统计指标
  "STID_obj_MS", "STID_obj_before", "is_agg_size", "all_agg_size", "spot_agg_index", "agg_ratio",
  "col_group", "row_group", "spot_agg_index",
  "median_edge", "median_edge_value", "cut_value", "cut_label", "i_interval",
  "Var1", "Freq",

  # Shiny 应用相关
  "currentData", "currentVars", "iris", "mtcars", "p2",

  # Seurat 特定对象/选项
  "Seurat.options", "LayerData",

  # 其他临时或特定变量
  "key", "attribute", "row_nm", "Mg_nm", "is_Mg", "topN", "topNShape",
  "clog_warning", "have_null", "Label", "i_agg_size", "pos_spot", "currentData", "currentVars",
  "org.Hs.eg.db", "org.Mm.eg.db"
))
