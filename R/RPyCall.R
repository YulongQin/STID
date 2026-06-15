

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# conda management
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Configure Python Conda Environment
#'
#' Configures the Python environment used by STID through \pkg{reticulate}.
#' If \code{conda_nm} is not provided, the function tries to use the default
#' Miniconda environment managed by \pkg{reticulate}; if it is unavailable,
#' Miniconda will be installed automatically. If \code{conda_nm} is provided,
#' the function checks whether the specified Conda environment exists and then
#' activates it.
#'
#' @param conda_nm Character, name of an existing Conda environment to use.
#' If \code{NULL}, the default Miniconda path returned by
#' \code{\link[reticulate]{miniconda_path}} will be used
#' (default: \code{NULL}).
#' @param conda_path Character, path to a Conda installation or environment.
#' Currently reserved for future use and not used directly in this function
#' (default: \code{NULL}).
#'
#' @return Invisibly returns \code{NULL}. The function is mainly called for its
#' side effects, including configuring the Python environment and printing
#' Python configuration information.
#'
#' @importFrom reticulate miniconda_path use_condaenv install_miniconda
#' @importFrom reticulate py_config py_available conda_list
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Use the default reticulate Miniconda environment
#' configure_conda()
#'
#' # Use an existing Conda environment
#' configure_conda(conda_nm = "scanpy")
#' }
configure_conda <- function(conda_nm = NULL,conda_path = NULL){
  clog_start()

  #> configure conda
  if(is.null(conda_nm)){
    clog_warn(paste0("No conda_nm, install and use the default miniconda: ",miniconda_path()))
    tryCatch({
      use_condaenv(miniconda_path())
      clog_normal(paste0("miniconda has been installed. Use it directly."))
    }, error = function(e){
      clog_warn("Failed to use the default miniconda, try to install it first.")
      install_miniconda(path = miniconda_path(),
                        update = TRUE,
                        force = FALSE)
      use_condaenv(miniconda_path(),required = TRUE)
      py_config()
    })
  }else{
    envs <- conda_list()
    if (conda_nm %in% envs$name) {
      clog_normal(paste0("detected existing conda environment '", conda_nm, "'. Use it directly."))
      use_condaenv(conda_nm,required = TRUE)
      py_config()
    }else{
      clog_error(paste0("The Conda environment '", conda_nm, "' is not found. Please use conda_list() to check existing environments and choose a existing name"))
    }
  }

  #> configure python
  if(py_available()){
    clog_normal("Sucessfully configured python, the following is your py_config: \n")
    print(py_config())
  }else{
    clog_error("Failed to configure python, please check your conda environment and python installation.")
  }
  clog_end()
}


#' Check Required Python Packages
#'
#' Checks whether the required Python packages are installed and importable in
#' the currently configured Python environment. This function is intended for
#' internal use after \code{\link{configure_conda}} has configured Python
#' through \pkg{reticulate}.
#'
#' @param pkgs Character vector, names of Python packages to check.
#'
#' @return Invisibly returns \code{NULL}. The function prints package checking
#' messages and warnings through STID logging functions.
#'
#' @importFrom reticulate py_available py_config py_list_packages py_module_available
#'
#' @keywords internal
#' @noRd
.check_py_pkgs <- function(pkgs = NULL){
  if(!py_available()){
    clog_error("Python is not configured, please run configure_conda() first.")
  }else{
    clog_normal("Python is configured, the following is your py_config: \n")
    print(py_config())
  }

  #>
  clog_normal("Check if the required python packages are installed in the current conda environment.")
  all_installed_pkgs <- py_list_packages()$package
  not_installed_pkgs <- setdiff(pkgs, all_installed_pkgs)
  if(length(not_installed_pkgs) > 0){
    clog_warn(paste0("The following python packages are not installed: ", paste(not_installed_pkgs, collapse = ", ")))
    clog_normal("Try to install the missing packages in the current conda environment by conda_install()")
  }else{
    clog_normal("All required python packages are installed.")
  }

  #>
  clog_normal("Check if the required python packages can be imported successfully.")
  not_available_pkgs <- c()
  for(pkg in pkgs){
    if(!py_module_available(pkg)){
      not_available_pkgs <- c(not_available_pkgs, pkg)
    }
  }
  if(length(not_available_pkgs) > 0){
    clog_warn(paste0("The following python packages can not be imported: ", paste(not_available_pkgs, collapse = ", ")))
    clog_normal("Use import() to check the import error.")
  }else{
    clog_normal("All required python packages can be imported successfully.")
  }
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Convert between RDS and H5AD files
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Convert AnnData h5ad file to Seurat rds file
#'
#' Converts a Python AnnData h5ad file to a Seurat rds file, preserving spatial
#' coordinates, dimensional reductions, and optionally SCT-transformed data.
#'
#' @param file_path Character, path to the input h5ad file
#' @param data_type Character, data type - "stRNA" (spatial) or "scRNA" (single-cell)
#'        (default: "stRNA")
#' @param convert_mode Character, conversion method - "scanpy" or "seurat"
#'        (default: "scanpy")
#' @param assay_id Character, assay name for the Seurat object (default: NULL,
#'        auto-detected based on data_type)
#' @param X_index Character, which matrix to use as counts - "X" or "rawX"
#'        (default: "rawX")
#' @param binsize Numeric, bin size for spatial coordinate scaling (default: 1)
#' @param SCT_index Logical, whether to convert SCT results from stereopy
#'        (default: FALSE)
#' @param reduction_index Logical, whether to preserve dimensional reductions
#'        (default: FALSE)
#' @param image_index Logical, whether to create spatial image object
#'        (default: NULL, auto-detected based on data_type)
#' @param return_object Logical, whether to return the Seurat object
#'        (default: TRUE)
#' @param grp_nm Character, group name for output organization (default: "sample1")
#' @param dir_nm Character, directory name for output (default: "M0_h5ad2rds")
#'
#' @return If return_object = TRUE, returns a Seurat object; otherwise returns NULL
#'
#' @import reticulate
#' @import Matrix
#' @import rjson
#' @import Seurat
#' @importFrom anndata read_h5ad
#' @importFrom SeuratDisk Convert LoadH5Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert spatial transcriptomics h5ad to Seurat
#' configure_conda()
#' seurat_obj <- h5ad2rds(
#'   file_path = "data/spatial_data.h5ad",
#'   data_type = "stRNA",
#'   binsize = 100,
#'   reduction_index = TRUE,
#'   image_index = TRUE
#' )
#' }
h5ad2rds <- function(
    file_path = NULL,
    data_type = "stRNA", convert_mode = "scanpy",
    assay_id = NULL, X_index = "rawX", binsize = 1,
    SCT_index = FALSE, reduction_index = FALSE, image_index = NULL, return_object = TRUE,
    grp_nm = NULL,dir_nm = "M0_h5ad2rds"
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check()
  match.arg(data_type, choices = .STID_globals$supported_types)
  match.arg(convert_mode, choices = c("scanpy","seurat"))
  match.arg(X_index, choices = c("X","rawX"))
  .check_py_pkgs(pkgs = c("scanpy","squidpy","pandas","numpy","seaborn","matplotlib"))
  if(data_type == "stRNA"){
    assay_id <- ifelse(is.null(assay_id),"Spatial",assay_id)
    image_index <- ifelse(is.null(image_index),TRUE,image_index)
  }else if(data_type == "scRNA"){
    assay_id <- ifelse(is.null(assay_id),"RNA",assay_id)
    image_index <- FALSE
  }
  clog_normal(paste0("Your binsize is ",binsize))
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir


  # >>> Start main pipeline
  clog_step("Read the h5ad file")
  clog_normal("Processing the h5ad")
  file_nm <- gsub(".h5ad", ".rds", basename(file_path))
  clog_normal(paste0("Converting ", basename(file_path)," -> ", file_nm))

  #
  if(convert_mode == "seurat"){
    clog_normal("convert h5ad to rds by seurat mode")
    Convert(file_path, dest = "h5seurat", assay = assay_id, overwrite = TRUE)

    h5file <- paste(paste(unlist(strsplit(file_path, "h5ad", fixed = TRUE)), collapse='h5ad'), "h5seurat", sep="")
    clog_normal(paste(c("Finished! Converting h5ad to h5seurat file at:", h5file), sep=" ", collapse=NULL))

    object <- LoadH5Seurat(h5file,meta.data = F, misc = F) # 不读入misc中的内容，默认只有count数据
    clog_normal(paste(c("Successfully load h5seurat:", h5file), sep=" ", collapse=NULL))

    if(image_index){
      if(!("x" %in% colnames(object@meta.data)) & !("y" %in% colnames(object@meta.data)) &
         !is.null(object@reductions[[assay_id]])){
        clog_normal("Add x and y to meta.data")
        object@meta.data[c("x","y")] <- object@reductions[[assay_id]]@cell.embeddings # seurat格式meta中的x和y一般来源于这里
      }
    }

  }else if(convert_mode == "scanpy"){
    clog_normal("Converting h5ad to rds by scanpy mode")
    ad <- read_h5ad(file_path)
    clog_normal('if an invalid class dgRMatrix object: x slot is not of type double error occurs,
          please use adata.X = adata.raw.X.astype("float32") in python')

    if(is.null(ad$raw) | X_index =="X") {
      clog_normal(paste0("your X_index is ",X_index))
      clog_warn("No raw data, use X as counts")
      X_index <- "X"
      counts <- ad$X
      X_type <- .determine_X_type(counts)
      clog_normal(paste("X matrix is ", X_type))
      if(X_type == "scale.data"){
        clog_normal("Convert scale.data to dgRMatrix")
        counts <- as(counts, "dgRMatrix")
      }
    }else if(X_index == "rawX"){
      counts <- ad$raw$X
      X_type <- .determine_X_type(counts)
      clog_normal(paste("rawX matrix is ", X_type))
      if(X_type == "scale.data"){
        clog_normal("Convert scale.data to dgRMatrix")
        counts <- as(counts, "dgRMatrix")
      }
    }

    ad_dim <- dim(counts)
    metadata <- ad$obs
    if(X_index == "X"){
      features <- ad$var
    }else if(X_index == "rawX"){
      features <- ad$raw$var
    }else{
      clog_error("X_index is not correct")
    }

    if(image_index){
      if(!("x" %in% metadata) & !("y" %in% metadata)){
        clog_normal("Add x and y to meta.data")
        coord <- ad$obsm[[str_to_lower(assay_id)]] %>%
          as.data.frame() %>%
          mutate(across(everything(),as.numeric))
        colnames(coord) <- c("x","y")
        metadata <- cbind(metadata,coord)
      }
    }
    # rm(ad);gc()

    counts_dgT <- Matrix::summary(counts)
    rm(counts);gc()
    counts_dgC <- sparseMatrix(i = counts_dgT$i, j = counts_dgT$j,  # 重构dgRMatrix -> dgCMatrix
                               x = counts_dgT$x, dims = ad_dim, repr = "C") %>%

      t()
    rm(counts_dgT);gc()
    dimnames(counts_dgC) <- list(rownames(features), rownames(metadata))
    object <- CreateSeuratObject(counts = counts_dgC,
                                 meta.data = metadata,
                                 assay = assay_id,
                                 names.field = 1,
                                 names.delim = "_")
    rownames(object) <- rownames(features)
    rm(counts_dgC);gc()
  }
  print(paste0("The gene name is : ",paste0(head(rownames(object)),collapse = " ")," ..."))

  #
  clog_step("Add reduction")
  if( reduction_index){
    if(convert_mode == "scanpy"){ # Seurat格式会自动转
      if(!is.null(ad$obsm[["X_umap"]])){
        X_umap <- ad$obsm[["X_umap"]]
        colnames(X_umap) <- paste0("umap_",1:ncol(X_umap))
        rownames(X_umap) <- rownames(ad)
        object@reductions$umap <- CreateDimReducObject(embeddings = X_umap, key = "umap_",assay = assay_id)
      }
      if(!is.null(ad$obsm[["X_tsne"]])){
        X_tsne <- ad$obsm[["X_tsne"]]
        colnames(X_tsne) <- paste0("tsne_",1:ncol(X_tsne))
        rownames(X_tsne) <- rownames(ad)
        object@reductions$tsne <- CreateDimReducObject(embeddings = X_tsne, key = "tsne_",assay = assay_id)
      }
      if(!is.null(ad$obsm[["X_pca"]])){
        X_pca <- ad$obsm[["X_pca"]]
        colnames(X_pca) <- paste0("pca_",1:ncol(X_pca))
        rownames(X_pca) <- rownames(ad)
        object@reductions$pca <- CreateDimReducObject(embeddings = X_pca, key = "pca_",assay = assay_id)
      }
    }
  }else{
    object@reductions[[assay_id]] <- NULL
  }

  #
  clog_step("Add misc")
  if(SCT_index){
    if (
      !is.null(object@misc$sct_counts) &&
      !is.null(object@misc$sct_data) &&
      !is.null(object@misc$sct_scale) &&
      !is.null(object@misc$sct_cellname) &&
      !is.null(object@misc$sct_genename) &&
      !is.null(object@misc$sct_top_features)
    ) {
      clog_normal("convert stereopy SCT result to seurat SCT result")
      sct.assay.out <- CreateAssayObject(counts=object[['Spatial']]@counts, check.matrix=FALSE)
      sct.assay.out <- SetAssayData(
        object = sct.assay.out,
        slot = "data",
        new.data = log1p(x=GetAssayData(object=sct.assay.out, slot="counts"))
      )
      sct.assay.out@scale.data <- as.matrix(object@misc$sct_scale)
      colnames(sct.assay.out@scale.data) <- object@misc$sct_cellname
      rownames(sct.assay.out@scale.data) <- object@misc$sct_scale_genename
      sct.assay.out <- Seurat:::SCTAssay(sct.assay.out, assay.orig='Spatial')
      Seurat::VariableFeatures(object = sct.assay.out) <- object@misc$sct_top_features
      object[['SCT']] <- sct.assay.out
      DefaultAssay(object=object) <- 'SCT'

      # TODO: tag the reductions as SCT, this will influence the find_cluster choice of data
      object@reductions$pca@assay.used <- 'SCT'
      object@reductions$umap@assay.used <- 'SCT'
      assay.used <- 'SCT'
      clog_normal("Got SCTransform result in object, create a new SCTAssay and set it as default assay.")
    }
  } else {
    assay.used <- assay_id
    clog_normal("Get raw counts only, auto create log-normalize data.")
    object <- NormalizeData(object, assay = assay.used)
    clog_normal("Create log-normalize data.")
  }

  #
  clog_step("Add image")
  if( image_index){
    clog_normal("Start add image to seurat object, This may take some minutes...")
    cell_coords <- unique(object@meta.data[, c('x', 'y')]) # 一般不会重复吧
    cell_coords['cell'] <- row.names(cell_coords)

    # meta中的坐标在seurat中是修改的，但是metadata中的是没有修改的
    clog_normal("Here we do the same three steps for the x and y coordinates:
                  1. Subtract the minimum value 2. Divide by binsize 3. Plus one")
    cell_coords$x <- (cell_coords$x - min(cell_coords$x))/binsize + 1 # 从1开始，rds的坐标都是从1开始的？
    cell_coords$y <- (cell_coords$y - min(cell_coords$y))/binsize + 1

    # object of images$slice1@image, all illustrated as 1 since no concrete pic
    tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))

    # object of images$slice1@coordinates, concrete coordinate of X and Y
    tissue_positions_list <- data.frame(row.names = cell_coords$cell,
                                        tissue = 1,
                                        row = cell_coords$y, col = cell_coords$x,
                                        imagerow = cell_coords$y, imagecol = cell_coords$x)
    # @images$slice1@scale.factors
    scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                     tissue_hires_scalef = 1,
                                     tissue_lowres_scalef = 1))

    # generate object @images$slice1
    .generate_stereo_spatial <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
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

    stereo_spatial <- .generate_stereo_spatial(image = tissue_lowres_image,
                                              scale.factors = fromJSON(scalefactors_json), # 啥作用？
                                              tissue.positions = tissue_positions_list)

    # can be thought of as a background of spatial
    # import image into seurat object
    object@images[['slice1']] <- stereo_spatial
    object@images$slice1@key <- "slice1_"
    object@images$slice1@assay <- assay.used

  }
  clog_normal("The seurat object info:")
  print(object)

  #
  clog_step("Save rds")
  outfile <- paste0(output_dir, "/", file_nm)
  clog_normal(paste0("Save rds to ", outfile))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  gc()
  saveRDS(object, outfile)
  # >>> End main pipeline

  # >>> Final
  .save_function_params("h5ad2rds", envir = environment(), file = paste0(output_dir,"Log_function_params_(h5ad2rds).log") )
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(h5ad2rds).log")) %>% invisible()
  if(return_object){
    return(object)
  }else{
    return( invisible(NULL))
  }
}

#' Determine the type of expression matrix from AnnData
#'
#' Internal function to classify the expression matrix type based on matrix class
#' and properties.
#'
#' @param counts Matrix object from AnnData (typically dgRMatrix)
#'
#' @return Character string indicating matrix type: "count", "data", "scale.data",
#'         or "unknown"
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' X_type <- determine_X_type(adata$X)
#' }
.determine_X_type <- function(counts = NULL){
  if(!inherits(counts, "dgRMatrix")){
    return("scale.data")
  }else{
    if(all(counts[1,] %% 1 == 0) & ((colSums(counts) %>% unique() %>% length())>1)){
      return("count")
    }else if((colSums(counts) %>% unique() %>% length()) == 1){
      return("data")
    }else{
      return("unknown")
    }
  }
}


#' Convert Seurat rds file to AnnData h5ad file
#'
#' Converts a Seurat rds file to a Python AnnData h5ad file for use with scanpy
#' and other Python-based tools.
#'
#' @param seurat_obj A Seurat object to convert
#' @param data_type Character, data type - "stRNA" (spatial) or "scRNA" (single-cell)
#'        (default: "stRNA")
#' @param convert_mode Character, conversion method - "seurat" (using SeuratDisk)
#'        (default: "seurat")
#' @param assay_id Character, assay name to convert (default: NULL, auto-detected)
#' @param grp_nm Character, group name for output organization (default: NULL)
#' @param dir_nm Character, directory name for output (default: "M0_rds2h5ad")
#'
#' @return NULL (invisible), saves h5ad file to disk
#'
#' @import reticulate
#' @import Matrix
#' @import dplyr
#' @import rjson
#' @import Seurat
#' @import ggplot2
#' @importFrom SeuratDisk SaveH5Seurat Convert LoadH5Seurat
#' @import SeuratObject
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert Seurat object to h5ad
#' configure_conda()
#' rds2h5ad(
#'   seurat_obj = seurat_object,
#'   data_type = "stRNA",
#'   grp_nm = "sample1"
#' )
#' }
rds2h5ad <- function(
    seurat_obj = NULL, data_type = "stRNA",
    convert_mode = "seurat", assay_id = NULL,
    grp_nm = NULL,dir_nm = "M0_rds2h5ad"
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check()
  match.arg(data_type, choices = .STID_globals$supported_types)
  match.arg(convert_mode, choices = c("scanpy","seurat"))
  .check_py_pkgs(pkgs = c("scanpy","squidpy","pandas","numpy","seaborn","matplotlib"))
  if(data_type == "stRNA"){
    assay_id <- "Spatial"
  }else if(data_type == "scRNA"){
    assay_id <- "RNA"
  }
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir

  # >>> Start main pipeline
  .rds2h5ad_pipeline(
    seurat_obj = seurat_obj,
    convert_mode = convert_mode, assay_id = assay_id,
    output_dir = output_dir
  )
  # >>> End main pipeline

  # >>> Final
  .save_function_params("rds2h5ad", envir = environment(), file = paste0(output_dir,"Log_function_params_(rds2h5ad).log") )
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(rds2h5ad).log")) %>% invisible()
  invisible(NULL)
}


.rds2h5ad_pipeline <- function(
    seurat_obj = NULL,
    convert_mode = "seurat", assay_id = NULL,
    output_dir = NULL
){
  clog_normal("Read the seurat object")
  clog_normal("Processing the seurat object")
  file_nm <- as.character(substitute(seurat_obj))
  clog_normal(paste0("Converting ", file_nm, ".rds -> ", file_nm, ".h5ad"))

  #
  if(convert_mode == "seurat"){
    clog_normal("convert rds to h5ad by seurat mode")
    seurat_obj <- suppressMessages({
      UpdateSeuratObject(seurat_obj)
    })
    seurat_obj@meta.data <- seurat_obj@meta.data %>%
      mutate(across(where(is.factor), as.character)) # 转为character，否则h5ad中只保留数字格式
    seurat_obj[[assay_id]] <- as(seurat_obj[[assay_id]], Class = "Assay")
    seurat_obj@assays[[assay_id]]["scale.data"] <- NULL
    clog_warn("Remove the scale.data")

    clog_normal("SaveH5Seurat")
    suppressMessages({
      SaveH5Seurat(seurat_obj,filename="test.h5seurat", overwrite = TRUE)
    })
    SaveH5Seurat(seurat_obj,filename="test.h5seurat", overwrite = TRUE)
    clog_normal("Convert to h5ad")
    outfile <- paste0(output_dir, "/", file_nm,".h5ad")
    clog_normal(paste0("Save h5ad to ", outfile))
    Convert("test.h5seurat",
            dest = outfile,
            overwrite = TRUE)
    file.remove("test.h5seurat")

  }else if(convert_mode == "scanpy"){
    clog_normal("Converting rds to h5ad by scanpy mode")
    clog_error("The scanpy mode of rds2h5ad is not implemented yet, please use seurat mode for now.")
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# .RPyCall_squidpy
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Internal pipeline for Squidpy spatial co-localization analysis
#'
#' Converts selected STID single-sample data to h5ad format, calls Python
#' Squidpy through reticulate, computes spatial neighborhood enrichment and
#' co-occurrence patterns, and saves spatial and co-occurrence plots for each
#' sample.
#'
#' @param STID_obj STID object containing spatial transcriptomics data and metadata
#' @param loop_single Character vector, sample identifiers to process
#' @param niche_key Character, niche key used to restrict the analysis to niche
#'   cells. If NULL, cells are selected from the metadata specified by meta_key
#' @param meta_key Character, metadata key used when niche_key is NULL
#'   (default: "coord")
#' @param group_by Character, column name used as the cell group or cluster label
#'   for Squidpy neighborhood enrichment and co-occurrence analysis
#' @param group_use Character vector, specific groups in group_by to include
#'   in the analysis. If NULL, all groups are used
#' @param coord_interval Numeric vector, spatial coordinate intervals for each
#'   sample, used as spot size and neighborhood radius
#' @param tmp_dir Character, temporary directory. Currently reserved for interface
#'   consistency with other internal pipelines
#' @param output_dir Character, directory for intermediate h5ad files and output
#'   data
#' @param photo_dir Character, directory for saving spatial and co-occurrence
#'   plots
#'
#' @return List of Squidpy neighborhood enrichment results per sample
#'
#' @import reticulate
#' @importFrom anndata read_h5ad
#'
#' @keywords internal
#' @noRd
.RPyCall_squidpy <- function(
    STID_obj = NULL,
    loop_single = NULL,
    niche_key = NULL,
    meta_key = NULL,
    group_by = NULL,
    group_use = NULL,
    coord_interval = NULL,
    tmp_dir = NULL,
    output_dir = NULL,
    photo_dir = NULL
){
  # Reference: https://kayla-morrell.github.io/SquidpyR/articles/C_reproduceSquidpyFigures.html
  results_list <- list()
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  data_format <- GetInfo(STID_obj, info_key = "data_info", sub_key = "data_format")[[1]]
  data_platform <- GetInfo(STID_obj, info_key = "data_info", sub_key = "data_platform")[[1]]

  #> check
  .check_py_pkgs(pkgs = c("scanpy","squidpy","pandas","numpy","seaborn","matplotlib"))
  if(!is.null(niche_key)){
    clog_normal("Using niche_key ...")
    Niche_cells <- GetSSNicheCells(STID_obj = STID_obj, niche_key = niche_key, loop_id = loop_single) %>%
      bind_rows()
    .check_column_exist(Niche_cells, "is_Niche")
    Niche_cells <- Niche_cells %>%
      filter(is_Niche) # only niche, not all
  }else{
    if(is.null(meta_key)){
      clog_warn("Both niche_key and meta_key are NULL, will use meta_key: coord ")
      meta_key <- "coord"
    }
    clog_normal("Using meta_key ...")
    Niche_cells <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]] %>%
      filter(!!sym(samp_colnm) %in% loop_single)
  }

  #>
  seurat_obj <- as.Seurat(STID_obj)
  # seurat_obj <- suppressMessages({
  #   UpdateSeuratObject(seurat_obj)
  # })
  seurat_obj <- subset(seurat_obj, cells = rownames(Niche_cells))
  clog_step("Start RDS to H5AD conversion")
  .rds2h5ad_pipeline(
    seurat_obj = seurat_obj,
    convert_mode = "seurat", assay_id = "Spatial",
    output_dir = output_dir
  )

  #>
  clog_step("Start RPyCall_squidpy analysis")
  clog_normal(paste0("Import python packages"))
  sc <- reticulate::import("scanpy")
  sq <- reticulate::import("squidpy")
  plt <- reticulate::import("matplotlib.pyplot")
  sns <- reticulate::import("seaborn")
  pd <- reticulate::import("pandas")
  np <- reticulate::import("numpy")
  sns$set_theme(style = "white")  # 设置白色背景主题
  # sc$logging$print_header()
  sc$set_figure_params(facecolor="white", figsize = c(8,8)) #
  sc$settings$verbosity <- 'hint'
  # plt$rcParams['pdf.fonttype'] = 42
  # plt$rcParams['ps.fonttype'] = 42
  # plt$rcParams['svg.fonttype'] = 'none'
  # plt$rcParams['figure.dpi'] = 900
  # plt$rcParams['savefig.dpi'] = 900

  py_run_string("
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.dpi'] = 900
plt.rcParams['savefig.dpi'] = 900
")
  # py_run_string("print(plt.rcParams['pdf.fonttype'])")

  #
  clog_normal(paste0("read and process h5ad data"))
  py_run_string("import scanpy as sc; import numpy as np")
  py_run_string(paste0("adata = sc.read_h5ad('",output_dir,"/seurat_obj.h5ad')"))
  py_run_string("adata.obsm['spatial'] = adata.obs[['x','y']].to_numpy(dtype=np.uint32)")
  py_run_string("adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})")
  py_run_string(paste0("adata.write_h5ad('",output_dir,"/adata.h5ad')"))

  #>
  # adata <- sq$datasets$seqfish() # demo
  adata <- sc$read_h5ad(paste0(output_dir,"/adata.h5ad")) # class(adata)
  for(i in seq_along(loop_single)){
    # i = 1
    i_single <- loop_single[i]
    clog_loop(paste0("Processing samp_id: ", i_single, " (", i, "/", length(loop_single), ")"))
    photo_subdir <- paste0(photo_dir,"/",i_single)
    dir.create(paste0(photo_subdir,"/co_occurrence"), recursive = TRUE,showWarnings = FALSE)
    sc$settings$figdir <- photo_subdir
    i_adata <- adata[adata$obs[samp_colnm] == i_single, ]
    if(!is.null(group_use)){
      i_adata <- i_adata[i_adata$obs[group_by] %in% group_use, ]
    }

    #> spatial plot
    sc$pl$spatial(
      i_adata,
      color = group_by,
      size = 1.15,
      img_key = "hires",
      show = FALSE,
      # save = "_spatial_plot.pdf",
      palette = "gist_ncar",
      spot_size = coord_interval[i],
      scale_factor = 1,
      title = paste0(i_single, " spatial plot")
    )
    plt$savefig(
      paste0(photo_subdir, "/spatial_plot_group.pdf"),
      bbox_inches = "tight",
      dpi = 900L
    )
    plt$close()

    #> nhood_enrichment
    clog_normal("Start nhood_enrichment")
    if(data_format == "square_grid") {
      sq$gr$spatial_neighbors(i_adata, coord_type="generic",radius=coord_interval[i]) # added to adata.obsp, adata.uns
      # sq$gr$spatial_neighbors(i_adata, n_rings=2, coord_type="grid", n_neighs=4)
    }else if(data_format == "hex_grid" & data_platform == "Visium") {
      sq$gr$spatial_neighbors(i_adata, n_rings=2, coord_type="grid", n_neighs=6) # added to adata.obsp, adata.uns
    }
    sq$gr$nhood_enrichment(i_adata, cluster_key = group_by,n_perms = 100) # added to adata.uns
    arr <- np$array(i_adata$uns[[paste0(group_by,"_nhood_enrichment")]][["zscore"]])
    results_list[[i_single]] <- arr %>% as.data.frame()
    vmax_value <- quantile(arr %>% as.vector(), probs = 0.9,na.rm = TRUE)
    # vmin_value <- quantile(arr %>% as.vector(), probs = 0.1,na.rm = TRUE)
    sq$pl$nhood_enrichment(
      i_adata,
      cluster_key = group_by,
      mode = "zscore",
      vmin = 0, vmax = vmax_value,
      # save = "nhood_enrichment.pdf",
      cmap = 'Blues'
    )
    plt$savefig(
      paste0(photo_subdir, "/nhood_enrichment.pdf"),
      bbox_inches = "tight",
      dpi = 900L
    )
    plt$close()

    #> co_occurrence
    clog_normal("Start co_occurrence")
    i_adata$obs[,group_by] <- make.names(i_adata$obs[,group_by]) %>%
      factor()
    loop_grps <- unique(i_adata$obs[,group_by]) %>% as.character() %>%  sort()
    sq$gr$co_occurrence(
      i_adata,
      cluster_key = group_by
    )
    for(j in seq_along(loop_grps)){
      # j = 1
      j_grp <- loop_grps[j]
      sq$pl$co_occurrence(
        i_adata,
        cluster_key = group_by,
        clusters = j_grp,
        # save = paste0("co_occurrence/", j_grp, "_co_occurrence.pdf"),
        dpi = 900
      )
      plt$savefig(
        paste0(photo_subdir, "/co_occurrence/", j_grp, "_co_occurrence.pdf"),
        bbox_inches = "tight",
        dpi = 900L
      )
      plt$close()
    }
  }
  return(results_list)
}






