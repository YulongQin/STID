#' @include zzz.R
#' @include Utilities.R
#' @importFrom methods setClass
#' @importFrom methods setMethod
#' @importFrom methods new
#' @importFrom methods getClass
#' @importFrom methods show
NULL


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SingleSampNiche
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The SingleSampNiche Class
#'
#' A container for niche analysis results from a single sample. Stores information
#' about detected niches, including their spatial coordinates, associated cells,
#' and gene expression patterns.
#'
#' @slot samp_info List containing sample information:
#'   \itemize{
#'     \item samp_id: Character, sample identifier
#'   }
#' @slot niche_info Data frame containing metadata for each niche:
#'   \itemize{
#'     \item niche_key: Character, unique identifier for the niche
#'     \item meta_key: Character, reference to metadata source
#'     \item ROI_type: Character, type of ROI ("ROI", "Spot", or "Unknown")
#'     \item pos_colnm: Character, column name for positive spot labels
#'     \item center_colnm: Character, column name for center spot labels
#'     \item edge_colnm: Character, column name for edge spot labels
#'     \item all_label_colnm: Character, column name for all spot labels
#'     \item all_dist_colnm: Character, column name for distance to center
#'     \item niche_cells_num: Integer, number of cells in the niche
#'     \item niche_genes_num: Integer, number of genes analyzed
#'     \item description: Character, description of the niche
#'   }
#' @slot niche_cells List containing data frames of cell-level information for
#'   each niche key
#' @slot niche_genes List containing data frames of gene-level information for
#'   each niche key
#'
#' @name SingleSampNiche-class
#' @rdname SingleSampNiche-class
#' @exportClass SingleSampNiche
#'
#' @family niche-classes
#'
#' @aliases SingleSampNiche
#'
SingleSampNiche <- setClass(
  Class = "SingleSampNiche",
  slots = list(
    samp_info = "list",
    niche_info = "data.frame",
    niche_cells = "list",
    niche_genes = "list"
  ),
  prototype = prototype(
    samp_info = list(
      samp_id = character()
    ),
    niche_info = data.frame(
      niche_key = character(),
      meta_key = character(),
      ROI_type = character(), # "ROI","Spot","Unknown"
      pos_colnm = character(),
      center_colnm = character(),
      edge_colnm = character(),
      all_label_colnm = character(),
      all_dist_colnm = character(),
      niche_cells_num = integer(),
      niche_genes_num = integer(),
      description = character(),
      stringsAsFactors = FALSE,
      row.names = NULL
    ),
    niche_cells = list(),
    niche_genes = list()
  )
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MultiSampNiche
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The MultiSampNiche Class
#'
#' A container for comparative niche analysis across multiple samples. Stores
#' integrated information from multiple SingleSampNiche objects for cross-sample
#' comparisons.
#'
#' @slot samp_info List containing multi-sample information:
#'   \itemize{
#'     \item multi_id: Character, unique identifier for the multi-sample analysis
#'     \item samp_id: Character vector, sample identifiers included
#'     \item samp_grp: Character, sample group assignments
#'     \item compare_mode: Character, comparison mode ("Comparative" or "Temporal")
#'   }
#' @slot niche_info List of data frames containing niche metadata for each
#'   niche key across samples
#' @slot niche_cells List of data frames containing cell-level information for
#'   each niche key across samples
#' @slot niche_genes List of data frames containing gene-level information for
#'   each niche key across samples
#'
#' @name MultiSampNiche-class
#' @rdname MultiSampNiche-class
#' @exportClass MultiSampNiche
#'
#' @family niche-classes
#'
#' @aliases MultiSampNiche
#'
MultiSampNiche <- setClass(
  Class = "MultiSampNiche",
  slots = list(
    samp_info = "list",
    niche_info = "list",
    niche_cells = "list",
    niche_genes = "list"
  ),
  prototype = prototype(
    samp_info = list(
      multi_id = character(), # multi_index
      samp_id = character(),
      samp_grp = character(),
      compare_mode = character() # "Comparative"/"Temporal"
    ),
    niche_info = list(),
    niche_cells = list(),
    niche_genes = list()
  )
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# STID_analysis
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The STID_analysis Class
#'
#' A container for all analysis results and metadata within an STID object.
#' Stores project information, metadata records, and both single-sample and
#' multi-sample niche analysis results.
#'
#' @slot STID_info List containing comprehensive project information:
#'   \itemize{
#'     \item samp_info: Sample-level metadata
#'     \item data_info: Data structure information
#'     \item project_info: Project-level metadata
#'     \item comment_info: Additional comments and notes
#'   }
#' @slot meta_data_record List containing metadata tracking information:
#'   \itemize{
#'     \item meta_data_info: Data frame tracking all metadata additions
#'     \item meta_data_list: List of all stored metadata data frames
#'   }
#' @slot SingleSampNiche List of SingleSampNiche objects, one per sample
#' @slot MultiSampNiche List of MultiSampNiche objects for cross-sample analyses
#'
#' @name STID_analysis-class
#' @rdname STID_analysis-class
#' @exportClass STID_analysis
#'
#' @family STID-classes
#'
#' @aliases STID_analysis
#'
STID_analysis <- setClass(
  Class = "STID_analysis",
  slots = list(
    STID_info = "list",
    meta_data_record = "list",
    SingleSampNiche = "list",
    MultiSampNiche  = "list"
  ),
  prototype = prototype(
    STID_info = list(
      samp_info = list(
        samp_id = character(),
        infect_time = character(),
        host_org = character(),
        pathogen_grp = character(),
        pathogen_org = character()
      ),
      data_info = list(
        samp_colnm = character(),
        samp_grp_colnm = character(),
        celltype_colnm = character(),
        x_colnm = character(),
        y_colnm = character(),
        pathogen_genes = character(),
        data_format = character(),
        data_platform = character(),
        binsize = numeric(),
        coord_interval = numeric(),
        base_unit = numeric() # * mm
      ),
      project_info = list(
        project_id = character(),
        description = character()
      ),
      comment_info = list()
    ),
    meta_data_record = list(
      meta_data_info = data.frame(
        meta_key = character(),
        time = character(),
        func_nm = character(),
        dir_nm = character(),
        grp_nm = character(),
        asso_key = character(),
        description = character(),
        stringsAsFactors = FALSE,
        row.names = NULL
      ),
      meta_data_list = list()
    ),
    SingleSampNiche = list(),
    MultiSampNiche  = list()
  )
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# STID Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' The STID Class
#'
#' The STID (Spatial Transcriptomics toolkit for Infectious Diseases) object extends the Seurat
#' class to add specialized functionality for analyzing infection and spatial
#' transcriptomics data. It inherits all functionality from Seurat objects while
#' adding additional slots for tracking infection-related analyses.
#'
#' @slot STID_analysis An \code{\link{STID_analysis}} object containing all
#'   infection-specific analysis results and metadata
#' @slot assays A list of assays for this project (inherited from Seurat)
#' @slot meta.data Contains meta-information about each cell (inherited from Seurat)
#' @slot active.assay Name of the active, or default, assay (inherited from Seurat)
#' @slot active.ident The active cluster identity (inherited from Seurat)
#' @slot graphs A list of \code{\link[SeuratObject]{Graph}} objects (inherited from Seurat)
#' @slot neighbors Neighbor graphs (inherited from Seurat)
#' @slot reductions A list of dimensional reduction objects (inherited from Seurat)
#' @slot images A list of spatial image objects (inherited from Seurat)
#' @slot project.name Name of the project (inherited from Seurat)
#' @slot misc A list of miscellaneous information (inherited from Seurat)
#' @slot version Version of Seurat this object was built under (inherited from Seurat)
#' @slot commands A list of logged commands (inherited from Seurat)
#' @slot tools A list of miscellaneous data generated by other tools (inherited from Seurat)
#'
#' @name STID-class
#' @rdname STID-class
#' @exportClass STID
#'
#' @family STID-classes
#'
#' @aliases STID
#'
#' @examples
#' \dontrun{
#' # Create an STID object from a Seurat object
#' ist_obj <- as.STID(
#'   seurat_obj = seurat_object,
#'   samp_colnm = "sample_id",
#'   celltype_colnm = "cell_type",
#'   host_org = "human",
#'   pathogen_grp = "virus",
#'   pathogen_org = "SARS-CoV-2",
#'   data_format = "StereoSeq"
#' )
#' }
#'
STID <- setClass(
  Class = 'STID',
  contains = 'Seurat', # 需要加载Seurat才行
  slots = list(
    STID_analysis = 'STID_analysis' # 默认第一位
  ),
  prototype = prototype(
    STID_analysis = STID_analysis()
    )
)
# getClass("STID") # view the class definition
setMethod("initialize", "STID", # 对STID使用seurat的函数，返回的是seurat对象，但是最终还会被转换成STID对象，所以会频繁触发初始化函数
          function(.Object,
                   STID_analysis = STID_analysis(),
                   assays = list(),
                   meta.data = NULL,
                   active.assay = character(length = 0L),
                   active.ident = NULL,
                   graphs = list(),
                   neighbors = list(),
                   reductions = list(),
                   images = list(),
                   project.name = getOption(
                     x = 'Seurat.object.project',
                     default = Seurat.options$Seurat.object.project
                   ),
                   misc = list(),
                   version = packageVersion(pkg = 'SeuratObject'),
                   commands = list(),
                   tools = list(),
                   ...) {
            # clog_normal("Initializing STID object...") # seurat有些函数会频繁触发
            # .Object <- callNextMethod(.Object, ...) # 自定义initialize就可以禁用父类的initialize方法
            .Object@STID_analysis <- STID_analysis
            .Object@assays <- assays
            .Object@meta.data <- meta.data
            .Object@active.assay <- active.assay
            .Object@active.ident <- active.ident
            .Object@graphs <- graphs
            .Object@neighbors <- neighbors
            .Object@reductions <- reductions
            .Object@images <- images
            .Object@project.name <- project.name
            .Object@misc <- misc
            .Object@version <- version
            .Object@commands <- commands
            .Object@tools <- tools
            return(.Object)
          }
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# as.STID method: Convert Seurat to STID
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert Seurat object to STID object
#'
#' Converts a Seurat object to an STID (Integrated Single-cell Transcriptomics)
#' object, adding infection-specific metadata and analysis structures.
#'
#' @param seurat_obj A Seurat object to be converted
#' @param host_org Character, host organism - "human", "mouse", or "unknown"
#' @param pathogen_grp Character, pathogen group - "virus", "bacteria", "parasite", or "unknown"
#' @param pathogen_org Character, specific pathogen organism name
#' @param samp_colnm Character, column name in metadata containing sample IDs
#' @param samp_grp_colnm Character, column name in metadata containing sample groups
#' @param celltype_colnm Character, column name in metadata containing cell types
#' @param x_colnm Character, column name for x coordinates
#' @param y_colnm Character, column name for y coordinates
#' @param pathogen_genes Character vector, names of pathogen genes
#' @param data_format Character, spatial data format.
#' One of "square_grid", "hex_grid", or "single_cell".
#' For "single_cell", the default spatial unit, binsize,
#' and coord_interval are 1 unless otherwise specified.
#' @param data_platform Character, data platform - "StereoSeq", "Visium", "VisiumHD", "SlideSeq", or "unknown"
#' @param binsize Numeric, bin size for spatial data
#' @param coord_interval Numeric, coordinate coord_interval
#' @param base_unit Numeric, base unit for spatial coordinates, and the unit is mm
#' @param project_id Character, project identifier
#' @param description Character, project description
#' @param ... Additional arguments passed to methods
#'
#' @return An STID object containing the converted data with infection-specific
#'   metadata and analysis structures
#'
#' @rdname as.STID
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert Seurat object to STID object
#' ist_obj <- as.STID(
#'   seurat_obj = seurat_object,
#'   samp_colnm = "sample_id",
#'   celltype_colnm = "cell_type",
#'   host_org = "human",
#'   pathogen_grp = "virus",
#'   pathogen_org = "SARS-CoV-2",
#'   data_format = "square_grid",
#'   data_platform = "StereoSeq",
#'   binsize = 1,
#'   coord_interval = 1,
#'   base_unit = 0.5,
#'   project_id = "STID_Project_001",
#'   description = "This is a test project for converting Seurat to STID."
#' )
#' }
#'
setGeneric("as.STID", function(seurat_obj, ...) standardGeneric("as.STID")) # needs to be defined manually

#' @rdname as.STID
#' @export
setMethod("as.STID", "Seurat",
          function(seurat_obj = NULL, host_org = NULL, pathogen_grp = NULL, pathogen_org = NULL,
                   samp_colnm = NULL, samp_grp_colnm = NULL, celltype_colnm = NULL,
                   x_colnm = NULL, y_colnm = NULL,pathogen_genes = NULL,
                   data_format = NULL, data_platform = NULL,
                   binsize = NULL, coord_interval = NULL, base_unit = NULL,
                   project_id = NULL, description = NULL,...
                   ) {

            # >>> Start pipeline
            clog_start()

            # >>> Check input patameter
            clog_check()
            if (!inherits(seurat_obj, "Seurat")) {
              clog_error("Input object is not a Seurat object.")
            }
            .check_not_any_null(samp_colnm, samp_grp_colnm, host_org, pathogen_grp,
                             pathogen_org, data_format, data_platform) # celltype_colnm/pathogen_genes can be NULL
            match.arg(host_org, choices = .STID_globals$supported_hosts, several.ok = FALSE)
            match.arg(pathogen_grp, choices = .STID_globals$supported_pathogens, several.ok = FALSE)
            match.arg(data_format, choices = .STID_globals$supported_formats, several.ok = FALSE)
            match.arg(data_platform, choices = .STID_globals$supported_platforms, several.ok = FALSE)
            if(!is.null(pathogen_genes)){
              if(!all(pathogen_genes %in% rownames(seurat_obj))){
                clog_error(paste0("Not all pathogen_genes found in seurat_obj rownames: ",
                                  paste(setdiff(pathogen_genes, rownames(seurat_obj)), collapse = ", ")))
              }
            }

            #>
            clog_normal("UpdateSeuratObject to ensure compatibility with the latest Seurat version")
            seurat_obj <- suppressMessages({
              UpdateSeuratObject(seurat_obj)
            })

            #> samp_colnm, samp_grp_colnm, celltype_colnm
            meta_data <- seurat_obj@meta.data
            meta_colnm <- colnames(meta_data)
            .check_column_exist(meta_data, samp_colnm, samp_grp_colnm, celltype_colnm)
            if(is.factor(seurat_obj@meta.data[[samp_colnm]])){
              clog_normal("Using existing sample factor levels")
              samp_id <- levels(seurat_obj@meta.data[[samp_colnm]])
            }else{
              clog_normal("Using the default sample order")
              samp_id <- seurat_obj@meta.data[[samp_colnm]] %>% as.character() %>% unique() # %>% sort()
            }

            #> x_colnm and y_colnm
            if("x" %in% meta_colnm & "y" %in% meta_colnm){
              clog_normal("Using existing 'x' and 'y' columns in meta.data for spatial coordinates")
            } else if(is.null(x_colnm)|is.null(y_colnm)){
              clog_normal("You didn't specify x_colnm or y_colnm,
                          adding spatial coordinates from images to meta.data as 'x' and 'y'")
              seurat_obj@meta.data <- .add_coord2metadata(seurat_obj, meta_data = NULL,
                                                          x_colnm = "x", y_colnm = "y")$meta_data
            }else if(x_colnm != "x" | y_colnm != "y"){
              if( !x_colnm %in% meta_colnm | !y_colnm %in% meta_colnm){
                clog_error(paste0("x_colnm or y_colnm not found in meta.data: ",
                                  paste(meta_colnm, collapse = ", ")))
              }else{
                seurat_obj@meta.data <- seurat_obj@meta.data %>%
                  dplyr::rename(x = all_of(x_colnm), y = all_of(y_colnm))
                clog_warn(paste0("Renaming x_colnm to 'x' and y_colnm to 'y' in meta.data")) # !!!!
              }
            }

            #> base_unit
            if(is.null(base_unit)){
              clog_warn("You didn't specify base_unit, using the default base_unit for your data_platform")
              if(data_platform == "StereoSeq"){
                base_unit <- 0.5 # mm
              }else if(data_platform == "VisiumHD"){
                base_unit <- 2
              }else if(data_platform == "SlideSeq"){
                base_unit <- 10
              }else if(data_platform == "Visium"){
                base_unit <- 65
              }else{
                clog_error(paste0("The default base_unit for data_platform ", data_platform, " is not defined.
                                  Please specify the base_unit for your data_platform."))

              }
            }

            # convert the coord_interval of x and y to 1
            if(data_format == "square_grid"){
              coord_interval_x <- seurat_obj@meta.data$x %>% unique() %>% sort() %>% diff() %>%
                table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% .[1,1] %>% droplevels() %>% as.numeric()
              coord_interval_y <- seurat_obj@meta.data$y %>% unique() %>% sort() %>% diff() %>%
                table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% .[1,1] %>% droplevels() %>% as.numeric()
              if(coord_interval_x != coord_interval_y){
                clog_error(paste0("Your x and y coordinates have inconsistent coord_intervals: ",
                                  coord_interval_x, " vs ", coord_interval_y, "."))
              }
              if(is.null(binsize)){
                clog_warn(paste0("You didn't specify binsize, using the coord_interval of x and y coordinates: ",
                                 coord_interval_x))
                binsize <- coord_interval_x
                if(coord_interval_x == 1){
                  clog_warn("Your x and y coordinates already have an coord_interval of 1.
                            Please check if your coordinates are already corrected.")
                }
              }
              clog_normal("For StereoSeq data, we will convert the coord_interval of x and y to 1")
              seurat_obj@meta.data[["x"]] <- seurat_obj@meta.data[["x"]]/coord_interval_x %>%
                `+`(0.0001) %>% round(digits = 0)
              seurat_obj@meta.data[["y"]] <- seurat_obj@meta.data[["y"]]/coord_interval_y %>%
                `+`(0.0001) %>% round(digits = 0)
              coord_interval_x <- seurat_obj@meta.data$x %>% unique() %>% sort() %>% diff() %>%
                table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% .[1,1] %>% as.numeric()
              coord_interval_y <- seurat_obj@meta.data$y %>% unique() %>% sort() %>% diff() %>%
                table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% .[1,1] %>% as.numeric()
              if(coord_interval_x != 1 | coord_interval_y != 1){
                clog_error(paste0("After binning, your x and y binsize are still different from 1: ",
                                  coord_interval_x, " vs ", coord_interval_y, "."))
              }
              coord_interval <- rep(1, length(samp_id))
              names(coord_interval) <- samp_id
            }else if(data_format == "hex_grid" & data_platform == "Visium"){
              if(is.null(binsize)){
                clog_warn(paste0("You didn't specify binsize, using the default binsize for Visium data: 200"))
                binsize <- 200
              }
              coord_interval <- seurat_obj@meta.data %>% group_by(!!sym(samp_colnm)) %>%
                summarise(scale.factors_spot = unique(scale.factors_spot)) %>%
                mutate(!!sym(samp_colnm) := factor(!!sym(samp_colnm), levels = samp_id)) %>%
                arrange(!!sym(samp_colnm)) %>%
                 pull(scale.factors_spot)
              coord_interval <- coord_interval/base_unit*100
              names(coord_interval) <- samp_id
            }else if(data_format == "single_cell"){
              if(is.null(binsize)){
                clog_warn(paste0("You didn't specify binsize, using the default binsize for single_cell data: 1"))
                binsize <- 1
              }
              coord_interval <- rep(1, length(samp_id))
              names(coord_interval) <- samp_id
            }else{
              clog_error(paste0("The plat_format ", data_format, " of data_platform ", data_platform,
                                " is not common in current STID analysis. Please check if your data_format and data_platform are correct."))
            }

            #>
            clog_normal(paste("Using samp_colnm:", samp_colnm, ", samp_grp_colnm:", samp_grp_colnm,
                              ", celltype_colnm:", celltype_colnm,", x_colnm: x , y_colnm: y"))
            clog_normal(paste("Using host_org:", host_org, ", pathogen_grp:", pathogen_grp, ", pathogen_org:", pathogen_org))
            clog_normal(paste("Using data_format:", data_format, ", data_platform:", data_platform,
                              ", binsize:", binsize,", coord_interval:", paste(coord_interval, collapse = ", "), "base_unit:", base_unit))
            # >>> End check

            # >>> Start main pipeline
            # Initialize STID_analysis container
            clog_step("Converting Seurat object to STID object...")
            STID_analysis_obj <- STID_analysis(
              STID_info = list(
                samp_info = list(
                  samp_id = samp_id,
                  infect_time = samp_id,
                  host_org = host_org,
                  pathogen_grp = pathogen_grp,
                  pathogen_org = pathogen_org
                ),
                data_info = list(
                  samp_colnm = samp_colnm,
                  samp_grp_colnm = samp_grp_colnm,
                  celltype_colnm = celltype_colnm,
                  x_colnm = "x",
                  y_colnm = "y",
                  pathogen_genes = pathogen_genes,
                  data_format = data_format,
                  data_platform = data_platform,
                  binsize = binsize,
                  coord_interval = coord_interval,
                  base_unit = base_unit
                ),
                project_info = list(
                  project_id = project_id,
                  description = description
                ),
                comment_info = list()
              ),
              meta_data_record = list(
                meta_data_info = data.frame(
                  meta_key = character(),
                  time = character(),
                  func_nm = character(),
                  dir_nm = character(),
                  grp_nm = character(),
                  asso_key = character(),
                  description = character(),
                  stringsAsFactors = FALSE,
                  row.names = NULL
                ),
                meta_data_list = list()
              ),
              SingleSampNiche = list(),
              MultiSampNiche  = list()
            )

            # Create STID object
            clog_normal("Conver the gene names: replace '_' with '-'")
            rownames(seurat_obj) <- gsub("_", "-", rownames(seurat_obj)) # !!! can modify searut, but cannot modify STID?
            obj <- new("STID",
                       STID_analysis = STID_analysis_obj,
                       assays = seurat_obj@assays,
                       meta.data = seurat_obj@meta.data,
                       active.assay = seurat_obj@active.assay,
                       active.ident = seurat_obj@active.ident,
                       graphs = seurat_obj@graphs,
                       neighbors = seurat_obj@neighbors,
                       reductions = seurat_obj@reductions,
                       images = seurat_obj@images,
                       project.name = seurat_obj@project.name,
                       misc = seurat_obj@misc,
                       version = seurat_obj@version,
                       commands = seurat_obj@commands,
                       tools = seurat_obj@tools
            )

            clog_end()
            return(obj)
          }
)

#' Add spatial coordinates from images to metadata
#'
#' Internal function to extract spatial coordinates from STID object images and
#' add them to the metadata.
#'
#' @param STID_obj An STID object containing spatial images
#' @param meta_data Data frame, existing metadata (if NULL, uses STID_obj@meta.data)
#' @param x_colnm Character, name for the x-coordinate column in output
#' @param y_colnm Character, name for the y-coordinate column in output
#' @param image_names Character vector, names of images to extract coordinates from
#'        (if NULL, uses all images in the STID object)
#'
#' @return List containing:
#'   \itemize{
#'     \item meta_data: Updated metadata with added coordinate columns
#'     \item spatial_coords: Data frame of extracted spatial coordinates
#'   }
#'
#' @keywords internal
#'
#' @noRd
.add_coord2metadata <- function(STID_obj = NULL, meta_data = NULL,
                                x_colnm = "x", y_colnm = "y", image_names = NULL) {
  if(is.null(image_names)){
    clog_normal("You didn't specify image_names, using all images in the STID object")
    image_names <- names(STID_obj@images)
    if(length(image_names) > 1){
      clog_warn(paste0("Your STID object contains multiple images: ", paste(image_names, collapse = ", "),
                       ". We will add coordinates from all images to the meta.data"))
    }
  }
  coord_list <- list()
  for (img_name in image_names) {
    img <- STID_obj@images[[img_name]]

    if(inherits(img, "VisiumV1")) {
      coords <- as.data.frame(img@coordinates)[, c("imagerow", "imagecol")]
      colnames(coords) <- c(y_colnm, x_colnm)
    }else if(inherits(img, "VisiumV2")){
      coords <- as.data.frame(img@boundaries[["centroids"]]@coords)
      rownames(coords) <- img@boundaries[["centroids"]]@cells
      colnames(coords) <- c( x_colnm, y_colnm)
    }else{
      clog_error(paste0("Image ", img_name, " is not of class VisiumV1 or VisiumV2"))
    }
    if(nrow(coords)==0){
      clog_warn(paste0("Image ", img_name, " has no coordinates, skipping..."))
      next
    }
    coords$img_name <- img_name
    coords$scale.factors_spot <- img@scale.factors$spot
    coords$cell_id <- rownames(coords)
    coord_list[[img_name]] <- coords
  }
  spatial_coords <- do.call(rbind, coord_list) %>%
    as.data.frame() %>%
    remove_rownames() %>%
    # rownames_to_column("merge_rownm") %>%
    column_to_rownames(var = "cell_id") %>%
    mutate(cell_id = rownames(.))
  if(missing(meta_data) | is.null(meta_data)){
    meta_data <- STID_obj@meta.data
  }
  if(!all(rownames(meta_data) %in% spatial_coords$cell_id)){
    clog_error("Not all cells in meta.data have spatial coordinates in the STID object images")
  }
  meta_data <- meta_data %>%
    mutate(cell_id = rownames(.)) %>%
    left_join(spatial_coords[c("img_name",x_colnm,y_colnm,"scale.factors_spot","cell_id")], by = "cell_id") %>%
    column_to_rownames(var = "cell_id") %>%
    as.data.frame()
  return(list(meta_data = meta_data, spatial_coords = spatial_coords))
}


#' Convert Seurat object to STID object
#'
#' Creates an STID object from a Seurat object with infection-specific metadata.
#'
#' @inheritParams as.STID
#' @param ... Additional arguments passed to \code{\link{as.STID}}
#'
#' @return An STID object
#'
#' @rdname CreateSTIDObject
#' @export
#'
#' @examples
#' \dontrun{
#' ist_obj <- CreateSTIDObject(
#'   seurat_obj = seurat_object,
#'   samp_colnm = "sample_id",
#'   celltype_colnm = "cell_type",
#'   host_org = "human",
#'   pathogen_grp = "virus",
#'   pathogen_org = "SARS-CoV-2",
#'   data_format = "square_grid",
#'   data_platform = "StereoSeq"
#' )
#' }
CreateSTIDObject <- function(
    seurat_obj = NULL,
    host_org = NULL,
    pathogen_grp = NULL,
    pathogen_org = NULL,
    samp_colnm = NULL,
    samp_grp_colnm = NULL,
    celltype_colnm = NULL,
    x_colnm = NULL,
    y_colnm = NULL,
    pathogen_genes = NULL,
    data_format = NULL,
    data_platform = NULL,
    binsize = NULL,
    coord_interval = NULL,
    base_unit = NULL,
    project_id = NULL,
    description = NULL,
    ...
) {
  as.STID(
    seurat_obj = seurat_obj,
    host_org = host_org,
    pathogen_grp = pathogen_grp,
    pathogen_org = pathogen_org,
    samp_colnm = samp_colnm,
    samp_grp_colnm = samp_grp_colnm,
    celltype_colnm = celltype_colnm,
    x_colnm = x_colnm,
    y_colnm = y_colnm,
    pathogen_genes = pathogen_genes,
    data_format = data_format,
    data_platform = data_platform,
    binsize = binsize,
    coord_interval = coord_interval,
    base_unit = base_unit,
    project_id = project_id,
    description = description,
    ...
  )
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# as.Seurat
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert STID object to Seurat object
#'
#' Converts an STID object back to a standard Seurat object, preserving all
#' assay data and basic metadata but removing STID-specific analysis structures.
#'
#' @param STID_obj An STID object to be converted
#'
#' @return A Seurat object containing the same assay data and basic metadata
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert STID object back to Seurat
#' seurat_obj <- as.Seurat.STID(ist_obj)
#' }
#'
as.Seurat.STID <- function(STID_obj){
  clog_normal("Converting STID object to Seurat object...")
  obj <- new("Seurat",
             assays = STID_obj@assays,
             meta.data = STID_obj@meta.data,
             active.assay = STID_obj@active.assay,
             active.ident = STID_obj@active.ident,
             graphs = STID_obj@graphs,
             neighbors = STID_obj@neighbors,
             reductions = STID_obj@reductions,
             images = STID_obj@images,
             project.name = STID_obj@project.name,
             misc = STID_obj@misc,
             version = STID_obj@version,
             commands = STID_obj@commands,
             tools = STID_obj@tools
  )
  return(obj)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CreateSingleSampNiche
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Create a SingleSampNiche object from analysis results
#'
#' Creates and stores a SingleSampNiche object within an STID object, organizing
#' niche detection results for a single sample. This function extracts relevant
#' cell and gene information based on specified niche parameters.
#'
#' @param STID_obj An STID object containing analysis results
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param niche_key Character, unique identifier for the niche analysis
#' @param meta_key Character, metadata key containing the source data
#' @param ROI_type Character, type of ROI - "ROI", "Spot", or "Unknown"
#' @param pos_colnm Character, column name containing positive spot labels
#' @param neg_value Character, value indicating negative spots (default: "neg")
#' @param center_colnm Character, column name indicating center spots (optional)
#' @param edge_colnm Character, column name indicating edge spots (optional)
#' @param all_label_colnm Character, column name for all spot labels (optional)
#' @param all_dist_colnm Character, column name for distance to center (optional)
#' @param other_colnm Character vector, additional columns to include (optional)
#' @param description Character, description of the niche analysis (optional)
#' @param ... Additional arguments passed to methods
#'
#' @return Modified STID object with added SingleSampNiche information
#'
#' @import Seurat
#' @import rlang
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create SingleSampNiche object
#' STID_obj <- CreateSingleSampNiche(
#'   STID_obj = ist_obj,
#'   niche_key = "niche_virulence",
#'   meta_key = "M2_NicheDetect_STS_20240101",
#'   ROI_type = "ROI",
#'   pos_colnm = "ROI_label",
#'   center_colnm = "ROI_center",
#'   edge_colnm = "ROI_edge",
#'   all_label_colnm = "All_ROI_label",
#'   all_dist_colnm = "All_Dist2ROIcenter"
#' )
#' }
#'
CreateSingleSampNiche <- function(
    STID_obj = NULL,
    loop_id = "LoopAllSamp",
    niche_key = NULL,
    meta_key = NULL,
    ROI_type = NULL,
    pos_colnm = NULL,
    neg_value = "neg",
    center_colnm = NULL,
    edge_colnm = NULL,
    all_label_colnm = NULL,
    all_dist_colnm = NULL,
    other_colnm = NULL,
    description = NULL
) {
  UseMethod("CreateSingleSampNiche", STID_obj)
}

#' @rdname CreateSingleSampNiche
#' @export
CreateSingleSampNiche.STID <- function(
    STID_obj = NULL,
    loop_id = "LoopAllSamp",
    niche_key = NULL,
    meta_key = NULL,
    ROI_type = NULL,
    pos_colnm = NULL,
    neg_value = "neg",
    center_colnm = NULL,
    edge_colnm = NULL,
    all_label_colnm = NULL,
    all_dist_colnm = NULL,
    other_colnm = NULL,
    description = NULL
) {
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_not_any_null(niche_key, meta_key, ROI_type,pos_colnm)
  if(!ROI_type %in% c("ROI","Spot","Unknown")){
    clog_error(paste0("ROI_type must be 'ROI','Spot' or 'Unknown', not '", ROI_type, "'"))
  }
  if(ROI_type %in% c("ROI","Spot")){
    if(is.null(center_colnm) | is.null(edge_colnm) | is.null(all_label_colnm) | is.null(all_dist_colnm)){
      clog_error("For ROI or Spot type, center_colnm, edge_colnm, all_label_colnm and all_dist_colnm cannot be NULL")
    }
  }
  now_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  loop_single <- .check_loop_single(STID_obj = STID_obj,loop_id = loop_id, mode = 1)
  # >>> End check

  # >>> Start main pipeline
  clog_step("Creating CreateSingleSampNiche object...")

  #> niche_cells
  meta_unlist <- meta_key %>% unlist()
  if(length(meta_unlist) != 1){
    clog_warn("The meta_key contains multiple values. GetMetaData will use all meta_key,
              but only the first meta_key will recorded in niche_info")
  }
  meta_data <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]]
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  coord_colnm <- GetCoordInfo(STID_obj) %>% na.omit()
  .check_column_exist(meta_data, samp_colnm, coord_colnm, pos_colnm, center_colnm,
                      edge_colnm, all_label_colnm, all_dist_colnm, other_colnm)
  for(i in seq_along(loop_single)){ # !!! not valid_single
    i_single <- loop_single[i]
    clog_loop(paste0("Processing samp: ", i_single, " (", i, "/", length(loop_single), ")"))

    # niche_cells
    i_meta_data <- meta_data %>% filter(!!sym(samp_colnm) == i_single)
    niche_cells <- i_meta_data[, c(coord_colnm, pos_colnm, center_colnm, edge_colnm,
                                   all_label_colnm, all_dist_colnm, other_colnm)]
    pos_meta_data <- i_meta_data %>%
      filter(.,!!sym(pos_colnm) != neg_value)
    len_pos_value <- nrow(pos_meta_data)
    if(len_pos_value == 0){
      niche_info <- data.frame()
      niche_cells <- list()
      niche_genes <- list()
    }else{
      niche_cells <- niche_cells %>%
        mutate(is_Niche = if_else(!!sym(pos_colnm) != neg_value, TRUE, FALSE), .before = 1)

      #> niche_genes
      count_all <- GetAssayData(STID_obj, layer = "count")
      count_pos <- count_all[,rownames(pos_meta_data)]
      niche_genes <- data.frame(gene = rownames(STID_obj),
                                avg_count_all = Matrix::rowMeans(count_all),
                                detection_rate_all = Matrix::rowMeans(count_all > 0),
                                avg_count_pos = Matrix::rowMeans(count_pos),
                                detection_rate_pos = Matrix::rowMeans(count_pos > 0),
                                stringsAsFactors = FALSE,
                                row.names = rownames(STID_obj))

      #> niche_info
      niche_info <- data.frame(
        niche_key = niche_key,
        meta_key = meta_unlist[1],
        ROI_type = ROI_type,
        pos_colnm = pos_colnm,
        center_colnm = center_colnm %>% .null_to_na(),
        edge_colnm = edge_colnm %>% .null_to_na(),
        all_label_colnm = all_label_colnm %>% .null_to_na(),
        all_dist_colnm = all_dist_colnm %>% .null_to_na(),
        niche_cells_num = nrow(pos_meta_data) %>% .null_to_0(),
        niche_genes_num = nrow(niche_genes) %>% .null_to_0(),
        description =  description %>% .null_to_na(),
        stringsAsFactors = FALSE,
        row.names = niche_key
      )
    }

    # SingleSampNiche_obj
    singlesamp_nms <- STID_obj@STID_analysis@SingleSampNiche %>% names()
    if(!i_single %in% singlesamp_nms){
      if(!is_empty(niche_cells)){
        clog_normal(paste0("Creating SingleSampNiche object for samp: ", i_single))
        clog_normal(paste0("Adding niche_key: ", niche_key))
        niche_cells <- list(niche_cells) %>% setNames(niche_key)
        niche_genes <- list(niche_genes) %>% setNames(niche_key)
      }else{
        clog_warn(paste0("No niche_cells found for samp: ", i_single, ", creating empty SingleSampNiche object..."))
      }
      SingleSampNiche_obj <- new("SingleSampNiche",
                                 samp_info = list(
                                   samp_id = i_single
                                 ),
                                 niche_info = niche_info,
                                 niche_cells = niche_cells,
                                 niche_genes = niche_genes
      )
      STID_obj@STID_analysis@SingleSampNiche[[i_single]] <- SingleSampNiche_obj
    }else{
      if(!is_empty(niche_cells)){
        clog_normal(paste0("SingleSampNiche object already exists for samp: ", i_single))
        STID_niche_info <- STID_obj@STID_analysis@SingleSampNiche[[i_single]]@niche_info
        if(niche_key %in% STID_niche_info$niche_key){
          clog_warn(paste0("niche_key: ", niche_key, " already exists, overwriting..."))
        }else{
          clog_normal(paste0("Adding niche_key: ", niche_key))
        }
        if(!is_empty(STID_niche_info)){
          STID_obj@STID_analysis@SingleSampNiche[[i_single]]@niche_info <-
            STID_niche_info %>%
            filter(., .data$niche_key != .env$niche_key) %>%
            bind_rows(niche_info)
        }else{
          STID_obj@STID_analysis@SingleSampNiche[[i_single]]@niche_info <- niche_info
        }
        STID_obj@STID_analysis@SingleSampNiche[[i_single]]@niche_cells[[niche_key]] <- niche_cells
        STID_obj@STID_analysis@SingleSampNiche[[i_single]]@niche_genes[[niche_key]] <- niche_genes
      }else{
        clog_warn(paste0("No niche_cells found for samp: ", i_single, ", skipping..."))
      }
    }
  }

  # >>> Final
  clog_end()
  return(STID_obj)
}

#' @rdname CreateSingleSampNiche
#' @export
CreateSSNiche <- function(...){
  CreateSingleSampNiche(...)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CreateMultiSampNiche
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a MultiSampNiche object for cross-sample comparison
#'
#' Creates and stores a MultiSampNiche object within an STID object, combining
#' niche analysis results from multiple samples for comparative analysis.
#'
#' @param STID_obj An STID object containing SingleSampNiche objects
#' @param multi_id Character, unique identifier for the multi-sample analysis
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param compare_mode Character, comparison mode - "Comparative" or "Temporal"
#' @param niche_key Character, niche key to combine across samples
#' @param description Character, description of the multi-sample analysis
#' @param ... Additional arguments passed to methods
#'
#' @return Modified STID object with added MultiSampNiche information
#'
#' @import Seurat
#' @import rlang
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create MultiSampNiche object for comparative analysis
#' STID_obj <- CreateMultiSampNiche(
#'   STID_obj = ist_obj,
#'   multi_id = "comparison_infected_vs_control",
#'   loop_id = c("infected1", "infected2", "control1", "control2"),
#'   compare_mode = "Comparative",
#'   niche_key = "niche_virulence"
#' )
#' }
#'
CreateMultiSampNiche <- function(
    STID_obj = NULL,
    multi_id = NULL,
    loop_id = "LoopAllSamp",
    compare_mode = NULL,
    niche_key = NULL,
    description = NULL
) {
  UseMethod("CreateMultiSampNiche", STID_obj)
}

#' @rdname CreateMultiSampNiche
#' @export
CreateMultiSampNiche.STID <- function(
    STID_obj = NULL,
    multi_id = NULL,
    loop_id = "LoopAllSamp",
    compare_mode = NULL,
    niche_key = NULL,
    description = NULL
) {
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_not_any_null(compare_mode,niche_key)
  match.arg(compare_mode, choices = c("Comparative", "Temporal"))
  now_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  valid_single <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "samp_id")[[1]] # !!! must need
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)

  #>
  meta_data <- GetMetaData(STID_obj = STID_obj, meta_key = "raw")[[1]]
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  samp_grp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_grp_colnm")[[1]]
  samp_grp_df <- meta_data %>% dplyr::select(all_of(samp_colnm), all_of(samp_grp_colnm)) %>% distinct()
  samp_grp <- meta_data[[samp_grp_colnm]][match(loop_single, meta_data[[samp_colnm]])] %>%
    factor(.,levels = .)
  if(length(loop_single) != length(samp_grp)){
    clog_error("loop_single and samp_grp must have the same length")
  }
  len_samp_grp <- length(unique(samp_grp))
  if(len_samp_grp < 2){
    clog_error("samp_grp must contain at least two groups")
  }
  if(compare_mode == "Comparative"){
    clog_warn("The first group in samp_grp will be used as the control group")
    if(len_samp_grp != 2){
      clog_error("For Comparative mode, samp_grp must contain exactly two groups")
    }
  }

  #>
  valid_single_index <- match(loop_single, valid_single ) %>% sort() # !!! must need
  if(is.null(multi_id)){
    # multi_id <- paste0("multi_", paste(loop_single, collapse = "_"))
    multi_id <- paste0(compare_mode,"_",valid_single_index %>% paste(collapse = "_"))
    clog_warn("Your multi_id is NULL, using default multi_id")
  }
  clog_normal(paste0("Your multi_id: ", multi_id))
  # >>> End check

  # >>> Start main pipeline
  clog_step("Creating CreateMultiSampNiche object...")

  #> samp_info
  samp_info <- list(
    multi_id = multi_id,
    samp_id = loop_single,
    samp_grp = samp_grp,
    compare_mode = compare_mode
  )

  #> niche_info
  # browser()
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  niche_info <- lapply(STID_obj@STID_analysis@SingleSampNiche, function(x){
    .niche_info <- x@niche_info
    if(!niche_key %in% .niche_info$niche_key){
      return(data.frame())
    }else{
      .niche_info <- .niche_info[niche_key, , drop = FALSE]
      return(.niche_info)
    }
  }) %>% bind_rows(.id = samp_colnm)
  if(nrow(niche_info) == 0){
    clog_error(paste0("niche_key: ", niche_key, " not found in any samp"))
  }
  rownames(niche_info) <- niche_info[,samp_colnm]
  if(!all(loop_single %in% rownames(niche_info))){
    clog_warn(paste0("niche_key: ", niche_key, " not found in samp: ",
                      paste(setdiff(loop_single, rownames(niche_info)), collapse = ", ")))
    niche_info <- list()
    niche_cells <- list()
    niche_genes <- list()
  }else{
    niche_info <- niche_info[loop_single, , drop = FALSE]
    # GetSSNicheCells?
    niche_cells <- lapply(STID_obj@STID_analysis@SingleSampNiche, function(x){
      x@niche_cells[[niche_key]]
    })[loop_single] %>%
      bind_rows()
    niche_genes <- data.frame(gene = rownames(STID_obj),
                              stringsAsFactors = FALSE,
                              row.names = rownames(STID_obj))

  }

  #> MultiSampNiche_obj
  valid_multi <- STID_obj@STID_analysis@MultiSampNiche %>% names()
  if(!multi_id %in% valid_multi){
    if(!is_empty(niche_cells)){
      clog_normal(paste0("Creating MultiSampNiche object for multi_id: ", multi_id))
      niche_info <- list(niche_info) %>% setNames(niche_key)
      niche_cells <- list(niche_cells) %>% setNames(niche_key)
      niche_genes <- list(niche_genes) %>% setNames(niche_key)
    }else{
      clog_warn(paste0("Some samps don't exist niche_cells for niche_key: ", niche_key, ", creating empty MultiSampNiche object..."))
    }
    MultiSampNiche_obj <- new("MultiSampNiche",
                              samp_info = samp_info,
                              niche_info = niche_info,
                              niche_cells = niche_cells,
                              niche_genes = niche_genes
    )
    STID_obj@STID_analysis@MultiSampNiche[[multi_id]] <- MultiSampNiche_obj
  }else{
    if(!is_empty(niche_cells)){
      clog_normal(paste0("MultiSampNiche object already exists for samp: ", multi_id))
      STID_niche_info <- STID_obj@STID_analysis@MultiSampNiche[[multi_id]]@niche_info
      if(niche_key %in% names(STID_niche_info)){
        clog_warn(paste0("niche_key: ", niche_key, " already exists, overwriting..."))
      }else{
        clog_normal(paste0("Adding niche_key: ", niche_key))
      }
      STID_obj@STID_analysis@MultiSampNiche[[multi_id]]@niche_info[[niche_key]] <- niche_info # is list, not data.frame
      STID_obj@STID_analysis@MultiSampNiche[[multi_id]]@niche_cells[[niche_key]] <- niche_cells
      STID_obj@STID_analysis@MultiSampNiche[[multi_id]]@niche_genes[[niche_key]] <- niche_genes
    }
  }

  # >>> Final
  clog_end()
  return(STID_obj)
}

#' @rdname CreateMultiSampNiche
#' @export
CreateMSNiche <- function(...){
  CreateMultiSampNiche(...)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Info getter and setter
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get information from STID object
#'
#' Retrieves specific information from the STID object's STID_info slot.
#'
#' @param STID_obj An STID object
#' @param info_key Character, main information category (e.g., "samp_info", "data_info")
#' @param sub_key Character vector, specific sub-keys to retrieve (optional)
#' @param ... Additional arguments passed to methods
#'
#' @return Requested information from the STID object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get sample column name
#' samp_colnm <- GetInfo(ist_obj, info_key = "data_info", sub_key = "samp_colnm")
#'
#' # Get all data information
#' data_info <- GetInfo(ist_obj, info_key = "data_info")
#' }
#'
GetInfo <- function(
    STID_obj,
    info_key,
    ...
) {
  UseMethod("GetInfo", STID_obj)
}

#' @rdname GetInfo
#' @export
GetInfo.STID <- function(
    STID_obj = NULL,
    info_key = NULL, # vector
    sub_key = NULL,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(info_key)
  if(length(info_key) > 1){
    clog_error("The length of info_key is greater than 1")
  }

  # >>> Start main pipeline
  # clog_step("Getting info...")
  info_keys <- names(STID_obj@STID_analysis@STID_info)
  res_list <- list()
  if(!info_key %in% info_keys){
    clog_error(paste0("info_key ", info_key, " not found in STID_info, skipping..."))
  }
  # clog_normal(paste0("Getting info for meta_key: ", info_key))
  if(is.null(sub_key)){
    res_list <- STID_obj@STID_analysis@STID_info[[info_key]]
  }else{
    if(!all(sub_key %in% names(STID_obj@STID_analysis@STID_info[[info_key]]))){
      clog_error(paste0("sub_key not found in STID_info$", info_key, ": ",
                        paste(setdiff(sub_key, names(STID_obj@STID_analysis@STID_info[[info_key]])), collapse = ", "),
                        ", skipping..."))
    }
    # clog_normal(paste0("Getting info for sub_key: ", paste(sub_key, collapse = ", ")))
    res_list <- STID_obj@STID_analysis@STID_info[[info_key]][sub_key]
  }
  return(res_list)
}

#' Get coordinate-related column names
#'
#' Convenience function to retrieve all coordinate-related column names
#' from an STID object.
#'
#' @param STID_obj An STID object
#'
#' @return Named vector of coordinate column names
#'
#' @export
#'
#' @examples
#' \dontrun{
#' coord_info <- GetCoordInfo(ist_obj)
#' }
#'
GetCoordInfo <- function(
    STID_obj
) {
  UseMethod("GetCoordInfo", STID_obj)
}

#' @rdname GetCoordInfo
#' @export
GetCoordInfo.STID <- function(STID_obj = NULL){
  coord_list <- GetInfo(
    STID_obj = STID_obj,
    info_key = "data_info",
    sub_key = c("samp_colnm", "samp_grp_colnm","x_colnm","y_colnm","celltype_colnm")
  )
  coord_list <- lapply(coord_list, .null_to_na)
  coord_nm <- coord_list %>% unlist(use.names = T) # must be unlist, not be [[1]]
  return(coord_nm)
}


#' Set information in STID object
#'
#' Modifies specific information in the STID object's STID_info slot.
#'
#' @param STID_obj An STID object
#' @param info_key Character, main information category
#' @param sub_key Character, specific sub-key to set (optional)
#' @param info_value Value to set
#' @param ... Additional arguments passed to methods
#'
#' @return Modified STID object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Update project description
#' STID_obj <- SetInfo(ist_obj,
#'                    info_key = "project_info",
#'                    sub_key = "description",
#'                    info_value = "Updated project description")
#' }
#'
SetInfo <- function(
    STID_obj,
    info_key,
    sub_key,
    info_value,
    ...
) {
  UseMethod("SetInfo", STID_obj)
}

#' @rdname SetInfo
#' @export
SetInfo.STID <- function(
    STID_obj = NULL,
    info_key = NULL, # string
    sub_key = NULL,
    info_value = NULL, # any
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(info_key,info_value)

  # >>> Start main pipeline
  clog_step("Setting info...")
  info_keys <- names(STID_obj@STID_analysis@STID_info)
  if(!info_key %in% info_keys){
    clog_error(paste0("info_key ", info_key, " not found in STID_info"))
  }
  clog_normal(paste0("Setting info for meta_key: ", info_key))
  if(is.null(sub_key)){
    STID_obj@STID_analysis@STID_info[[info_key]] <- info_value
  }else{
    clog_normal(paste0("Setting info for sub_key: ", sub_key))
    STID_obj@STID_analysis@STID_info[[info_key]][[sub_key]] <- info_value
  }
  return(STID_obj)
}

#' Add information to STID object
#'
#' Appends information to existing slots in the STID_info, primarily for
#' adding comments and notes.
#'
#' @param STID_obj An STID object
#' @param info_key Character, information category (default: "comment_info")
#' @param sub_key Character, specific sub-key to add to
#' @param info_value Value to append
#' @param ... Additional arguments passed to methods
#'
#' @return Modified STID object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Add a comment about analysis parameters
#' STID_obj <- AddInfo(ist_obj,
#'                    info_key = "comment_info",
#'                    sub_key = "parameters",
#'                    info_value = "Used DBSCAN with eps=5, minPts=10")
#' }
#'
AddInfo <- function(
    STID_obj,
    info_key,
    sub_key,
    info_value,
    ...
) {
  UseMethod("AddInfo", STID_obj)
}

#' @rdname AddInfo
#' @export
AddInfo.STID <- function(
    STID_obj = NULL,
    info_key = "comment_info", # string
    sub_key = NULL,
    info_value = NULL, # list or vector
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(info_key,info_value)
  if(info_key != "comment_info"){
    clog_warn("Adding info is only recommended for 'comment_info' meta_key")
  }

  # >>> Start main pipeline
  clog_step("Adding info...")
  info_keys <- names(STID_obj@STID_analysis@STID_info)
  if(!info_key %in% info_keys){
    clog_error(paste0("info_key ", info_key, " not found in STID_info"))
  }
  clog_normal(paste0("Adding info for meta_key: ", info_key))
  info_list <- STID_obj@STID_analysis@STID_info[[info_key]]
  if(is.null(sub_key)){
    STID_obj@STID_analysis@STID_info[[info_key]] <- c(info_list, info_value)
  }else{
    if(!sub_key %in% names(info_list)){
      clog_error(paste0("sub_key ", sub_key, " not found in STID_info$", info_key))
    }
    clog_normal(paste0("Adding info for sub_key: ", sub_key))
    STID_obj@STID_analysis@STID_info[[info_key]][[sub_key]] <- c(info_list[[sub_key]], info_value)
  }
  return(STID_obj)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MetaData getter and setter
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get metadata from STID object
#'
#' Retrieves metadata from the STID object's meta_data_record, optionally
#' combining multiple metadata sources and adding coordinate information.
#'
#' @param STID_obj An STID object
#' @param meta_key Character vector, metadata keys to retrieve
#' @param add_coord Logical, whether to add coordinate columns (default: TRUE)
#' @param ... Additional arguments passed to methods
#'
#' @return List of metadata data frames, one for each meta_key
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get raw metadata
#' raw_meta <- GetMetaData(ist_obj, meta_key = "raw")
#'
#' # Get multiple metadata sources combined
#' combined_meta <- GetMetaData(ist_obj,
#'                              meta_key = list(c("raw", "coord"),
#'                                             c("M1_SpotDetect_Gene_20240101")))
#' }
#'
GetMetaData <- function(
    STID_obj,
    meta_key,
    ...
) {
  UseMethod("GetMetaData", STID_obj)
}

#' @rdname GetMetaData
#' @export
GetMetaData.STID <- function(
    STID_obj = NULL,
    meta_key = NULL, # vector, more meta data
    add_coord = TRUE,
    ...
) {

  # >>> Start pipeline
  # clog_step("Getting meta data...")
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(meta_key)

  # >>> Start main pipeline
  meta_keys <- names(STID_obj@STID_analysis@meta_data_record$meta_data_list)
  res_list <- list()
  coord_colnm <- GetCoordInfo(STID_obj) %>% na.omit()
  coord_data <- STID_obj@meta.data[coord_colnm]

  # >
  meta_key <- meta_key %>% as.list()
  for (i in seq_along(meta_key)){
    # clog_loop(paste0("Processing meta_key: ", i, "/", length(meta_key)))
    i_key <- meta_key[[i]]
    res_df_list <- list()
    for(j in seq_along(i_key)){
      j_key <- i_key[[j]]
      if( j_key == "raw"){
        # clog_normal(paste0("Getting meta data for meta_key: ", j_key))
        res_df_list[[j_key]] <- STID_obj@meta.data
      }else if(j_key == "coord"){
        # clog_normal(paste0("Getting meta data for meta_key: ", j_key))
        res_df_list[[j_key]] <- coord_data
      }else{
        if(!j_key %in% meta_keys){
          clog_error(paste0("meta_key ", j_key, " not found in STID meta_data_record, skipping..."))
          next
        }
        # clog_normal(paste0("Getting meta data for meta_key: ", j_key))
        res_df_list[[j_key]] <- STID_obj@STID_analysis@meta_data_record$meta_data_list[[j_key]]
      }
    }
    # check rownames
    if(length(unique(sapply(res_df_list, function(df) paste(rownames(df), collapse = "\r")))) > 1){
      clog_error(paste0("Meta data rownames do not match for meta_key: ", paste(i_key, collapse = ", ")))
    }
    # res_df <- bind_cols(res_df_list)
    res_df <- .bind_cols_keep_first(res_df_list)

    if(add_coord){
      select_colnm <- setdiff(coord_colnm, colnames(res_df))
      if(length(select_colnm) == 0){
        clog_warn(paste0("Meta data already contains coordinate columns, skipping adding coordinates..."))
        res_list[[i]] <- res_df
        next
      }
      # clog_normal(paste0("Adding coordinate columns to meta data : ",paste0(select_colnm, collapse = ", ")))
      # res_list[[i]] <- bind_cols(coord_data[select_colnm],res_df,.name_repair = "unique") # trigger new names
      res_list[[i]] <- cbind(coord_data[select_colnm],res_df)
    }else{
      res_list[[i]] <- res_df
    }
  }
  return(res_list)
}

#' Bind data frames while keeping first occurrence of duplicate columns
#'
#' Internal function to combine multiple data frames, skipping columns that
#' already exist in the result.
#'
#' @param dfs List of data frames to bind
#'
#' @return Combined data frame with unique column names
#'
#' @keywords internal
#'
#' @noRd
#'
.bind_cols_keep_first <- function(dfs = NULL) {
  if (length(dfs) == 0) {
    clog_error("No data frames provided to bind_cols_keep_first.")
  }
  if (length(dfs) == 1) {
    return(bind_cols(dfs))
    }
  result <- dfs[[1]]
  skipped_cols <- character()
  for (i in 2:length(dfs)) {
    df_i <- dfs[[i]]
    if (is.null(df_i) || ncol(df_i) == 0) {
      # clog_warn(paste0("Data frame ", i, " is NULL or has no columns, skipping..."))
      next
    }

    existing_names <- names(result)
    new_cols <- setdiff(names(df_i), existing_names)
    dup_cols  <- intersect(names(df_i), existing_names)
    if (length(dup_cols) > 0) {
      skipped_cols <- c(skipped_cols, dup_cols)
    }
    if (length(new_cols) > 0) {
      result <- dplyr::bind_cols(result, df_i[new_cols])
    }
  }
  skipped_cols <- unique(skipped_cols)
  if (length(skipped_cols) > 0) {
    # clog_warn(paste0(
    #   "The following columns were skipped (already present in earlier data frames): ",
    #   paste(sQuote(skipped_cols), collapse = ", ")
    # ))
  }
  return(result)
}

#' Get metadata information
#'
#' Retrieves the metadata tracking information data frame from an STID object.
#'
#' @param STID_obj An STID object
#' @param ... Additional arguments passed to methods
#'
#' @return Data frame containing metadata tracking information
#'
#' @export
#'
#' @examples
#' \dontrun{
#' meta_info <- GetMetaInfo(ist_obj)
#' }
#'
GetMetaInfo <- function(
    STID_obj,
    ...
) {
  UseMethod("GetMetaInfo", STID_obj)
}

#' @rdname GetMetaInfo
#' @export
GetMetaInfo.STID <- function(
    STID_obj = NULL,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")

  # >>> Start main pipeline
  meta_data_info <- STID_obj@STID_analysis@meta_data_record$meta_data_info
  return(meta_data_info)
}


#' Add metadata to STID object
#'
#' Adds a new metadata data frame to the STID object's meta_data_record and
#' updates the tracking information.
#'
#' @param STID_obj An STID object
#' @param meta_key Character, unique identifier for the metadata
#' @param add_data Data frame, metadata to add
#' @param dir_nm Character, directory name for output (optional)
#' @param grp_nm Character, group name for output (optional)
#' @param asso_key Character, associated metadata key (optional)
#' @param description Character, description of the metadata (optional)
#' @param ... Additional arguments passed to methods
#'
#' @return Modified STID object with added metadata
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Add custom metadata
#' STID_obj <- AddMetaData(ist_obj,
#'                        meta_key = "custom_annotation",
#'                        add_data = custom_metadata,
#'                        description = "Manual cell annotations")
#' }
#'
AddMetaData <- function(
    STID_obj,
    ...
) {
  UseMethod("AddMetaData", STID_obj)
}

#' @rdname AddMetaData
#' @export
AddMetaData.STID <- function(
    STID_obj = NULL,
    meta_key = NULL, # string
    add_data = NULL, # data.frame
    dir_nm = NA,
    grp_nm = NA,
    asso_key = NULL,
    description = NULL,
    ...
) {

  # >>> Check input patameter
  clog_step("Adding meta data...")
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(meta_key,add_data)
  clog_normal(paste0("Your meta_key: ", meta_key))
  # >>> End check

  # >>> Start main pipeline
  asso_key <- asso_key %>% .null_to_na() %>% as.list()
  STID_obj <- .AddMetaData_Info(STID_obj = STID_obj,
                               meta_key = meta_key,
                               add_data = add_data,
                               dir_nm = dir_nm,
                               grp_nm = grp_nm,
                               asso_key = asso_key,
                               description = description
  )
  STID_obj <- .AddMetaData_List(STID_obj = STID_obj,
                                  meta_key = meta_key,
                                  add_data = add_data
  )
  return(STID_obj)
}

#' Add metadata to STID object info record
#'
#' Internal function to update the metadata tracking information.
#'
#' @param STID_obj An STID object
#' @param meta_key Character, metadata key
#' @param add_data Data frame, metadata to add
#' @param asso_key Character, associated metadata key
#' @param description Character, description
#' @param dir_nm Character, directory name
#' @param grp_nm Character, group name
#'
#' @return Modified STID object
#'
#' @keywords internal
#'
#' @noRd
#'
.AddMetaData_Info <- function(STID_obj = NULL,meta_key = NULL, add_data = NULL,
                              asso_key = NULL, description = NULL,
                              dir_nm = NA, grp_nm = NA
){
  now_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  top_function <- sys.calls()[[1]][[1]] %>% deparse()
  meta_data_info <- STID_obj@STID_analysis@meta_data_record$meta_data_info
  add_metadata_info <- data.frame(
    meta_key = meta_key,
    time = now_time,
    func_nm = top_function,
    dir_nm = dir_nm,
    grp_nm = grp_nm,
    asso_key = I(asso_key), # keep list structure
    description = description %>% .null_to_na(),
    stringsAsFactors = FALSE,
    row.names = meta_key
  )
  if(meta_key %in% rownames(meta_data_info)){
    clog_warn("The current meta_key already exists in meta_data_info, overwrite it.")
    meta_data_info[meta_key, ] <- add_metadata_info
  }else{
    meta_data_info <- rbind(meta_data_info, add_metadata_info)
  }
  STID_obj@STID_analysis@meta_data_record$meta_data_info <- meta_data_info
  return(STID_obj)
}

#' Add metadata to STID object list record
#'
#' Internal function to store metadata in the meta_data_list.
#'
#' @param STID_obj An STID object
#' @param meta_key Character, metadata key
#' @param add_data Data frame, metadata to add
#'
#' @return Modified STID object
#'
#' @keywords internal
#'
#' @noRd
#'
.AddMetaData_List <- function(STID_obj = NULL, meta_key = NULL, add_data = NULL){
  if(all(rownames(add_data) != rownames(STID_obj@meta.data))){
    clog_error("Row names of add_data must match row names of STID_obj meta.data")
  }
  STID_obj@STID_analysis@meta_data_record$meta_data_list[[meta_key]] <- add_data
  return(STID_obj)
}

#' Remove metadata from STID object
#'
#' Removes specified metadata from the STID object's meta_data_record.
#'
#' @param STID_obj An STID object
#' @param meta_key Character vector, metadata keys to remove
#' @param ... Additional arguments passed to methods
#'
#' @return Modified STID object with removed metadata
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Remove temporary metadata
#' STID_obj <- RemoveMetaData(ist_obj, meta_key = "temp_analysis")
#' }
#'
RemoveMetaData <- function(
    STID_obj,
    meta_key,
    ...
) {
  UseMethod("RemoveMetaData", STID_obj)
}

#' @rdname RemoveMetaData
#' @export
RemoveMetaData.STID <- function(
    STID_obj = NULL,
    meta_key = NULL, # vector
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(meta_key)
  clog_normal(paste0("Your meta_key: ", paste0(meta_key, collapse = ", ")))

  # >>> Start main pipeline
  clog_step("Removing meta data...")
  meta_data_info <- STID_obj@STID_analysis@meta_data_record$meta_data_info
  for (i_key in meta_key){
    if(!i_key %in% meta_data_info$meta_key){
      clog_warn(paste0("meta_key ", i_key, " not found in STID meta_data_record, skipping..."))
      next
    }
    clog_normal(paste0("Removing meta data for meta_key: ", i_key))
    meta_data_info <- meta_data_info %>% filter(meta_key != i_key)
    STID_obj@STID_analysis@meta_data_record$meta_data_list[[i_key]] <- NULL
  }
  STID_obj@STID_analysis@meta_data_record$meta_data_info <- meta_data_info
  return(STID_obj)
}

#' Add columns to existing metadata
#'
#' Adds new columns to an existing metadata data frame in the STID object.
#'
#' @param STID_obj An STID object
#' @param meta_key Character, metadata key to modify
#' @param add_data Data frame, new columns to add
#' @param ignore_rownm Logical, whether to ignore row name matching (default: FALSE)
#' @param ... Additional arguments passed to methods
#'
#' @return Modified STID object with added metadata columns
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Add new columns to raw metadata
#' new_cols <- data.frame(new_score = runif(ncol(ist_obj)))
#' STID_obj <- AddMetaColumn(ist_obj,
#'                          meta_key = "raw",
#'                          add_data = new_cols)
#' }
#'
AddMetaColumn <- function(
    STID_obj,
    ...
) {
  UseMethod("AddMetaColumn", STID_obj)
}

#' @rdname AddMetaColumn
#' @export
AddMetaColumn.STID <- function(
    STID_obj = NULL,
    meta_key = NULL, # string
    add_data = NULL, # data.frame
    ignore_rownm = FALSE,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(add_data,meta_key)

  # >>> Start main pipeline
  clog_step("Adding meta columns...")
  meta_keys <- names(STID_obj@STID_analysis@meta_data_record$meta_data_list)
  if(meta_key == "raw"){
    meta_data <- GetMetaData(STID_obj, meta_key = "raw")[[1]]
  }else{
    if(!meta_key %in% meta_keys){
      clog_error(paste0("meta_key ", meta_key, " not found in STID meta_data_record"))
    }
    meta_data <- GetMetaData(STID_obj, meta_key = meta_key)[[1]]
  }
  if(!ignore_rownm){
    if(!all(rownames(add_data) == rownames(meta_data))){
      clog_error("Row names of add_data must match row names of existing meta data")
    }
  }else{
    if(nrow(add_data) != nrow(meta_data)){
      clog_error("Number of rows in add_data must match number of rows in existing meta data")
    }
  }
  new_meta_data <- bind_cols(meta_data, add_data)
  if(meta_key == "raw"){
    STID_obj@meta.data <- new_meta_data # AddMetaData???
  }else{
    STID_obj@STID_analysis@meta_data_record$meta_data_list[[meta_key]] <- new_meta_data
  }
  return(STID_obj)
}

#' Remove columns from existing metadata
#'
#' Removes specified columns from an existing metadata data frame in the STID object.
#'
#' @param STID_obj An STID object
#' @param meta_key Character, metadata key to modify
#' @param remove_colnm Character vector, column names to remove
#' @param ... Additional arguments passed to methods
#'
#' @return Modified STID object with removed metadata columns
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Remove temporary columns
#' STID_obj <- RemoveMetaColumn(ist_obj,
#'                             meta_key = "raw",
#'                             remove_colnm = c("temp_score", "temp_label"))
#' }
#'
RemoveMetaColumn <- function(
    STID_obj,
    ...
) {
  UseMethod("RemoveMetaColumn", STID_obj)
}

#' @rdname RemoveMetaColumn
#' @export
RemoveMetaColumn.STID <- function(
    STID_obj = NULL,
    meta_key = NULL, # string
    remove_colnm = NULL, # vector
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(meta_key,remove_colnm)

  # >>> Start main pipeline
  clog_step("Removing meta columns...")
  meta_keys <- names(STID_obj@STID_analysis@meta_data_record$meta_data_list)
  if(meta_key == "raw"){
    meta_data <- GetMetaData(STID_obj, meta_key = "raw")[[1]]
  }else{
    if(!meta_key %in% meta_keys){
      clog_error(paste0("meta_key ", meta_key, " not found in STID meta_data_record"))
    }
    meta_data <- GetMetaData(STID_obj, meta_key = meta_key)[[1]]
  }
  if(!all(remove_colnm %in% colnames(meta_data))){
    clog_error(paste0("Some remove_colnm not found in existing meta data: ",
                      paste(setdiff(remove_colnm, colnames(meta_data)), collapse = ", ")))
  }
  new_meta_data <- meta_data %>% dplyr::select(-all_of(remove_colnm))
  if(meta_key == "raw"){
    STID_obj@meta.data <- new_meta_data
  }else{
    STID_obj <- .AddMetaData_List(STID_obj = STID_obj,
                                    meta_key = meta_key,
                                    add_data = new_meta_data
    )
  }
  return(STID_obj)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SingSampNiche Info getter and setter
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get SingleSampNiche information
#'
#' Retrieves niche information from SingleSampNiche objects for specified samples.
#'
#' @param STID_obj An STID object
#' @param loop_id Character vector, sample identifiers
#' @param niche_key Character, niche key to retrieve
#' @param ... Additional arguments passed to methods
#'
#' @return List of niche information data frames, one per sample
#'
#' @export
#'
#' @examples
#' \dontrun{
#' niche_info <- GetSSNicheInfo(ist_obj,
#'                              loop_id = c("sample1", "sample2"),
#'                              niche_key = "niche_virulence")
#' }
#'
GetSSNicheInfo <- function(
    STID_obj,
    loop_id,
    niche_key,
    ...
) {
  UseMethod("GetSSNicheInfo", STID_obj)
}

#' @rdname GetSSNicheInfo
#' @export
GetSSNicheInfo.STID <- function(
    STID_obj = NULL,
    loop_id = "LoopAllSamp", # vector
    niche_key = NULL,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(loop_id,niche_key)
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)

  # >>> Start main pipeline
  clog_normal("Getting SingleSampNiche info...")
  niche_info_list <- list()
  for(i_single in loop_single){
    logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
    if(!logic_niche){
      clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
      next
    }
    niche_info <- STID_obj@STID_analysis@SingleSampNiche[[i_single]]@niche_info %>%
      filter(niche_key == niche_key)
    niche_info_list[[i_single]] <- niche_info
  }
  return(niche_info_list)
}

#' Get SingleSampNiche cells
#'
#' Retrieves cell-level data from SingleSampNiche objects for specified samples.
#'
#' @param STID_obj An STID object
#' @param loop_id Character vector, sample identifiers
#' @param niche_key Character, niche key to retrieve
#' @param ... Additional arguments passed to methods
#'
#' @return List of cell data frames, one per sample
#'
#' @export
#'
#' @examples
#' \dontrun{
#' niche_cells <- GetSSNicheCells(ist_obj,
#'                                loop_id = c("sample1", "sample2"),
#'                                niche_key = "niche_virulence")
#' }
#'
GetSSNicheCells <- function(
    STID_obj,
    loop_id,
    niche_key,
    ...
) {
  UseMethod("GetSSNicheCells", STID_obj)
}

#' @rdname GetSSNicheCells
#' @export
GetSSNicheCells.STID <- function(
    STID_obj = NULL,
    loop_id = "LoopAllSamp", # vector
    niche_key = NULL,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(loop_id,niche_key)
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)

  # >>> Start main pipeline
  clog_normal("Getting SingleSampNiche cells...")
  niche_cells_list <- list()
  for(i_single in loop_single){
    logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
    if(!logic_niche){
      clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
      next
    }
    niche_cells <- STID_obj@STID_analysis@SingleSampNiche[[i_single]]@niche_cells[[niche_key]]
    niche_cells_list[[i_single]] <- niche_cells
  }
  return(niche_cells_list)
}

#' Add columns to SingleSampNiche cells
#'
#' Adds new columns to the cell data in SingleSampNiche objects.
#'
#' @param STID_obj An STID object
#' @param loop_id Character vector, sample identifiers
#' @param meta_key Character, metadata key containing new data
#' @param select_colnm Character vector, columns to add
#' @param niche_key Character, niche key to modify
#' @param ... Additional arguments passed to methods
#'
#' @return Modified STID object with updated niche cell data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Add new annotation columns to niche cells
#' STID_obj <- AddSSNicheCells(ist_obj,
#'                            loop_id = "LoopAllSamp",
#'                            meta_key = "custom_annotation",
#'                            select_colnm = c("cell_type", "confidence"),
#'                            niche_key = "niche_virulence")
#' }
#'
AddSSNicheCells <- function(
    STID_obj,
    loop_id,
    meta_key,
    select_colnm,
    niche_key,
    ...
) {
  UseMethod("AddSSNicheCells", STID_obj)
}

#' @rdname AddSSNicheCells
#' @export
AddSSNicheCells.STID <- function(
    STID_obj = NULL,
    loop_id = "LoopAllSamp", # vector
    meta_key = NULL, # string
    select_colnm = NULL, # vector
    niche_key = NULL,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(loop_id,meta_key,select_colnm,niche_key)
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)

  # >>> Start main pipeline
  clog_normal("Adding SingleSampNiche cells...")
  if(meta_key == "raw"){
    meta_data <-STID_obj@meta.data
  }else if(meta_key == "coord"){
    meta_data <- GetMetaData(STID_obj, meta_key = "coord")[[1]]
  }else if(!meta_key %in% names(STID_obj@STID_analysis@meta_data_record$meta_data_list)){
    clog_error(paste0("meta_key ", meta_key, " not found in STID meta_data_record"))
  }else{
    meta_data <- STID_obj@STID_analysis@meta_data_record$meta_data_list[[meta_key]]
    # meta_data <- GetMetaData(STID_obj, meta_key = meta_key)[[1]] # the same
  }
  .check_column_exist(meta_data,select_colnm)
  new_niche_cells <- meta_data %>%
    dplyr::select(all_of(select_colnm))
  for(i_single in loop_single){
    logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
    if(!logic_niche){
      clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
      next
    }
    raw_niche_cells <- STID_obj@STID_analysis@SingleSampNiche[[i_single]]@niche_cells[[niche_key]]
    i_new_niche_cells <- new_niche_cells[rownames(raw_niche_cells), , drop = FALSE]
    niche_cells <- bind_cols(raw_niche_cells, i_new_niche_cells)
    STID_obj@STID_analysis@SingleSampNiche[[i_single]]@niche_cells[[niche_key]] <- niche_cells
  }
  return(STID_obj)
}

#' Get SingleSampNiche genes
#'
#' Retrieves gene-level data from SingleSampNiche objects for specified samples.
#'
#' @param STID_obj An STID object
#' @param loop_id Character vector, sample identifiers
#' @param niche_key Character, niche key to retrieve
#' @param ... Additional arguments passed to methods
#'
#' @return List of gene data frames, one per sample
#'
#' @export
#'
#' @examples
#' \dontrun{
#' niche_genes <- GetSSNicheGenes(ist_obj,
#'                                loop_id = c("sample1", "sample2"),
#'                                niche_key = "niche_virulence")
#' }
#'
GetSSNicheGenes <- function(
    STID_obj,
    loop_id,
    niche_key,
    ...
) {
  UseMethod("GetSSNicheGenes", STID_obj)
}

#' @rdname GetSSNicheGenes
#' @export
GetSSNicheGenes.STID <- function(
    STID_obj = NULL,
    loop_id = "LoopAllSamp", # vector
    niche_key = NULL,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(loop_id,niche_key)
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)

  # >>> Start main pipeline
  clog_normal("Getting SingleSampNiche genes...")
  niche_genes_list <- list()
  for(i_single in loop_single){
    logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
    if(!logic_niche){
      clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
      next
    }
    niche_genes <- STID_obj@STID_analysis@SingleSampNiche[[i_single]]@niche_genes[[niche_key]]
    niche_genes_list[[i_single]] <- niche_genes
  }
  return(niche_genes_list)
}

#' Add annotations to SingleSampNiche genes
#'
#' Adds new columns to the gene data in SingleSampNiche objects.
#'
#' @param STID_obj An STID object
#' @param gene Character vector, gene names to annotate
#' @param label Character vector, labels to assign to genes (optional)
#' @param add_colnm Character, name of new column to add
#' @param loop_id Character vector, sample identifiers
#' @param niche_key Character, niche key to modify
#' @param ... Additional arguments passed to methods
#'
#' @return Modified STID object with updated niche gene data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Add gene family annotations
#' STID_obj <- AddSSNicheGenes(ist_obj,
#'                            gene = c("gene1", "gene2", "gene3"),
#'                            label = c("toxin", "toxin", "adhesin"),
#'                            add_colnm = "gene_family",
#'                            niche_key = "niche_virulence")
#' }
#'
AddSSNicheGenes <- function(
    STID_obj,
    ...
) {
  UseMethod("AddSSNicheGenes", STID_obj)
}

#' @rdname AddSSNicheGenes
#' @export
AddSSNicheGenes.STID <- function(
    STID_obj = NULL,
    gene = NULL,
    label = NULL,
    add_colnm = NULL,
    loop_id = "LoopAllSamp", # vector
    niche_key = NULL,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(gene,add_colnm,loop_id,niche_key)
  if(sum(duplicated(gene)) > 0){
    clog_error("gene contains duplicated gene names, please check.")
  }
  if(!all(gene %in% rownames(STID_obj))){
    clog_error(paste0("Some gene not found in STID_obj rownames: ",
                      paste(setdiff(gene, rownames(STID_obj)), collapse = ", ")))
  }
  loop_single <- .check_loop_single(STID_obj = STID_obj, loop_id = loop_id)

  # >>> Start main pipeline
  clog_normal("Adding SingleSampNiche genes...")
  for(i_single in loop_single){
    logic_niche <- .check_niche_exist(STID_obj, i_single, niche_key)
    if(!logic_niche){
      clog_warn(paste0("niche_key ", niche_key, " not found in SingleSampNiche of samp_id ", i_single, ", skipping..."))
      next
    }
    niche_genes <- GetSSNicheGenes(STID_obj, loop_id = i_single, niche_key = niche_key)[[1]]
    if(add_colnm %in% colnames(niche_genes)){
      clog_warn(paste0("add_colnm ", add_colnm, " already exists in niche_genes of samp_id ", i_single, ", overwrite it."))
    }
    niche_genes[[add_colnm]] <- NA
    if(is.null(label)){
      niche_genes[gene,add_colnm] <- gene
    }else{
      niche_genes[gene,add_colnm]  <- label
    }
    STID_obj@STID_analysis@SingleSampNiche[[i_single]]@niche_genes[[niche_key]] <- niche_genes
  }
  return(STID_obj)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MultiSampNiche Info getter and setter
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get MultiSampNiche information
#'
#' Retrieves niche information from MultiSampNiche objects for specified multi-sample analyses.
#'
#' @param STID_obj An STID object
#' @param loop_id Character vector, multi-sample analysis identifiers
#' @param niche_key Character, niche key to retrieve
#' @param ... Additional arguments passed to methods
#'
#' @return List of niche information data frames, one per multi-sample analysis
#'
#' @export
#'
#' @examples
#' \dontrun{
#' multi_info <- GetMSNicheInfo(ist_obj,
#'                              loop_id = c("comparison1", "comparison2"),
#'                              niche_key = "niche_virulence")
#' }
#'
GetMSNicheInfo <- function(
    STID_obj,
    loop_id,
    niche_key,
    ...
) {
  UseMethod("GetMSNicheInfo", STID_obj)
}

#' @rdname GetMSNicheInfo
#' @export
GetMSNicheInfo.STID <- function(
    STID_obj = NULL,
    loop_id = "LoopAllMulti", # vector
    niche_key = NULL,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(loop_id,niche_key)
  loop_multi <- .check_loop_multi(STID_obj = STID_obj, loop_id = loop_id)

  # >>> Start main pipeline
  clog_normal("Getting MultiSampNiche info...")
  niche_info_list <- list()
  for(i_multi in loop_multi){
    logic_niche <- .check_niche_exist(STID_obj, i_multi, niche_key, samp_mode  = "MS")
    if(!logic_niche){
      clog_warn(paste0("niche_key: ", niche_key, " not found in MultiSampNiche of multi_id: ", i_multi, ", skipping..."))
      next
    }
    niche_info <- STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@niche_info[[niche_key]]
    niche_info_list[[i_multi]] <- niche_info
  }
  return(niche_info_list)
}

#' Get MultiSampNiche cells
#'
#' Retrieves cell-level data from MultiSampNiche objects for specified multi-sample analyses.
#'
#' @param STID_obj An STID object
#' @param loop_id Character vector, multi-sample analysis identifiers
#' @param niche_key Character, niche key to retrieve
#' @param ... Additional arguments passed to methods
#'
#' @return List of cell data frames, one per multi-sample analysis
#'
#' @export
#'
#' @examples
#' \dontrun{
#' multi_cells <- GetMSNicheCells(ist_obj,
#'                                loop_id = c("comparison1", "comparison2"),
#'                                niche_key = "niche_virulence")
#' }
#'
GetMSNicheCells <- function(
    STID_obj,
    loop_id,
    niche_key,
    ...
) {
  UseMethod("GetMSNicheCells", STID_obj)
}


#' @rdname GetMSNicheCells
#' @export
GetMSNicheCells.STID <- function(
    STID_obj = NULL,
    loop_id = "LoopAllMulti", # vector
    niche_key = NULL,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(loop_id,niche_key)
  loop_multi <- .check_loop_multi(STID_obj = STID_obj, loop_id = loop_id)

  # >>> Start main pipeline
  clog_normal("Getting MultiSampNiche cells...")
  niche_cells_list <- list()
  for(i_multi in loop_multi){
    logic_niche <- .check_niche_exist(STID_obj, i_multi, niche_key, samp_mode  = "MS")
    if(!logic_niche){
      clog_warn(paste0("niche_key: ", niche_key, " not found in MultiSampNiche of multi_id: ", i_multi, ", skipping..."))
      next
    }
    niche_cells <- STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@niche_cells[[niche_key]]
    niche_cells_list[[i_multi]] <- niche_cells
  }
  return(niche_cells_list)
}

#' Add columns to MultiSampNiche cells
#'
#' Adds new columns to the cell data in MultiSampNiche objects.
#'
#' @param STID_obj An STID object
#' @param loop_id Character vector, multi-sample analysis identifiers
#' @param meta_key Character, metadata key containing new data
#' @param select_colnm Character vector, columns to add
#' @param niche_key Character, niche key to modify
#' @param ... Additional arguments passed to methods
#'
#' @return Modified STID object with updated multi-sample niche cell data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Add new annotations to multi-sample niche cells
#' STID_obj <- AddMSNicheCells(ist_obj,
#'                            loop_id = "LoopAllMulti",
#'                            meta_key = "custom_annotation",
#'                            select_colnm = c("cell_type", "confidence"),
#'                            niche_key = "niche_virulence")
#' }
#'
AddMSNicheCells <- function(
    STID_obj,
    loop_id,
    meta_key,
    select_colnm,
    niche_key,
    ...
) {
  UseMethod("AddMSNicheCells", STID_obj)
}

#' @rdname AddMSNicheCells
#' @export
AddMSNicheCells.STID <- function(
    STID_obj = NULL,
    loop_id = "LoopAllMulti", # vector
    meta_key = NULL, # string
    select_colnm = NULL, # vector
    niche_key = NULL,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(loop_id,meta_key,select_colnm,niche_key)
  loop_multi <- .check_loop_multi(STID_obj = STID_obj, loop_id = loop_id)

  # >>> Start main pipeline
  clog_normal("Adding MultiSampNiche cells...")
  if(meta_key == "raw"){
    meta_data <-STID_obj@meta.data
  }else if(meta_key == "coord"){
    meta_data <- GetMetaData(STID_obj, meta_key = "coord")[[1]]
  }else if(!meta_key %in% names(STID_obj@STID_analysis@meta_data_record$meta_data_list)){
    clog_error(paste0("meta_key ", meta_key, " not found in STID meta_data_record"))
  }else{
    meta_data <- STID_obj@STID_analysis@meta_data_record$meta_data_list[[meta_key]]
    # meta_data <- GetMetaData(STID_obj, meta_key = meta_key)[[1]] # the same
  }
  .check_column_exist(meta_data,select_colnm)
  new_niche_cells <- meta_data %>%
    dplyr::select(all_of(select_colnm))
  for(i_multi in loop_multi){
    logic_niche <- .check_niche_exist(STID_obj, i_multi, niche_key, samp_mode  = "MS")
    if(!logic_niche){
      clog_warn(paste0("niche_key: ", niche_key, " not found in MultiSampNiche of multi_id: ", i_multi, ", skipping..."))
      next
    }
    raw_niche_cells <- STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@niche_cells[[niche_key]]
    i_new_niche_cells <- new_niche_cells[rownames(raw_niche_cells), , drop = FALSE]
    niche_cells <- bind_cols(raw_niche_cells, i_new_niche_cells)
    STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@niche_cells[[niche_key]] <- niche_cells
  }
  return(STID_obj)
}

#' Get MultiSampNiche genes
#'
#' Retrieves gene-level data from MultiSampNiche objects for specified multi-sample analyses.
#'
#' @param STID_obj An STID object
#' @param loop_id Character vector, multi-sample analysis identifiers
#' @param niche_key Character, niche key to retrieve
#' @param ... Additional arguments passed to methods
#'
#' @return List of gene data frames, one per multi-sample analysis
#'
#' @export
#'
#' @examples
#' \dontrun{
#' multi_genes <- GetMSNicheGenes(ist_obj,
#'                                loop_id = c("comparison1", "comparison2"),
#'                                niche_key = "niche_virulence")
#' }
#'
GetMSNicheGenes <- function(
    STID_obj,
    loop_id,
    niche_key,
    ...
) {
  UseMethod("GetMSNicheGenes", STID_obj)
}

#' @rdname GetMSNicheGenes
#' @export
GetMSNicheGenes.STID <- function(
    STID_obj = NULL,
    loop_id = "LoopAllMulti", # vector
    niche_key = NULL,
    ...
) {

  # >>> Start pipeline
  if (!is(STID_obj, "STID")) clog_error("STID_obj must be an STID STID_obj")
  .check_not_any_null(loop_id,niche_key)
  loop_multi <- .check_loop_multi(STID_obj = STID_obj, loop_id = loop_id)

  # >>> Start main pipeline
  clog_normal("Getting MultiSampNiche genes...")
  niche_genes_list <- list()
  for(i_multi in loop_multi){
    logic_niche <- .check_niche_exist(STID_obj, i_multi, niche_key, samp_mode  = "MS")
    if(!logic_niche){
      clog_warn(paste0("niche_key: ", niche_key, " not found in MultiSampNiche of multi_id: ", i_multi, ", skipping..."))
      next
    }
    niche_genes <- STID_obj@STID_analysis@MultiSampNiche[[i_multi]]@niche_genes[[niche_key]]
    niche_genes_list[[i_multi]] <- niche_genes
  }
  return(niche_genes_list)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# show method: STID object overview
# 这里就用原始的数据结构吧，不然这里的要用的setter和getter函数太多了
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Show method for STID objects
#'
#' Displays a concise summary of an STID object, including Seurat components
#' and STID-specific analysis information.
#'
#' @param object An STID object
#'
#' @return Invisibly returns the object
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' ist_obj
#' }
#'
setMethod(
  f = "show",
  signature = "STID",
  # signature = "Seurat",
  definition = function(object) {
    callNextMethod() # show Seurat part, please see the Seurat::show method
    cat('\n')

    # STID
    cat("STID_analysis:\n")
    STID_analysis <- object@STID_analysis # show must be object, not STID_obj
    cat("STID_info:\n")
    cat(" - samp_id:", paste0(STID_analysis@STID_info$samp_info$samp_id,collapse = ", "), "\n")
    invisible(capture.output(
      coord_colnm <- GetCoordInfo(object)
    ))
    # cat(" - coord_colnm:", paste(names(coord_colnm), coord_colnm, sep = " = ", collapse = ", "), "\n")
    cat(" - coord_colnm: \n"); print(coord_colnm)
    cat("meta_data_record:", nrow(STID_analysis@meta_data_record$meta_data_info), "\n")
    cat(" - meta_key:", paste(STID_analysis@meta_data_record$meta_data_info$meta_key, collapse = ", "), "\n")
    STID_SingleSamp <- STID_analysis@SingleSampNiche
    cat("SingleSampNiche:", length(STID_SingleSamp), "\n")
    if (length(STID_SingleSamp) > 0) {
      for (ss in STID_SingleSamp) {
        cat(" - ", ss@samp_info$samp_id, ": ",
            nrow(ss@niche_info)," niches", "\n", sep = "")
      }
    }
    STID_MultiSamp <- STID_analysis@MultiSampNiche
    cat("MultiSampNiche:", length(STID_MultiSamp), "\n")
    if (length(STID_MultiSamp) > 0) {
      for (ms in STID_MultiSamp) {
        cat(" - ", ms@samp_info$multi_id, ": ",
            length(ms@samp_info$samp_id), " samples", "; ",
            length(ms@niche_info)," niches", "\n", sep = "")
      }
    }
    cat("\n")
    # browser()
  }
)

#' Print method for STID objects
#'
#' Provides a detailed printout of STID object contents, including STID_info,
#' meta_data_record, SingleSampNiche, and MultiSampNiche information.
#'
#' @param x An STID object
#' @param ... Additional arguments passed to print
#'
#' @return Invisibly returns the object
#'
#' @export
#' @method print STID
#'
#' @examples
#' \dontrun{
#' print(ist_obj)
#' }
#'
print.STID <- function(x, ...) {
  cat("\n")
  cat("STID_analysis:\n")
  STID_analysis <- x@STID_analysis
  cat("\n")
  cat(">>> STID_info:\n")
  info_list <- c(STID_analysis@STID_info$samp_info,
                 STID_analysis@STID_info$data_info,
                 STID_analysis@STID_info$project_info)
  len_pathogen_genes <- length(info_list$pathogen_genes)
  if(len_pathogen_genes > 3){
    info_list$pathogen_genes <- paste0(paste0(info_list$pathogen_genes[1:3], collapse = ", "),
                                       ", ... (total ", len_pathogen_genes, " genes)")
  }else{
    info_list$pathogen_genes <-  paste0(info_list$pathogen_genes, collapse = ", ")
  }
  info_list <- lapply(info_list, function(x) {
    if(is.null(x)) {
      return("NULL")  # 或 return(NA) 或 return("")
    } else {
      return(x)
    }
  })
  print(as_tibble(info_list))
  cat("\n")
  cat(">>> meta_data_record:\n")
  print(as_tibble(STID_analysis@meta_data_record$meta_data_info))
  cat("\n")
  STID_SingleSamp <- STID_analysis@SingleSampNiche
  cat(">>> SingleSampNiche:", length(STID_SingleSamp), "\n")
  if (length(STID_SingleSamp) > 0) {
    for (ss in STID_SingleSamp) {
      cat(" - ", ss@samp_info$samp_id, ": ",
          nrow(ss@niche_info), " niches", "(", paste(ss@niche_info$niche_key, collapse = ", "), ")",
          "\n", sep = "")
    }
  }
  STID_MultiSamp <- STID_analysis@MultiSampNiche
  cat("\n")
  cat(">>> MultiSampNiche:", length(STID_MultiSamp), "\n")
  if (length(STID_MultiSamp) > 0) {
    for (ms in STID_MultiSamp) {
      cat(" - ", ms@samp_info$multi_id, ": ",
          length(ms@samp_info$samp_id), " samples", "; ",
          length(ms@niche_info), " niches", "(", paste(names(ms@niche_info), collapse = ", "), ")", # !!! length(ms@niche_info)
          "\n", sep = "")
    }
  }

  cat("\n")
  invisible(x)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# check_loop_single/check_loop_multi
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check and process single-sample loop identifiers
#'
#' Internal function to validate and process loop_id for single-sample operations.
#'
#' @param STID_obj An STID object
#' @param valid_single Character vector, valid sample identifiers (optional)
#' @param loop_id Character, loop identifier(s) to check
#' @param mode Integer, mode for determining valid samples (1: all samples, 2: SingleSampNiche samples)
#'
#' @return Character vector of validated sample identifiers
#'
#' @keywords internal
#'
#' @noRd
#'
.check_loop_single <- function(STID_obj = NULL,valid_single = NULL,loop_id = NULL,mode = 2){
  #> loop_single
  if(is.null(valid_single)){
    if(mode == 1){
      valid_single <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "samp_id")[[1]]
    }else if(mode == 2){
      valid_single <- STID_obj@STID_analysis@SingleSampNiche %>% names()
    }
  }

  if(!all(loop_id %in% c("LoopAllSamp", valid_single))){
    clog_error(paste0("loop_id must be 'LoopAllSamp' or in: ", paste(valid_single, collapse = ", ")))
  }
  # clog_normal(paste0("Your loop_id: ", paste(loop_id, collapse = ", ")))
  if("LoopAllSamp" %in% loop_id){
    loop_single <- valid_single
  }else{
    loop_single <- loop_id
  }
  if(!all(loop_single %in% valid_single)){
    clog_error(paste0("Some samp_id not found in SingleSampNiche: ",
                      paste(setdiff(loop_single, valid_single), collapse = ", ")))
  }
  return(loop_single)
}

#' Check and process multi-sample loop identifiers
#'
#' Internal function to validate and process loop_id for multi-sample operations.
#'
#' @param STID_obj An STID object
#' @param loop_id Character, loop identifier(s) to check
#'
#' @return Character vector of validated multi-sample identifiers
#'
#' @keywords internal
#'
#' @noRd
#'
.check_loop_multi <- function(STID_obj,loop_id){
  #> loop_multi
  valid_multi <- STID_obj@STID_analysis@MultiSampNiche %>% names()
  if(!all(loop_id %in% c("LoopAllMulti", valid_multi))){
    clog_error(paste0("loop_id must be 'LoopAllMulti' or in: ", paste(valid_multi, collapse = ", ")))
  }
  # clog_normal(paste0("Your loop_id: ", paste(loop_id, collapse = ", ")))
  if("LoopAllMulti" %in% loop_id){
    loop_multi <- valid_multi
  }else{
    loop_multi <- loop_id
  }
  if(!all(loop_multi %in% valid_multi)){
    clog_error(paste0("Some loop_id not found in MultiSampNiche: ",
                      paste(setdiff(loop_multi, valid_multi), collapse = ", ")))
  }
  return(loop_multi)
}

.check_niche_exist <- function(STID_obj, samp_id, niche_key, samp_mode  = "SS"){
  if(samp_mode  == "SS"){
    valid_niche <- STID_obj@STID_analysis@SingleSampNiche[[samp_id]]@niche_cells %>% names()
  }else if(samp_mode  == "MS"){
    valid_niche <- STID_obj@STID_analysis@MultiSampNiche[[samp_id]]@niche_cells %>% names()
  }else{
    clog_error("samp_mode  must be 'SS' or 'MS'")
  }
  if(niche_key %in% valid_niche){
    return(TRUE)
  }else{
    return(FALSE)
  }
}




