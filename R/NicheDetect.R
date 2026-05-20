

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NicheDetect_Lasso
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Interactive lasso selection for ROI (Region of Interest) detection
#'
#' This function launches an interactive Shiny application that allows users to
#' manually draw lasso selections to define regions of interest (ROIs) on spatial
#' transcriptomics data. Supports multiple ROIs per sample and saves the results
#' for downstream analysis.
#'
#' @param STID_obj An STID object containing spatial transcriptomics data
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param meta_key Character, metadata key containing spatial coordinates and grouping
#'        information (default: "coord")
#' @param group_by Character, column name for coloring points in the interactive plot
#' @param lasso_para List, parameters for lasso selection when STID_obj is NULL, including:
#'       - meta_data: Data frame containing spatial coordinates and metadata (required if STID_obj is NULL)
#'       - samp_colnm: Column name for sample identifiers in meta_data (required if STID_obj is NULL)
#'       - x_colnm: Column name for x coordinates in meta_data (default: "x")
#'       - y_colnm: Column name for y coordinates in meta_data (default: "y")
#'       - data_format: Format of spatial data, either "square_grid" or "hex_grid" (default: "square_grid")
#'       - data_platform: Platform of spatial data, either "StereoSeq" or "Visium" (default: "StereoSeq")
#'       - binsize: Bin size for spatial data (default: 1)
#'       - coord_interval: Coordinate interval for spatial data (default: 1)
#' @param col Color palette for visualization (default: COLOR_LIST$PALETTE_WHITE_BG)
#' @param description Character, description of the analysis (default: NULL)
#' @param grp_nm Character, group name for output organization (default: NULL, uses timestamp)
#' @param dir_nm Character, directory name for output (default: "M2_NicheDetect_Lasso")
#'
#' @return Returns the modified STID object with added metadata containing ROI labels,
#'         distances to ROI centers, and edge information
#'
#' @import shiny
#' @import ggplot2
#' @importFrom sp point.in.polygon
#' @import dplyr
#' @import shinycssloaders
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch interactive ROI selection
#' STID_obj <- NicheDetect_Lasso(
#'   STID_obj = ist_object,
#'   group_by = "cell_type",
#'   meta_key = "coord"
#' )
#' }
NicheDetect_Lasso <- function(STID_obj = NULL,
                              loop_id = "LoopAllSamp", # must be
                              meta_key = "coord",
                              group_by = NULL,
                              lasso_para = list(
                                meta_data = NULL,
                                samp_colnm = NULL,
                                x_colnm = "x",
                                y_colnm = "y",
                                data_format = "square_grid",
                                data_platform = "StereoSeq",
                                binsize = 1,
                                coord_interval = 1
                              ),
                              col = COLOR_LIST$PALETTE_WHITE_BG,
                              description = NULL,
                              grp_nm = NULL, dir_nm = "M2_NicheDetect_Lasso"
                              ){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_check()
  if (!inherits(STID_obj, "STID")) {
    if(is.null(lasso_para$meta_data)){
      clog_error("The input STID_obj is not a valid STID object, and the lasso_para$meta_data is also NULL.")
    }else{
      clog_warn("The input STID_obj is not a valid STID object, but the lasso_para$meta_data is provided.
      When you use the lasso_para$meta_data to run the lasso selection, you must provide group_by and lasso_para parameters.")
    }
  }
  clog_normal(paste0("Your ggplot2 version is: ", packageVersion("ggplot2")))
  .check_null_args(meta_key, group_by)

  #>
  if(inherits(STID_obj, "STID")){
    loop_single <- .check_loop_single(STID_obj = STID_obj,loop_id = loop_id, mode = 1)
    all_single <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "samp_id")[[1]]
  }else{
    meta_data <- lasso_para$meta_data
    loop_single <- lasso_para$meta_data[[lasso_para$samp_colnm]] %>% unique() %>% as.character()
    all_single <- lasso_para$meta_data[[lasso_para$samp_colnm]] %>% unique() %>% as.character()
  }
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  if(inherits(STID_obj, "STID")){
    meta_data <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]]
    samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
    x_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "x_colnm")[[1]]
    y_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "y_colnm")[[1]]
    data_format <- GetInfo(STID_obj, info_key = "data_info",sub_key = "data_format")[[1]]
    data_platform <- GetInfo(STID_obj, info_key = "data_info",sub_key = "data_platform")[[1]]
    binsize <- GetInfo(STID_obj, info_key = "data_info",sub_key = "binsize")[[1]]
    coord_interval <- GetInfo(STID_obj, info_key = "data_info",sub_key = "coord_interval")[[1]]
  }else{
    samp_colnm <- lasso_para$samp_colnm
    x_colnm <- lasso_para$x_colnm
    y_colnm <- lasso_para$y_colnm
    data_format <- lasso_para$data_format
    data_platform <- lasso_para$data_platform
    binsize <- lasso_para$binsize
    coord_interval <- lasso_para$coord_interval
  }
  if(data_format  == "square_grid"){
    k <- 8; edge_thres <- 1
  }else if(data_format == "hex_grid" & data_platform == "Visium"){
    k <- 6; edge_thres <- 1
  }
  meta2STID_list <- list()
  n <- 1
  for(i in seq_along(all_single)){
    i_single <- all_single[i]
    i_coord_interval <- coord_interval[i]
    i_meta_data <- meta_data %>% filter(!!sym(samp_colnm) == i_single)
    if(!i_single %in% loop_single){
      meta2STID_list[[i_single]] <- i_meta_data
      next
    }
    clog_loop(paste0("Processing samp_id: ", i_single, " (", n, "/", length(loop_single), ")"))
    n <- n+1
    i_now_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
    file_nm <- paste0("NicheDetect_",i_single,"_",i_now_time, ".zip")
    i_output_dir <- paste0(output_dir,"/",i_single,"/")
    i_photo_dir <- paste0(photo_dir,"/",i_single,"/")
    dir.create(i_output_dir,recursive = TRUE,showWarnings = FALSE)
    dir.create(i_photo_dir,recursive = TRUE,showWarnings = FALSE)
    clog_warn("Please make sure to finish the lasso selection and save the results before closing the app.")
    print(
      .NicheDetect_lasso_shiny(meta_data = i_meta_data, samp_id = paste0(i_single," (", i, "/", length(loop_single), ")"),
                               group_by = group_by,x_colnm = x_colnm, y_colnm = y_colnm,
                               col = col,
                               file_nm =  file_nm ,output_dir = i_output_dir, photo_dir = i_photo_dir)
    )

    # >
    ROI_meta_path <- paste0(i_output_dir,"/","NicheDetect_Lasso_points.txt")
    if(!file.exists(ROI_meta_path)){
      clog_warn(paste0("The ROI meta file not found for samp_id: ", i_single))
      clog_warn("Skip this samp and  continue to the next samp.")
      meta2STID_list[[i_single]] <- i_meta_data
      next
    }
    clog_normal(paste0("Read ROI meta data from: ", ROI_meta_path))
    ROI_meta <- read.table(ROI_meta_path,
                           row.names = 1,header = TRUE,sep = "\t",stringsAsFactors = FALSE)

    clog_normal(paste0("The Results of Point in hull detection:"))
    print(
      table(ROI_meta["ROI_label"])
    )
    if(max(ROI_meta$ROI_count,na.rm = TRUE) >1){
      clog_warn("Detected points belong to multiple ROIs, these points will be removed in final results.")
    }
    res_meta_data <- ROI_meta %>%
      filter(ROI_count <2)

    #> CalDist2Center
    common_ids <- intersect(rownames(i_meta_data), rownames(res_meta_data))
    i_meta_data$ROI_label   <- "neg" # NA -> "neg"
    i_meta_data$ROI_count   <- NA_integer_
    i_meta_data[common_ids, "ROI_label"]   <- res_meta_data[common_ids, "ROI_label"]
    i_meta_data[common_ids, "ROI_count"]   <- res_meta_data[common_ids, "ROI_count"]
    i_meta_data <- .CalROIedge(meta_data = i_meta_data,
                               ROI_label_colnm = "ROI_label",
                               ROI_edge_colnm = "ROI_edge",
                               k = k,edge_thres = edge_thres)
    i_meta_data <- .CalROIcenter(meta_data = i_meta_data,
                                 ROI_label_colnm = "ROI_label",
                                 ROI_center_colnm = "ROI_center",
                                 coord_interval = i_coord_interval)
    i_meta_data <- .CalDist2Center(meta_data = i_meta_data,
                                   ROI_label_colnm = "ROI_label",
                                   ROI_center_colnm = "ROI_center",
                                   CalMode = "ROI")
    meta2STID_list[[i_single]] <- i_meta_data
  }

  #> AddMetaData
  meta2STID <- bind_rows(meta2STID_list) %>%
    as.data.frame()
  meta2STID <- meta2STID[rownames(meta_data), , drop = FALSE] # keep order
  new_meta_key <- paste0(dir_nm,"_",grp_nm)
  if (inherits(STID_obj, "STID")) {
    STID_obj <- AddMetaData(STID_obj = STID_obj,
                            meta_key = new_meta_key,
                            add_data = meta2STID,
                            dir_nm = dir_nm,
                            grp_nm = grp_nm,
                            asso_key = meta_key,
                            description = description)
  }else{
    STID_obj <- meta2STID
  }


  # >>> Final
  .save_function_params("NicheDetect_Lasso", envir = environment(), file = paste0(output_dir,"Log_function_params_(NicheDetect_Lasso).log") )
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(NicheDetect_Lasso).log")) %>% invisible()
  return(STID_obj)
}

#' Internal Shiny application for lasso-based ROI selection
#'
#' This function creates and runs a Shiny app for interactive ROI selection.
#' It is called internally by NicheDetect_Lasso and should not be used directly.
#'
#' @param meta_data Data frame containing spatial coordinates and metadata
#' @param samp_id Character, sample identifier for display
#' @param x_colnm Character, column name for x coordinates
#' @param y_colnm Character, column name for y coordinates
#' @param group_by Character, column name for coloring points
#' @param col Color palette for visualization
#' @param file_nm Character, filename for saving results
#' @param output_dir Character, output directory path
#' @param photo_dir Character, directory path for saving plots
#'
#' @return Runs a Shiny application (no return value)
#'
#' @import shiny
#' @import ggplot2
#' @importFrom sp point.in.polygon
#' @import dplyr
#' @importFrom shinycssloaders withSpinner
#' @import zip
#'
#' @keywords internal
#'
#' @noRd
.NicheDetect_lasso_shiny <- function(meta_data = NULL,samp_id = NULL,
                                     x_colnm = NULL, y_colnm = NULL,group_by = NULL,
                                     col = NULL,
                                     file_nm = NULL,output_dir = NULL, photo_dir = NULL){

  clog_normal("Launch Multi-ROI Lasso Shiny App...")
  ui <- fluidPage(
    titlePanel(strong("NicheDetect Lasso Tool"), windowTitle = "NicheDetect Lasso Tool"),

    sidebarLayout(
      sidebarPanel(
        # 数据集选择wellPanel
        wellPanel(
          h4(strong("Input data"), style = "margin-top: 0px;"),

          fluidRow(
            column(6,
                   selectInput("dataset",
                               label = tags$label("Select Dataset:", style = "font-weight: normal;"),
                               choices = c("Input data","upload files", "iris(Demo)", "mtcars(Demo)"),
                               selected = "Input data"), # 首次input$dataset不会触发reactive和observe,但是有值,必须设置默认select
            ),
            column(6,
                   h5("Sample name", style = "margin-top: 3px;"),
                   textOutput("samp_id")
            )
          ),

          #> conditionalPanel
          conditionalPanel(
            condition = "input.dataset != 'upload files'",
            fluidRow(
              column(4,
                     selectInput("X_axis",
                                 tags$label("X-axis:", style = "font-weight: normal;"), choices = NULL)
              ),
              column(4,
                     selectInput("Y_axis",
                                 tags$label("Y-axis:", style = "font-weight: normal;"), choices = NULL)
              ),
              column(4,
                     selectInput("Value",
                                 tags$label("Value:", style = "font-weight: normal;"),choices = NULL)
              ),
            )
          ),
          conditionalPanel(
            condition = "input.dataset == 'upload files'",
            fluidRow(
              column(12,
                     fileInput("custom_file",
                               label = tags$label("Upload your data (TXT or CSV):",
                                                  style = "font-weight: normal;"),
                               buttonLabel = "Browse",
                               accept = c(".txt", ".csv"))
              )
            ),
            fluidRow(
              column(4,
                     textInput("custom_x",
                               tags$label("X-axis:", style = "font-weight: normal;"),
                               value = "")
              ),
              column(4,
                     textInput("custom_y",
                               tags$label("Y-axis:", style = "font-weight: normal;"),
                               value = "")
              ),
              column(4,
                     textInput("custom_value",
                               tags$label("Value:", style = "font-weight: normal;"),
                               value = "")
              )
            )
          ),
          fluidRow(
            column(6, align = "right",
                   actionButton("start_plot", "Plot", class = "btn-info" , style = "width: 90%; margin-top: 10px;")
            ),
            column(6, align = "left",
                   actionButton("Next Sample", "Next Sample", class = "btn-danger", style = "width: 90%; margin-top: 10px;")
            )
          )
        ),

        #> Graphic parameter
        wellPanel(
          h4(strong("Graphic parameter"), style = "margin-top: 0px;"),
          fluidRow(
            column(4,
                   numericInput("point_size",
                                tags$label("Point size:", style = "font-weight: normal;"),
                                value = 1.5,min = 0,max = NA,step = 0.5),
                   textInput("ROI_edge_col",
                             tags$label("ROI edge color:", style = "font-weight: normal;"),
                             value = "black")
            ),
            column(4,
                   numericInput("point_alpha",
                                tags$label("Point alpha:", style = "font-weight: normal;"),
                                value = 0.6,min = 0,max = 1,step = 0.1),
                   numericInput("ROI_alpha",
                                tags$label("ROI alpha:", style = "font-weight: normal;"),
                                value = 0.5,min = 0,max = 1,step = 0.1)
            ),
            column(4,align = "left", # right
                   style = "display: block; margin-top: 30px;",
                   checkboxInput("coord_fixed", tags$label("coord fixed", style = "font-weight: normal;"),
                                 value = TRUE, width = NULL),
                   checkboxInput("y_reverse", tags$label("y reverse", style = "font-weight: normal;"),
                                 value = TRUE, width = NULL),
                   checkboxInput("plot_grid", tags$label("plot grid", style = "font-weight: normal;"),
                                 value = FALSE, width = NULL)
            ),
          )
        ),

        #> Lasso 操作
        wellPanel(
          h4(strong("Lasso selection"), style = "margin-top: 0px;"),
          fluidRow(
            style = "width: 100%; margin:5px 0px 5px 0px;padding:0px 0px 0px 0px",
            column(6, actionButton("new_lasso", "New", class = "btn-primary", style = "width: 100%; margin-top: 10px;")),
            column(6, actionButton("finish_lasso", "Finish", class = "btn-success", style = "width: 100%; margin-top: 10px;")),
            column(6, actionButton("clear_all", "Clear All", class = "btn-warning", style = "width: 100%; margin-top: 10px;")),
            column(6, downloadButton("download_data", "Save", class = "btn-info ", style = "width: 100%; margin-top: 10px;"))
          ),
          br(),
          # 视图控制
          wellPanel(
            h5(strong("Zoom & Pan"), style = "margin-top: 0px;"),
            fluidRow(
              column(1),
              column(1, actionButton("zoom_in", "", icon = icon("search-plus"),
                                     class = "btn-large-icon", title = "Zoom In")),
              column(1, actionButton("zoom_out", "", icon = icon("search-minus"),
                                     class = "btn-large-icon", title = "Zoom Out")),
              column(1),
              column(1, actionButton("pan_left", "", icon = icon("arrow-left"),
                                     class = "btn-large-icon")),
              column(1, actionButton("pan_right", "", icon = icon("arrow-right"),
                                     class = "btn-large-icon")),
              column(1, actionButton("pan_up", "", icon = icon("arrow-up"),
                                     class = "btn-large-icon")),
              column(1, actionButton("pan_down", "", icon = icon("arrow-down"),
                                     class = "btn-large-icon")),
              column(1),
              column(1, actionButton("reset_view", "", icon = icon("sync"),
                                     class = "btn-large-icon", title = "Reset"))

            )
          ),
          fluidRow(
            column(
              width = 5,
              h5("Label:", style = "margin-top: 3px;"),
              textInput("label", NULL, value = "ROI 1"),
              # br(),
              h5("Point coord:", style = "margin-top: 3px;"),
              verbatimTextOutput("Point_coord", placeholder = TRUE)
            ),

            column(
              width = 7,
              h5("ROIs info:", style = "margin-top: 2px;"),
              verbatimTextOutput("ROIs_info", placeholder = TRUE)
            )
          )
        ),
      ),

      # 主面板：绘图
      mainPanel(
        tabsetPanel(
          tabPanel(
            title = "Input Data Preview",
            wellPanel(
              h3(strong("Input Data Preview"),
                 style = "margin: 10px 0 10px 0px; font-weight: 500; color: #333; text-align: center;"),
              DT::DTOutput("head_table"),
              br(),
              fluidRow(
                column(12,
                       h4("Vars info:", style = "margin-top: 2px; font-weight: bold;"),
                       verbatimTextOutput("Vars_info", placeholder = TRUE))
              ),

              style = "
              padding: 20px;
              background-color: #ffffff;
              border-radius: 12px;
              box-shadow: 0 4px 12px rgba(0,0,0,0.1);
              border: 1px solid #e0e0e0;
              margin: 10px;
              height: 1000px;
              width: 97%;
              overflow: auto;

              /* Flex */
              display: flex;
              flex-direction: column;
              justify-content: start;
              align-items: stretch;
              "
            )
          ),
          tabPanel(
            title = "Data Visualization",
            wellPanel(
              h3(strong("Input Data Visualization"),
                 style = "margin: 10px 0 10px 0px; font-weight: 500; color: #333; text-align: center;"),
              h4("Click to add a point; Double-click or Click the finish to finish the polygon; Click new to remove the current polygon.",
                 style = "margin-left: 20px;"),

              plotOutput("plot",
                         width = "100%",
                         height = "875px",
                         # height = "100%",
                         click = "plot_click",
                         dblclick = "plot_dblclick",
                         brush = brushOpts(id = "plot_brush", resetOnNew = TRUE)) %>%
                withSpinner(hide.ui = FALSE), # hide.ui = FALSE 会让 spinner 叠加在旧图表上

              style = "
              padding: 10px;
              background-color: #ffffff;
              border-radius: 12px;
              box-shadow: 0 4px 12px rgba(0,0,0,0.1);
              border: 1px solid #e0e0e0;
              margin: 10px;
              height: 1000px;
              width: 97%;
              overflow: auto;

              /* Flex */
              display: flex;
              flex-direction: column;
              justify-content: start;
              align-items: stretch;
              "
            )
          )
        )  # end of tabsetPanel
      )  # end of mainPanel
    ),

    # 自定义内联CSS
    tags$head(
      tags$style(HTML("
      .well {
        border-radius: 8px;
        border: 1px solid #ddd;
        padding: 12px;
        margin-bottom: 10px;
        background-color: #F4F9FC;
      }
      .btn {
        margin: 2px 0;
      }
      .btn-large-icon {
        width: 40px;
        height: 40px;
        padding: 8px;
        font-size: 18px;
        border-radius: 8px;
      }
      .btn-small-icon {
        width: 30px;
        height: 30px;
        padding: 5px;
        font-size: 12px;
      }
      h4 {
        margin-bottom: 10px;
        color: #343a40;
        font-weight: 500;
      }
      .sidebar-panel {
        font-size: 14px;
      }
      #ROIs_info {
        height: 130px;
        overflow-y: auto;
        font-family: 'Courier New', monospace;
        font-size: 0.9em;
        background-color: #ffffff;
        padding: 8px;
        border-radius: 4px;
      }
      #Vars_info {
        height: 200px;
        overflow-y: auto;
        font-family: 'Courier New', monospace;
        font-size: 0.9em;
        background-color: #ffffff;
        padding: 8px;
        border-radius: 4px;
      }
      #Point_coord {
        height: 50px;
        overflow-y: auto;
        font-family: 'Courier New', monospace;
        font-size: 0.9em;
        background-color: #ffffff;
        padding: 8px;
        border-radius: 4px;
      }
      #samp_id {
        height: 35px;
        weight: bold;
        overflow-y: auto;
        font-family: 'Courier New', monospace;
        font-size: 0.9em;
        background-color: #ffffff;
        padding: 8px;
        border-radius: 4px;
      }
      .waiter {
        z-index: 9999 !important;
      }
      .waiter > div {
        z-index: 9999 !important;
      }
    "))
    )
  )

  server <- function(input, output, session) {
    # 存储所有点数据
    values <- reactiveValues(
      all_points = data.frame(x = numeric(), y = numeric(), group = integer()),
      current_points = data.frame(x = numeric(), y = numeric()),
      polygons = list(),
      current_polygon = NULL,
      polygon_counter = 1,
      x_range = NULL,
      y_range = NULL,
      plot_data = NULL
    )
    currentData <- NULL
    currentVars <- NULL

    #>
    observe({
      req(input$dataset) # 值非空就执行，req(input$dataset)发生了变化，所以observe自然就执行了
      currentData <<- switch(input$dataset,
                             "Input data" = {
                               req(meta_data, message = "Input meta_data is NULL.")
                               updateSelectInput(session, "X_axis", choices = names(meta_data), selected = x_colnm)
                               updateSelectInput(session, "Y_axis", choices = names(meta_data), selected = y_colnm)
                               updateSelectInput(session, "Value", choices = names(meta_data), selected = group_by)
                               meta_data
                             },
                             "iris(Demo)" = {
                               updateSelectInput(session, "X_axis", choices = names(iris), selected = "Sepal.Length")
                               updateSelectInput(session, "Y_axis", choices = names(iris), selected = "Sepal.Width")
                               updateSelectInput(session, "Value", choices = names(iris), selected = "Species")
                               iris
                             },
                             "mtcars(Demo)" = {
                               updateSelectInput(session, "X_axis", choices = names(mtcars), selected = "wt")
                               updateSelectInput(session, "Y_axis", choices = names(mtcars), selected = "mpg")
                               updateSelectInput(session, "Value", choices = names(mtcars), selected = "cyl")
                               mtcars
                             },
                             "upload files" = {
                               if(is.null(input$custom_file)) {
                                 showNotification("Please upload a file.", type = "error", duration = 5)
                                 return()
                                 req(input$custom_file, message = "Please upload a file.")
                               }
                               ext <- tools::file_ext(input$custom_file$name)
                               if (ext == "csv") {
                                 read.csv(input$custom_file$datapath, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
                               } else if (ext == "txt") {
                                 read.table(input$custom_file$datapath, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
                               } else {
                                 stop("Unsupported file type. Please upload a CSV or TXT file.")
                               }
                             },
                             stop("Unknown dataset")
      )

      #
      output$head_table <- DT::renderDT({
        if(nrow(currentData)>200){
          showNotification("Data has more than 200 rows, showing only the first 200 rows.", type = "warning")
          DT::datatable(head(currentData,n = 200), caption = paste0("Warnings! This data has ",nrow(currentData), " rows, but only shows 200 rows."))
        }else{
          DT::datatable(currentData, caption = "")
        }
      })
    })

    #> Next Sample
    observeEvent(input$`Next Sample`, {
      showNotification("App is closing...", type = "message")
      stopApp()
    })

    #> plot
    observeEvent(input$start_plot, {
      if (input$dataset == "upload files" & is.null(input$custom_file)) {
        showNotification("❌ ERROR: Please upload a file.", type = "error", duration = 5)
        return()
        req(input$custom_file, message = "Please upload a file.") # req() 的 message 只有在被 render* 函数（如 renderPlot, renderTable）直接调用时才会显示在 UI 上。
      }

      currentVars <<- if(input$dataset != "upload files"){
        list(x = input$X_axis, y = input$Y_axis, value = input$Value, title_nm = input$dataset)
      }else{
        if(is.null(input$custom_file)) {
          showNotification("Please upload a file.", type = "error", duration = 5)
          return()
          req(input$custom_file, message = "Please upload a file.")
        }
        title_nm <- tools::file_path_sans_ext(basename(input$custom_file$name))
        list(x = input$custom_x, y = input$custom_y, value = input$custom_value, title_nm = title_nm)
      }

      #
      # validate(
      #   need(nrow(currentData) > 0, showNotification("The data has 0 rows. Please check your file.")),
      #   need(length(currentVars$x) > 0 && !is.null(currentVars$x), 0, showNotification("X variable is not specified.")),
      #   need(length(currentVars$y) > 0 && !is.null(currentVars$y), 0, showNotification("Y variable is not specified.")),
      #   need(length(currentVars$value) > 0 && !is.null(currentVars$value), 0, showNotification("Value variable is not specified.")),
      #   need(currentVars$x %in% names(currentVars), 0, showNotification(paste("Column not found: x =", currentVars$x))),
      #   need(currentVars$y %in% names(currentVars), 0, showNotification(paste("Column not found: y =", currentVars$y))),
      #   need(currentVars$value %in% names(currentVars), 0, showNotification(paste("Column not found: value =", currentVars$value)))
      # )
      if (nrow(currentData) == 0) {
        showNotification("❌ ERROR: The dnput ata has 0 rows. Please check your file.", type = "error", duration = 5)
        return()
      }
      if (!(currentVars$x %in% names(currentData))) {
        showNotification(paste("❌ ERROR: Column not found: x =", currentVars$x), type = "error", duration = 5)
        return()
      }
      if (!(currentVars$y %in% names(currentData))) {
        showNotification(paste("❌ ERROR: Column not found: y =", currentVars$y), type = "error", duration = 5)
        return()
      }
      if (!(currentVars$value %in% names(currentData))) {
        showNotification(paste("❌ ERROR: Column not found: value =", currentVars$value), type = "error", duration = 5)
        return()
      }

      #
      values$plot_data <- data.frame(
        x = currentData[[currentVars$x]],
        y = currentData[[currentVars$y]],
        value = currentData[[currentVars$value]]
      )
      values$x_range <- range(values$plot_data$x, na.rm = TRUE)
      values$y_range <- range(values$plot_data$y, na.rm = TRUE)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)  # 不在初始化时触发

    # 主绘图函数：目前适用于4.0.2以上？其它版本没测试过
    make_plot <- function() {
      req(values$plot_data)


      p <- ggplot() +
        geom_point(data = values$plot_data, aes(x = x, y = y, color = value),
                   size = input$point_size,alpha = input$point_alpha,
                   shape = 16,stroke = 0) +
        theme_test() +
        labs(title = currentVars$title_nm, x = currentVars$x, y = currentVars$y, color = currentVars$value) +
        scale_color_manual(values = col) +
        theme(
          # legend.position = "none",
          # panel.grid = element_blank(),
          plot.margin = margin(0.4,0.4,0.4,0.4,'cm'),
          plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
          axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
          # axis.text.x = element_text(size = 10,face = 'bold',hjust = 1,angle = 45),
          axis.text.x = element_text(size = 10,face = 'bold',hjust = 0.5),
          axis.text.y  = element_text(size = 10,face = 'bold',hjust = 0.5),
          legend.title = element_text(size = 11,face = 'bold',hjust = 0),
          legend.text = element_text(size = 10,face = 'bold',hjust = 0),
          legend.key.size = unit(20, "pt")
        )
      # browser()

      #
      if(input$plot_grid){
        p <- p + theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
                       panel.grid.minor = element_line(color = "grey90", linewidth = 0.25))
      }
      if (!is.null(values$x_range) && !is.null(values$y_range)) {
        x_left_pos <- values$x_range[1] - 0.1 * abs(values$x_range[1])
        x_right_pos <- values$x_range[2] + 0.1 * abs(values$x_range[2])
        y_bottom_pos <- values$y_range[1] - 0.1 * abs(values$y_range[1])
        y_top_pos <- values$y_range[2] + 0.1 * abs(values$y_range[2])

        if(packageVersion("ggplot2") > 3.5){
          if(input$y_reverse){ # 旧版本<3.4这个要放在coord_fixed前面有效？
            # 这个3.4版本无效，3.5版本可用，但是又出现点击点坐标不对问题。因为
            # ≥ 3.4.4新版提供了 dimension（像素空间）和 continuous_range（数据空间）？但是4.0.2又解决了
            p <- p + scale_y_reverse()
          }
          if(input$coord_fixed){
            p <- p + coord_fixed(ratio = 1,
                                 xlim = c(x_left_pos, x_right_pos),
                                 ylim = c(y_bottom_pos, y_top_pos))
          }
        }else{
          if(input$coord_fixed){
            if(input$y_reverse){ # 3.4.4版本目前可以这么设置
              p <- p + coord_fixed(ratio = 1,
                                   xlim = c(x_left_pos, x_right_pos),
                                   ylim = c( y_top_pos,y_bottom_pos))
            }else{
              p <- p + coord_fixed(ratio = 1,
                                   xlim = c(x_left_pos, x_right_pos),
                                   ylim = c(y_bottom_pos, y_top_pos))
            }
          }else{
            if(input$y_reverse){
              p <- p +
                coord_cartesian(
                  xlim = c(x_left_pos, x_right_pos),
                  ylim = c(y_top_pos,y_bottom_pos))
            }else{
              p <- p +
                coord_cartesian(
                  xlim = c(x_left_pos, x_right_pos),
                  ylim = c(y_bottom_pos, y_top_pos))
            }
          }
        }

      }else{
        showNotification("The x_range or y_range is NULL, please click the Plot button to set the ranges.", type = "error")
      }

      # 绘制所有已完成的套索区域
      if (length(values$polygons) > 0) {
        for (i in seq_along(values$polygons)) { # list每个区域分开画
          poly_df <- values$polygons[[i]]$points
          label <- values$polygons[[i]]$label
          if(length(values$polygons)<10){
            p <- p +
              geom_polygon(data = poly_df, aes(x = x, y = y),
                           fill = rainbow(10)[i], alpha = input$ROI_alpha,
                           color = input$ROI_edge_col, linewidth = 0.5)
          }else{
            p <- p +
              geom_polygon(data = poly_df, aes(x = x, y = y),
                           fill = rainbow(length(values$polygons))[i], alpha = input$ROI_alpha)
          }
          p <- p +
            geom_text(data = data.frame(x = mean(poly_df$x), y = mean(poly_df$y)),
                      aes(x = x, y = y), label = label, size = 6)
        }
      }

      # 绘制当前正在绘制的套索区域
      num_current_points <- nrow(values$current_points)
      if (num_current_points > 1) {
        p <- p +
          geom_path(data = values$current_points, aes(x = x, y = y),
                    color = input$ROI_edge_col, linewidth = 1) +
          geom_point(data = values$current_points, aes(x = x, y = y),
                     color = input$ROI_edge_col, size = 3)

        if (num_current_points > 2) {
          closed_path <- rbind(values$current_points, values$current_points[1,])
          p <- p + geom_path(data = closed_path, aes(x = x, y = y),
                             color = input$ROI_edge_col, linetype = "dashed", linewidth = 1)
        }
      }else if (num_current_points > 0){
        p <- p +
          geom_point(data = values$current_points, aes(x = x, y = y),
                     color = input$ROI_edge_col, size = 3)
      }
      return(p)  # 返回图形对象
    }
    output$plot <- renderPlot({
      make_plot()
    })

    #> 处理点击事件 - observeEvent
    # 添加点到当前套索
    # ggplot2 <3.4和>=3.4的方法不同
    observeEvent(input$plot_click, {
      if (!is.null(input$plot_click)) {
        if(1){ # bug未解决
          # browser()
          new_point <- data.frame(x = input$plot_click$x, y = input$plot_click$y) # 旧版
          new_point <- data.frame(x = input$plot_click$x, y = input$plot_click$y)
          values$current_points <- rbind(values$current_points, new_point)
        }else{
          # browser()
          ggb <- ggplot2::ggplot_build(make_plot())
          panel <- ggb$layout$panel_params
          x_pos <- panel[[1]]$x.range[1]+(panel[[1]]$x.range[2]-panel[[1]]$x.range[1])*input$plot_click$x
          y_pos <- panel[[1]]$y.range[1]+(panel[[1]]$y.range[2]-panel[[1]]$y.range[1])*input$plot_click$y
          if(input$y_reverse){
            y_pos <- -y_pos
          }
          new_point <- data.frame(x = x_pos, y = y_pos)
          values$current_points <- rbind(values$current_points, new_point)
        }
      }
    })

    # 处理双击事件 - 完成当前套索
    observeEvent(input$plot_dblclick, {
      if (nrow(values$current_points) > 2) {
        closed_points <- rbind(values$current_points, values$current_points[1,]) # 闭合多边形
        values$polygons[[length(values$polygons) + 1]] <- list( # 保存多边形
          points = closed_points,
          label = input$label
        )
        values$current_points <- data.frame(x = numeric(), y = numeric()) # 重置当前点
        updateTextInput(session, "label", value = paste0("ROI ", length(values$polygons) + 1)) # 更新标签输入框
      }else{
        showNotification("The lasso must have at least 3 points to form a polygon.", type = "error")
      }
    })

    # 开始新套索区域
    observeEvent(input$new_lasso, {
      values$current_points <- data.frame(x = numeric(), y = numeric())
    })

    # 完成当前套索区域
    observeEvent(input$finish_lasso, {
      if (nrow(values$current_points) > 2) {
        # 闭合多边形
        closed_points <- rbind(values$current_points, values$current_points[1,])
        # 保存多边形
        values$polygons[[length(values$polygons) + 1]] <- list(
          points = closed_points,
          label = input$label
        )

        # 重置当前点
        values$current_points <- data.frame(x = numeric(), y = numeric())

        # 更新标签输入框
        updateTextInput(session, "label", value = paste0("ROI ", length(values$polygons) + 1))
      }else{
        showNotification("The lasso must have at least 3 points to form a polygon.", type = "error")
      }
    })

    # 清除所有区域
    observeEvent(input$clear_all, {
      values$polygons <- list()
      values$current_points <- data.frame(x = numeric(), y = numeric())
      values$polygon_counter <- 1
      updateTextInput(session, "label", value = "ROI 1")
    })

    # 缩放和导航功能
    observeEvent(input$zoom_in, {
      if (!is.null(values$x_range)) {
        x_mid <- mean(values$x_range)
        y_mid <- mean(values$y_range)
        x_range <- diff(values$x_range) * 0.8 # 缩小0.2，需要增加0.25
        y_range <- diff(values$y_range) * 0.8

        values$x_range <- c(x_mid - x_range/2, x_mid + x_range/2)
        values$y_range <- c(y_mid - y_range/2, y_mid + y_range/2)
      }
    })

    observeEvent(input$zoom_out, {
      if (!is.null(values$x_range)) {
        x_mid <- mean(values$x_range)
        y_mid <- mean(values$y_range)
        x_range <- diff(values$x_range) * 1.25
        y_range <- diff(values$y_range) * 1.25

        values$x_range <- c(x_mid - x_range/2, x_mid + x_range/2)
        values$y_range <- c(y_mid - y_range/2, y_mid + y_range/2)
      }
    })

    observeEvent(input$pan_left, {
      if (!is.null(values$x_range)) {
        x_shift <- diff(values$x_range) * 0.2
        values$x_range <- values$x_range - x_shift
      }
    })

    observeEvent(input$pan_right, {
      if (!is.null(values$x_range)) {
        x_shift <- diff(values$x_range) * 0.2
        values$x_range <- values$x_range + x_shift
      }
    })

    observeEvent(input$pan_up, {
      if (!is.null(values$y_range)) {
        y_shift <- diff(values$y_range) * 0.2
        values$y_range <- values$y_range - y_shift
      }
    })

    observeEvent(input$pan_down, {
      if (!is.null(values$y_range)) {
        y_shift <- diff(values$y_range) * 0.2
        values$y_range <- values$y_range + y_shift
      }
    })

    observeEvent(input$reset_view, {
      if (!is.null(values$plot_data)) {
        values$x_range <- range(values$plot_data$x, na.rm = TRUE)
        values$y_range <- range(values$plot_data$y, na.rm = TRUE)
      }
    })

    # 显示区域信息
    output$ROIs_info <- renderPrint({
      if (length(values$polygons) == 0) {
        cat("No ROIs drawn")
      } else {
        for (i in seq_along(values$polygons)) {
          cat(sprintf("Num %d (%s): %d points\n",
                      i,
                      values$polygons[[i]]$label,
                      nrow(values$polygons[[i]]$points) - 1))
        }
      }
    })
    output$Vars_info <- renderPrint({
      if (is.null(values$plot_data) || nrow(values$plot_data) == 0) {
        cat("No input data")
      } else {
        paste0(currentVars$x,": min(",min(values$plot_data$x),"), max(",max(values$plot_data$x),")") %>%
          print()
        paste0(currentVars$y,": min(",min(values$plot_data$y),"), max(",max(values$plot_data$y),")") %>%
          print()
        if(is.numeric(values$plot_data$value)){
          paste0(currentVars$value,": min(",min(values$plot_data$value),"), max(",max(values$plot_data$value),")") %>%
            print()
        }else{
          print(paste0(currentVars$value,":"))
          print(table(values$plot_data$value))
        }
      }
    })
    output$Point_coord <- renderPrint({
      if (nrow(values$current_points) == 0) {
        cat("No points selected")
      } else {
        latest_point <- tail(values$current_points, 1)
        cat(sprintf("x: %.4f, y: %.4f", latest_point$x, latest_point$y))
      }
    })
    output$samp_id <- renderText({
      samp_id
    })

    # 下载区域数据
    output$download_data <- downloadHandler(
      filename = function() {
        file_nm
      },
      content = function(file) { # 下载这个file，文件名在filename了
        if (is.null(values$polygons) || length(values$polygons) == 0) {
          message_path <- file.path(output_dir, "NicheDetect_Lasso_error.txt")
          writeLines("No lasso ROIs available for download. Please create some ROIs first.", message_path)

          zip::zip(
            zipfile = file,
            files = "NicheDetect_Lasso_error.txt",
            root = output_dir
          )
          return()
        }

        # polygon_df
        polygon_df <- data.frame()
        for (i in seq_along(values$polygons)) {
          poly_data <- values$polygons[[i]]$points
          poly_data$ROI_label <- values$polygons[[i]]$label
          polygon_df <- rbind(polygon_df, poly_data)
        }
        txt_path <- file.path(output_dir, "NicheDetect_Lasso_ROIcoords.txt")
        write.table(polygon_df, txt_path, sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

        # slect_point_df
        spot_name <- rownames(currentData) # meta_data = currentData
        x_all <- currentData[[currentVars$x]] # currentVars$x
        y_all <- currentData[[currentVars$y]]
        polygon_list <- split(polygon_df, polygon_df["ROI_label"])
        list_inside <- lapply(polygon_list, function(poly) {
          px <- poly$x # 提取多边形 x, y 坐标
          py <- poly$y
          inside <- point.in.polygon(x_all, y_all, px, py) # 0=outside, 1,2=inside/border,3=vertex
          return(inside)
        })
        logic_inside <- lapply(list_inside, function(x) {
          return(x > 0)
        })
        logic_edge <- lapply(list_inside, function(x) {
          return(x > 1) # TRUE if inside，根据逻辑值取点
        })
        logic_inside_matrix <- bind_cols(logic_inside)
        logic_edge_matrix <- bind_cols(logic_edge)
        colnm <- colnames(logic_inside_matrix) # 多边形名字
        ROI_label <- apply(logic_inside_matrix, 1, function(row) {
          included_rois <- colnm[which(row)]
          if (length(included_rois) == 0) {
            return(NA_character_)
          } else {
            return(paste(included_rois, collapse = "; "))
          }
        })
        ROI_count <- rowSums(logic_inside_matrix)
        ROI_edge <- rowSums(logic_edge_matrix) > 0

        # >
        meta_data <- currentData %>% # !!!
          # cbind(data.frame(ROI_label,ROI_count,ROI_edge)) %>%
          cbind(data.frame(ROI_label,ROI_count)) %>%
          filter(!is.na(ROI_label)) # is.na ir true, not neg
        write.table(meta_data,
                    file = file.path(output_dir, "NicheDetect_Lasso_points.txt"),
                    sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

        # plot
        pdf_path <- file.path(output_dir, "NicheDetect_Lasso_plot.pdf")
        pdf(pdf_path,width = 12,height = 10)
        print(make_plot())
        dev.off()
        png_path <- file.path(output_dir, "NicheDetect_Lasso_plot.png")
        png(png_path, width = 12, height = 10, units = "in",res = 90)
        print(make_plot())
        dev.off()

        # summary
        summary_path <- file.path(output_dir, "NicheDetect_Lasso_summary.txt")
        summary_text <- c(
          "Lasso Analysis Summary",
          "======================",
          paste("Number of polygons:", length(values$polygons)),
          paste("Total ROI coords:", nrow(polygon_df)),
          paste("Total ROI points:", nrow(meta_data)),
          paste("Generated on:", Sys.time())
        )
        writeLines(summary_text, summary_path)
        zip::zip(
          zipfile = file,
          files = c("NicheDetect_Lasso_ROIcoords.txt", "NicheDetect_Lasso_points.txt",
                    "NicheDetect_Lasso_plot.pdf","NicheDetect_Lasso_plot.png", "NicheDetect_Lasso_summary.txt"),  # 只写文件名
          root = output_dir  # 文件所在目录指定根目录
        )

        file.rename(pdf_path, file.path(photo_dir, "NicheDetect_Lasso_plot.pdf"))
        file.rename(png_path, file.path(photo_dir, "NicheDetect_Lasso_plot.png"))
      }
    )
  }
  shinyApp(ui, server)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# .CalROIedge and .CalROIcenter
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Calculate ROI edge spots based on k-nearest neighbors
#'
#' Identifies spots on the edge of ROIs by analyzing the labels of neighboring
#' spots. A spot is considered an edge if it has a certain number of neighbors
#' belonging to different ROIs.
#'
#' @param meta_data Data frame containing spatial coordinates and ROI labels
#' @param ROI_label_colnm Character, column name containing ROI labels
#' @param ROI_edge_colnm Character, column name to store edge status
#' @param k Integer, number of nearest neighbors to consider
#' @param edge_thres Integer, minimum number of different-ROI neighbors to be
#'        considered an edge spot
#'
#' @return Modified data frame with added ROI_edge column
#'
#' @import FNN
#'
#' @keywords internal
#'
#' @noRd
.CalROIedge <- function(meta_data = NULL,
                        ROI_label_colnm = NULL,ROI_edge_colnm = NULL,
                        k = NULL, edge_thres = NULL){
  clog_normal(paste0("Calculate ROI edges ..."))
  coords <- as.matrix(meta_data[, c("x", "y")])
  roi_labels <- meta_data[[ROI_label_colnm]]
  knn_result <- get.knn(coords, k = k)
  neighbor_indices <- knn_result$nn.index
  diff_roi_count <- sapply(1:nrow(neighbor_indices), function(i) { # 计算每个点周围不同ROI的邻居数量
    i_indices <- neighbor_indices[i, ]
    neighbor_rois <- roi_labels[i_indices]
    sum( neighbor_rois != roi_labels[i], na.rm = TRUE)  # 不同ROI的数量
  })
  is_edge <- (roi_labels != "neg") & (diff_roi_count >= edge_thres)
  meta_data[[ROI_edge_colnm]] <- is_edge
  clog_normal("ROI edge label distribution:")
  print(
    table(meta_data[c(ROI_edge_colnm,ROI_label_colnm)])
  )
  return(meta_data)
}

#' Calculate ROI center spots
#'
#' Identifies the central spot for each ROI based on median coordinates and
#' distance to centroid. May adjust coordinates if the nearest spot is too far
#' from the centroid.
#'
#' @param meta_data Data frame containing spatial coordinates and ROI labels
#' @param ROI_label_colnm Character, column name containing ROI labels
#' @param ROI_center_colnm Character, column name to store center status
#' @param binsize Numeric, bin size of the spatial transcriptomics data
#'
#' @return Modified data frame with added ROI_center column and potentially
#'         adjusted coordinates for center spots
#'
#' @importFrom proxy dist
#' @importFrom purrr map map2 map_dfr list_rbind
#'
#' @keywords internal
#'
#' @noRd
.CalROIcenter <- function(meta_data = NULL,ROI_label_colnm = NULL, ROI_center_colnm = NULL,
                          coord_interval = NULL){
  clog_normal("Calculate ROI centers...")
  # browser()
  ROI_labels <- meta_data[[ROI_label_colnm]] %>% na.omit() %>% unique()
  ROI_center <- meta_data %>%
    filter(!!sym(ROI_label_colnm) != "neg") %>%
    mutate(cell_id = rownames(.)) %>%
    group_split(!!sym(ROI_label_colnm)) %>%
    map(function(df){
      roi_label <- df[[ROI_label_colnm]] %>% unique()
      roi_coord <- df[,c("x","y")]
      centroid_x <- median(df[["x"]], na.rm = TRUE) # not mean
      centroid_y <- median(df[["y"]], na.rm = TRUE)
      centroid_coord <-  data.frame(centroid_x = centroid_x,
                                    centroid_y = centroid_y,
                                    row.names = roi_label)
      elu.dist <- proxy::dist(roi_coord, centroid_coord) %>% # 这个耗时比较少，可以不替换
        matrix(byrow = F, nrow = nrow(roi_coord))
      rownames(elu.dist) <- df$cell_id
      colnames(elu.dist) <- roi_label
      dist_vec <- elu.dist[,1]
      min_dist <- min(dist_vec)
      cellnm_center <- rownames(elu.dist)[which.min(dist_vec)]
      return_data <- c("cell_nm" = cellnm_center, "roi_label" = roi_label,
                       "min_dist" = min_dist,"centroid_x" = round(centroid_x),"centroid_y" = round(centroid_y))
      return(return_data)
    })
  res_data <- bind_rows(ROI_center) %>%
    mutate(across(c(min_dist,centroid_x,centroid_y),as.numeric))
  meta_data <- mutate(meta_data,!!ROI_center_colnm := FALSE,is_adjcoord = FALSE)
  for(i in 1:nrow(res_data)){
    roi_label <- res_data$roi_label[i]
    min_dist <- res_data$min_dist[i]
    cell_nm <- res_data$cell_nm[i]
    if(min_dist>5*coord_interval){
      clog_warn(paste0("The min_dist to center for ROI '",roi_label,"' is ",round(min_dist,2),
                       ", which is lt 5 times of coord_interval (",coord_interval,"). The center coord wil be adjusted."))
      centroid_x <- res_data$centroid_x[i]
      centroid_y <- res_data$centroid_y[i]
      clog_warn(paste0("The coord before adjustment: (",meta_data[cell_nm,"x"],",",meta_data[cell_nm,"y"],")"))
      clog_warn(paste0("The coord after adjustment: (",centroid_x,",",centroid_y,")"))
      meta_data[cell_nm,c(ROI_center_colnm,"is_adjcoord")] <-  c(TRUE,TRUE)
      meta_data[cell_nm,c("x","y")] <-  c(centroid_x,centroid_y)
    }else{
      meta_data[cell_nm,ROI_center_colnm] <- TRUE
    }
  }
  return(meta_data)
}

#' Calculate distances to ROI centers
#'
#' Computes the distance from each spot to the nearest ROI center. Supports
#' two calculation modes: ROI-specific or general spot-based.
#'
#' @param meta_data Data frame containing spatial coordinates and ROI information
#' @param ROI_label_colnm Character, column name containing ROI labels
#' @param ROI_center_colnm Character, column name indicating center spots
#' @param CalMode Character, calculation mode - "ROI" for ROI-specific distances,
#'        "Spot" for general spot-based distances (default: "ROI")
#'
#' @return Modified data frame with added columns:
#'         \item{All_centernm}{Name of nearest center spot}
#'         \item{All_ROI_label}{ROI label of nearest center}
#'         \item{All_Dist2ROIcenter}{Distance to nearest ROI center}
#'         \item{ROI_region}{Categorized region type (neg/ROIs/edge/center)}
#'
#' @import FNN
#' @import dplyr
#'
#' @keywords internal
#'
#' @noRd
.CalDist2Center <- function(meta_data = NULL,ROI_label_colnm = NULL,ROI_center_colnm = NULL,
                            CalMode = "ROI"){

  clog_normal("Calculate distance to ROI centers...")
  ROI_centers <- meta_data %>%
    filter(!!sym(ROI_center_colnm) == TRUE) %>%
    dplyr::select(x, y)
  all_coord <- meta_data[,c("x","y")]
  clog_normal(paste0("Calculating distance matrix : ", nrow(all_coord), " X ", nrow(ROI_centers)))
  center2ROIlabel <- meta_data %>%
    filter(!!sym(ROI_center_colnm) == TRUE) %>%
    dplyr::select( !!sym(ROI_label_colnm))
  if(CalMode == "ROI"){
    # 较慢，但是对于ROI必须使用这个
    elu.dist <- proxy::dist(all_coord, ROI_centers) %>%
      matrix(byrow = F, nrow = nrow(all_coord))
    clog_normal("Distance matrix calculated.")
    rownames(elu.dist) <- rownames(meta_data)
    colnames(elu.dist) <- rownames(ROI_centers)
    min_dist_idx <- max.col(-elu.dist, ties.method = "first")
    meta_data <- meta_data %>%
      mutate(All_centernm = rownames(ROI_centers)[min_dist_idx], # ROI模式需要调整
             All_ROI_label = center2ROIlabel[All_centernm,ROI_label_colnm],
             All_ROI_label = ifelse(!!sym(ROI_label_colnm) != "neg", !!sym(ROI_label_colnm),All_ROI_label), # 防止聚类标签的点由于距离被调整了
             All_centernm = ifelse(!!sym(ROI_label_colnm) != "neg", rownames(center2ROIlabel)[match(All_ROI_label,center2ROIlabel[,ROI_label_colnm])],All_centernm),
             All_Dist2ROIcenter = elu.dist[as.matrix(data.frame(rownames(meta_data),All_centernm))]
      )

  }else if(CalMode == "Spot"){
    result_knnx <- get.knnx(data = ROI_centers, query = all_coord, k = 1)
    min_dist_idx <- result_knnx$nn.index[,1]
    min_dist_values <- result_knnx$nn.dist[,1]
    meta_data <- meta_data %>%
      mutate(All_centernm = rownames(ROI_centers)[min_dist_idx], # Spot模式不需要调整
             All_ROI_label = center2ROIlabel[All_centernm,ROI_label_colnm],
             All_Dist2ROIcenter = min_dist_values)
  }else{
    clog_error("CalMode must be 'ROI' or 'Spot'")
  }

  # >
  meta_data <- meta_data %>%
    mutate(ROI_region = if_else(ROI_label == "neg" ,"neg","ROIs"),
           .after = "ROI_label") %>%
    mutate(ROI_region = if_else(ROI_edge,"edge",ROI_region,missing = ROI_region),
           ROI_region = if_else(ROI_center,"center",ROI_region,missing = ROI_region),
           All_ROI_label2 = if_else(ROI_label == "neg",paste0("light_",All_ROI_label),All_ROI_label),.after = "All_ROI_label") %>%
    mutate(ROI_region = factor(ROI_region,levels = c("neg","ROIs","edge","center")))
  return(meta_data)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NicheDetect_STS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Automated ROI detection using density-based spatial clustering
#'
#' This function automatically detects regions of interest (ROIs) using
#' density-based clustering (DBSCAN) on spatial transcriptomics data. It includes
#' preprocessing steps like density filtering and iterative spot updating.
#'
#' @param STID_obj An STID object containing spatial transcriptomics data
#' @param meta_key Character, metadata key containing positive spot information
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param spatial_scale_method Character, detection level - "region" or "spot" (default: "region")
#' @param region_detect_method Character, hull method - "convex" or "concave" (default: "convex")
#' @param concavity Numeric, concavity parameter for concave hull (default: 2)
#' @param update_spots Logical, whether to iteratively update positive spots
#'        (default: TRUE)
#' @param pos_colnm Character, column name containing positive spot labels
#' @param neg_value Character, value indicating negative/non-ROI spots (default: "neg")
#' @param density_thres Numeric, density threshold for filtering (0-1, default: 0.9)
#' @param ROI_size Integer, minimum number of spots to form an ROI (default: 10)
#' @param minPts Integer, minimum points parameter for DBSCAN (default: NULL,
#'        automatically determined based on data format)
#' @param k_kNNdist Integer, k value for kNN distance calculation (default: NULL,
#'        uses minPts if not specified)
#' @param description Character, description of the analysis (default: NULL)
#' @param grp_nm Character, group name for output organization (default: NULL, uses timestamp)
#' @param dir_nm Character, directory name for output (default: "M2_NicheDetect_STS")
#'
#' @return Returns the modified STID object with added metadata containing ROI labels,
#'         distances to ROI centers, edge information, and region classifications
#'
#' @import dbscan
#' @import Seurat
#' @import ggdensity
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Automatically detect ROIs based on positive spots
#' STID_obj <- NicheDetect_STS(
#'   STID_obj = ist_object,
#'   meta_key = "M1_SpotDetect_Gene_20240101",
#'   pos_colnm = "Label_geneA",
#'   density_thres = 0.9,
#'   ROI_size = 10
#' )
#' }
NicheDetect_STS <- function(STID_obj = NULL, meta_key = NULL,
                            loop_id = "LoopAllSamp",
                            pos_colnm = NULL, neg_value = "neg",
                            spatial_scale_method = "region", region_detect_method = "convex", # convex or concave
                            concavity = 2,
                            update_spots = TRUE,
                            density_thres = 0.9, ROI_size = 10,
                            minPts = NULL, k_kNNdist = NULL,
                            description = NULL,
                            grp_nm = NULL, dir_nm = "M2_NicheDetect_STS"
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
  .check_null_args(meta_key, pos_colnm)
  match.arg(region_detect_method, c('convex', 'concave'))
  match.arg(spatial_scale_method, c("region", "spot"))

  # >
  loop_single <- .check_loop_single(STID_obj = STID_obj,loop_id = loop_id, mode = 1)
  all_single <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "samp_id")[[1]]
  if(density_thres < 0 | density_thres > 1){
    clog_error("density_thres must be between 0 and 1")
  }
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  meta_data <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]]
  .check_column_exist(meta_data, pos_colnm)
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  x_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "x_colnm")[[1]]
  y_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "y_colnm")[[1]]
  data_format <- GetInfo(STID_obj, info_key = "data_info",sub_key = "data_format")[[1]]
  data_platform <- GetInfo(STID_obj, info_key = "data_info",sub_key = "data_platform")[[1]]
  binsize <- GetInfo(STID_obj, info_key = "data_info",sub_key = "binsize")[[1]]
  coord_interval <- GetInfo(STID_obj, info_key = "data_info",sub_key = "coord_interval")[[1]]
  meta2STID_list <- list()
  n <- 1
  raw_ROI_size <- ROI_size
  for(i in seq_along(all_single)){
    i_single <- all_single[i]
    i_coord_interval <- coord_interval[i]
    i_meta_data <- meta_data %>%
      filter(!!sym(samp_colnm) == i_single) %>%
      dplyr::select(!!sym(samp_colnm), !!sym(x_colnm), !!sym(y_colnm), !!sym(pos_colnm))
    if(!i_single %in% loop_single){
      meta2STID_list[[i_single]] <- i_meta_data
      next
    }
    clog_loop(paste0("Processing samp_id: ", i_single, " (", n, "/", length(loop_single), ")"))
    n <- n+1
    i_output_dir <- paste0(output_dir,"/",i_single,"/")
    i_photo_dir <- paste0(photo_dir,"/",i_single,"/")
    dir.create(i_output_dir,recursive = TRUE,showWarnings = FALSE)
    dir.create(i_photo_dir,recursive = TRUE,showWarnings = FALSE)

    #> pos_meta_data
    pos_meta_data <- i_meta_data %>%
      filter(.,!!sym(pos_colnm) != neg_value)
    len_pos_value <- nrow(pos_meta_data)
    cells_pct <- len_pos_value/nrow(i_meta_data)
    clog_normal(paste0("Positive cells percentage: ",round(cells_pct*100,2),"% (",
                       len_pos_value,"/", nrow(i_meta_data),")"))
    if(is.null(raw_ROI_size)){
      ROI_size <- max(nrow(pos_meta_data)*0.02,10) # 2% of positive cells or 10, whichever is larger
      clog_normal(paste0("ROI_size is set to ", ROI_size))
    }
    if(len_pos_value < ROI_size){
      clog_warn(paste0("The number of positive spots (",len_pos_value,") is less than ROI_size (",ROI_size,") , skiping..."))
      meta2STID_list[[i_single]] <- i_meta_data
      next
    }

    #>>> STS_ROI
    clog_step("Preprocessing for STS_ROI detection ...")

    #> .remove_low_density
    clog_title(paste0("Remove low-density points ..."))
    pos_meta_data <- .remove_low_density(meta_data = pos_meta_data, # 过滤了一些阳性细胞数
                                         density_thres = density_thres,
                                         photo_dir = i_photo_dir,
                                         samp_id = i_single)
    if(nrow(pos_meta_data) < ROI_size){
      clog_warn(paste0("The number of positive spots after density filtering (",nrow(pos_meta_data),") is less than ROI_size (",ROI_size,") , skiping..."))
      meta2STID_list[[i_single]] <- i_meta_data
      next
    }

    #> .updateSigSpots
    if(data_format  == "square_grid"){
      k_update <- 8; neg2pos_thres <- 5; pos2neg_thres <- 1; edge_thres <- 1; diff_thres <- len_pos_value*0.001
    }else if(data_format == "hex_grid" & data_platform == "Visium"){
      k_update <- 6; neg2pos_thres <- 4; pos2neg_thres <- 1; edge_thres <- 1; diff_thres <- 0
    }

    if(update_spots){
      clog_title(paste0("Update positive spots iteratively ..."))
      n <- 1
      raw_meta_data <- i_meta_data
      i_meta_data <- i_meta_data %>% # processing the i_meta_data according to the pos_meta_data
        mutate(!!sym(pos_colnm) := if_else(rownames(.) %in% rownames(pos_meta_data), !!sym(pos_colnm), neg_value, missing = neg_value)) %>%
        mutate(pos_spot = if_else(!!sym(pos_colnm) != neg_value, 1, 0, missing = 0))
      pos_value <- i_meta_data %>% filter(pos_spot == 1) %>% pull(!!sym(pos_colnm)) %>% unique()
      while(TRUE){
        clog_normal(paste0("Iteration: ", n))
        new_meta_data <- .updateSigSpots(meta_data = i_meta_data,
                                         pos_colnm = "pos_spot",
                                         k = k_update,
                                         pos2neg_thres = pos2neg_thres,
                                         neg2pos_thres = neg2pos_thres)
        raw_rownm <- i_meta_data[i_meta_data[["pos_spot"]] == 1,] %>% rownames()
        new_rownm <- new_meta_data[new_meta_data[["pos_spot"]] == 1,] %>% rownames()
        len_diff <- length(setdiff(raw_rownm, new_rownm)) + length(setdiff(new_rownm, raw_rownm))
        if(len_diff <= diff_thres){
          clog_normal(paste0("Converged: changed positive cells = ", len_diff, " (diff thres = ", diff_thres, ")"))
          break
        }else{
          clog_normal(paste0("Not converged: changed positive cells = ", len_diff))
          i_meta_data <- new_meta_data
        }
        n <- n+1
      }
      new_meta_data <- new_meta_data %>%
        mutate(!!sym(pos_colnm) := if_else(pos_spot == 1, pos_value, neg_value, missing = neg_value))
      pos_meta_data <- new_meta_data %>%
        filter(.,!!sym(pos_colnm) != neg_value)
      p1 <- Plot_Spatial(plot_data =new_meta_data,x_colnm = "x",y_colnm = "y",group_by = pos_colnm,
                         datatype = "discrete",
                         col = list(dis = c("grey95","#EB1E2C"),con = NULL),
                         pt_size = 0.25,vmin = NULL, vmax = "p99",
                         title = "Updated positive spots",
                         subtitle = NULL,black_bg = F)
      ggsave(filename = paste0(i_photo_dir,"/",i_single,"_Updated_PosSpots(NicheDetect_STS).pdf"),plot = p1,width = 6,height = 5)
    }
    if(nrow(pos_meta_data) < ROI_size){
      clog_warn(paste0("The number of positive spots after iterative updating (",nrow(pos_meta_data),") is less than ROI_size (",ROI_size,") , skiping..."))
      meta2STID_list[[i_single]] <- i_meta_data
      next
    }

    #> spatial_scale_method
    if(spatial_scale_method == "region"){
      clog_title("Region level detection ...")
      if(data_format  == "square_grid"){
        # minPts <- 5; eps <- 3.9
        if(is.null(minPts)){
          minPts <- 5 # 8*0.6 = 5
        }
      }else if(data_format == "hex_grid" & data_platform == "Visium"){
        if(is.null(minPts)){
          minPts <- 4 # 6*0.6 = 4
        }
      }
      if(is.null(k_kNNdist)){
        k_kNNdist <- minPts
      }
      eps <- mean(kNNdist(pos_meta_data[c("x","y")], k = k_kNNdist))

      #>
      clog_normal("DBSCAN clustering for ROI detection ...")
      clog_normal(paste0("DBSCAN parameters: eps = ", round(eps,2), ", minPts = ", minPts))
      DBSCAN_meta_data <- .DBSCAN_clustering(meta_data = pos_meta_data, # 含有ID
                                             eps = eps,
                                             minPts = minPts,
                                             ROI_size = ROI_size,
                                             samp_id = i_single,
                                             photo_dir = i_photo_dir)
      if(region_detect_method == "convex"){
        clog_normal("Using convex hull for ROI detection ...")
        hull_data <- .calculate_hull(meta_data = DBSCAN_meta_data,
                                     ROI_colnm = "cluster",
                                     region_detect_method = "convex",
                                     samp_id = i_single,
                                     photo_dir = i_photo_dir)
      }else if(region_detect_method == "concave"){
        clog_normal("Using concave hull for ROI detection ...")
        hull_data <- .calculate_hull(meta_data = DBSCAN_meta_data,
                                     ROI_colnm = "cluster",
                                     region_detect_method = "concave",
                                     concavity = concavity,
                                     samp_id = i_single,
                                     photo_dir = i_photo_dir)
      }
      if(is.null(hull_data)){
        clog_warn("No ROIs detected after hull calculation, skipping...")
        meta2STID_list[[i_single]] <- i_meta_data
        next
      }
      write.table(hull_data, file = paste0(i_output_dir,"/",i_single,"_DBSCAN_",region_detect_method,"_(NicheDetect_STS).txt"),
                  sep = "\t",col.names = NA, row.names = TRUE, quote = FALSE)
      res_meta_data <- .point_in_hull(meta_data = i_meta_data,
                                      hull_data = hull_data,
                                      ROI_colnm = "cluster",
                                      samp_id = i_single,
                                      photo_dir = i_photo_dir)


    }else if(spatial_scale_method == "spot"){
      clog_title("Spot level detection ...")
      if(data_format  == "square_grid"){
        adj_matrix <- .as_adjMatrix(coords = pos_meta_data[c("x","y")], k_neighbors = 8, dist_thres = i_coord_interval*1.25)
      }else if(data_format == "hex_grid" & data_platform == "Visium"){
        adj_matrix <- .as_adjMatrix(coords = pos_meta_data[c("x","y")], k_neighbors = 6, dist_thres = i_coord_interval*1.25)
      }
      res_meta_data <- .CalAdjMatrixCluster(adj_matrix = adj_matrix, pos_meta_data = pos_meta_data, ROI_size = ROI_size)
      if(nrow(res_meta_data) == 0){
        clog_warn("No ROIs detected after adjacency matrix clustering, skipping...")
        meta2STID_list[[i_single]] <- i_meta_data
        next
      }
    }

    #> CalDist2Center
    common_ids <- intersect(rownames(i_meta_data), rownames(res_meta_data))
    i_meta_data$ROI_label   <- "neg" # NA -> "neg"
    i_meta_data$ROI_count   <- NA_integer_
    i_meta_data[common_ids, "ROI_label"]   <- res_meta_data[common_ids, "ROI_label"]
    i_meta_data[common_ids, "ROI_count"]   <- res_meta_data[common_ids, "ROI_count"]
    i_meta_data <- .CalROIedge(meta_data = i_meta_data,
                               ROI_label_colnm = "ROI_label",
                               ROI_edge_colnm = "ROI_edge",
                               k = k_update,edge_thres = edge_thres)
    i_meta_data <- .CalROIcenter(meta_data = i_meta_data,
                                 ROI_label_colnm = "ROI_label",
                                 ROI_center_colnm = "ROI_center",
                                 coord_interval = i_coord_interval)
    i_meta_data <- .CalDist2Center(meta_data = i_meta_data,
                                   ROI_label_colnm = "ROI_label",
                                   ROI_center_colnm = "ROI_center",
                                   CalMode = "ROI")
    meta2STID_list[[i_single]] <- i_meta_data
  }

  #> AddMetaData
  meta2STID <- bind_rows(meta2STID_list) %>%
    as.data.frame()
  meta2STID <- meta2STID[rownames(meta_data), , drop = FALSE] # keep order
  new_meta_key <- paste0(dir_nm,"_",grp_nm)
  STID_obj <- AddMetaData(STID_obj = STID_obj,
                         meta_key = new_meta_key,
                         add_data = meta2STID,
                         dir_nm = dir_nm,
                         grp_nm = grp_nm,
                         asso_key = meta_key,
                         description = description)

  # >>> Final
  .save_function_params("NicheDetect_STS", envir = environment(), file = paste0(output_dir,"Log_function_params_(NicheDetect_STS).log"))
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(NicheDetect_STS).log")) %>% invisible()
  return(STID_obj)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# STS_ROI: preprocessing
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Remove low-density points based on kernel density estimation
#'
#' Filters out points in low-density regions using kernel density estimation.
#' Points below the specified density threshold are removed.
#'
#' @param meta_data Data frame containing spatial coordinates
#' @param density_thres Numeric, density threshold (0-1) - points with density
#'        below this quantile are removed
#' @param samp_id Character, sample identifier for plot titles
#' @param photo_dir Character, directory path for saving density plots
#'
#' @return Filtered data frame containing only high-density points
#'
#' @import ggdensity
#' @import ggplot2
#'
#' @keywords internal
#'
#' @noRd
.remove_low_density <- function(meta_data = NULL,density_thres = NULL,
                                samp_id = NULL,photo_dir = NULL){
  p1 <- ggplot(meta_data,aes(x,y))+
    theme_bw() +
    geom_point(shape=16,size=0.5,color = "#CAB2D6",stroke = 0,alpha = 0.6)+
    geom_hdr_lines(aes(color = after_stat(probs)), # fter_stat(probs)
                   method = "kde", # "kde","histogram", "freqpoly"，kde默认是高斯核函数
                   probs = c(0.9,0.7,0.5,0.3,0.1),
                   alpha = 1,
                   linewidth = 0.75 # 对应color的宽度
    ) +
    scale_color_manual(values = c('#0c3383', '#0a88ba', '#f2d338', '#f28f38', '#d91e1e')) +
    labs(x = "Column", y = "Row", title = samp_id) +
    coord_fixed(ratio = 1) +
    scale_y_reverse() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_rect(fill="white",color = "white")
    )
  ggsave(filename = paste0(photo_dir,"/",samp_id,"_PosDensity_(NicheDetect_STS).pdf"),
         plot = p1, width = 6, height = 5)

  # >
  res_df <- get_hdr(meta_data[c("x","y")], method = "kde",
                    n = 200,
                    # probs = c(0.9,0.7,0.5,0.3,0.1),
                    probs = seq(0.05, 0.95, by = 0.05) %>% rev(),
                    hdr_membership = TRUE)
  res_data <- res_df$data %>%
    filter(hdr_membership <= density_thres)
  table(res_data$hdr_membership)
  meta_data <- meta_data[rownames(res_data),]
  return(meta_data)
}

#' Iteratively update positive spot labels based on neighbor information
#'
#' Uses k-nearest neighbors to update positive spot labels in two steps:
#' 1. Negative to positive: spots with many positive neighbors become positive
#' 2. Positive to negative: positive spots with few positive neighbors become negative
#'
#' @param meta_data Data frame containing spatial coordinates and spot labels
#' @param pos_colnm Character, column name containing positive spot indicators (0/1)
#' @param k Integer, number of nearest neighbors to consider
#' @param pos2neg_thres Integer, threshold for positive to negative conversion
#'        (positive spots with <= this many positive neighbors become negative)
#' @param neg2pos_thres Integer, threshold for negative to positive conversion
#'        (negative spots with >= this many positive neighbors become positive)
#'
#' @return Modified data frame with updated positive spot labels
#'
#' @import FNN
#'
#' @keywords internal
#'
#' @noRd
.updateSigSpots <- function(meta_data = NULL, pos_colnm = NULL, k = NULL,
                            pos2neg_thres = NULL, neg2pos_thres = NULL) {
  clog_normal("Finding positive points using FNN method...")
  raw_meta_data <- meta_data
  coords <- as.matrix(meta_data[, c("x", "y")])
  current_status <- meta_data[[pos_colnm]]

  #> clog_normal("Step 1: Negative → Positive")
  knn_result1 <- get.knn(coords, k = k)
  neighbor_indices1 <- knn_result1$nn.index
  positive_count1 <- apply(neighbor_indices1, 1, function(indices) { # 每个点周围pos点数
    sum(meta_data[[pos_colnm]][indices])
  })
  neg2pos <- current_status == 0 & positive_count1 >= neg2pos_thres
  current_status[neg2pos] <- 1

  #> clog_normal("Step 2: Positive → Negative")
  knn_result2 <- get.knn(coords, k = k)
  neighbor_indices2 <- knn_result2$nn.index
  positive_count2 <- apply(neighbor_indices2, 1, function(indices) {
    sum(current_status[indices])  # 使用第一步更新后的状态
  })
  pos2neg <- current_status == 1 & positive_count2 <= pos2neg_thres
  current_status[pos2neg] <- 0

  # >
  meta_data[[pos_colnm]] <- current_status
  # meta_data$neighbors_count_step1 <- positive_count1
  # meta_data$neighbors_count_step2 <- positive_count2
  # meta_data$converted_neg2pos <- neg2pos
  # meta_data$converted_pos2neg <- pos2neg
  clog_normal(paste("Raw pos:", sum(raw_meta_data[[pos_colnm]] == 1), " neg → pos:", sum(neg2pos),
                    " pos → neg:", sum(pos2neg)," Final pos:", sum(current_status == 1)))
  return(meta_data)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# STS_ROI: Region_level
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Perform DBSCAN clustering on spatial coordinates
#'
#' Applies DBSCAN clustering to identify ROIs based on spatial density.
#' Includes visualization of k-distance plot and clustering results.
#'
#' @param meta_data Data frame containing spatial coordinates
#' @param eps Numeric, DBSCAN epsilon parameter (maximum distance between points)
#' @param minPts Integer, DBSCAN minimum points parameter
#' @param ROI_size Integer, minimum cluster size to be considered an ROI
#' @param samp_id Character, sample identifier for plot titles
#' @param photo_dir Character, directory path for saving plots
#'
#' @return Data frame with added cluster labels, filtered for minimum ROI size
#'
#' @import dbscan
#' @import patchwork
#'
#' @keywords internal
#'
#' @noRd
.DBSCAN_clustering <- function(meta_data = NULL,eps = NULL, minPts = NULL, ROI_size = NULL,
                               samp_id = NULL,photo_dir = NULL) {
  clog_normal(paste0("DBSCAN clustering ..."))

  #> kNNdistplot
  pdf(file = paste0(photo_dir,"/",samp_id,"_DBSCAN_kNNdistplot_(NicheDetect_STS).pdf"),
      width = 6,height = 5)
  kNNdistplot(meta_data[c("x","y")], k = minPts - 1)
  abline(h = eps, col = "red", lty = 2)
  dev.off()

  #> dbscan
  res_df <- dbscan(meta_data[c("x","y")],
                   eps = eps,
                   minPts = minPts)
  res_df$cluster <- paste0("ROI_",res_df$cluster)
  meta_data <- cbind(meta_data, cluster = res_df$cluster)
  meta_data_fil <- meta_data %>%
    group_by(cluster) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    filter(cluster != "ROI_0") %>%  # 去除噪声点
    filter(n >= ROI_size) # 去除噪声点
  clog_normal(paste0("The Results of DBSCAN clustering after filtering:"))
  print(
    table(meta_data_fil["cluster"])
  )

  #> plot
  p1 <- ggplot(meta_data, aes(x, y, color = cluster)) +
    geom_point(size = 1) +
    labs(x = "Column", y = "Row", title = samp_id) +
    theme_bw() +
    coord_fixed(ratio = 1) +
    scale_y_reverse() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_rect(fill="white",color = "white")
    )
  p2 <- ggplot(meta_data_fil, aes(x, y, color = cluster)) +
    geom_point(size = 1) +
    labs(x = "Column", y = "Row", title = samp_id) +
    theme_bw() +
    coord_fixed(ratio = 1) +
    scale_y_reverse() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_rect(fill="white",color = "white")
    )
  pdf(file = paste0(photo_dir,"/",samp_id,"_DBSCAN_ClusteringPlot_(NicheDetect_STS).pdf"),
      width = 8,height = 10)
  print(p1/p2)
  dev.off()

  return(meta_data_fil)
}

#' Calculate convex or concave hull for ROI clusters
#'
#' Computes hull polygons around clustered points to define ROI boundaries.
#' Supports both convex and concave hull methods.
#'
#' @param meta_data Data frame containing clustered points
#' @param ROI_colnm Character, column name containing cluster/ROI labels
#' @param region_detect_method Character, hull method - "convex" or "concave"
#' @param concavity Numeric, concavity parameter for concave hull (only used if region_detect_method is "concave")
#' @param samp_id Character, sample identifier for plot titles
#' @param photo_dir Character, directory path for saving plots
#'
#' @return Data frame containing hull polygon coordinates for each ROI
#'
#' @import concaveman
#' @import ggplot2
#'
#' @keywords internal
#'
#' @noRd
.calculate_hull <- function(meta_data = NULL,ROI_colnm = NULL,
                            region_detect_method = NULL, concavity = NULL,
                            samp_id = NULL,photo_dir = NULL) {
  meta_data <- meta_data[c("x","y",ROI_colnm)]
  hull_list <- list()
  for (cl in unique(meta_data[[ROI_colnm]])) {
    pts <- meta_data[meta_data[[ROI_colnm]] == cl, ]
    if (nrow(pts) >= 3) {
      if(region_detect_method == "convex"){
        idx <- chull(pts$x, pts$y)
        hull_pts <- pts[idx, , drop = FALSE]
      }else if(region_detect_method == "concave"){
        pts_matrix <- as.matrix(pts[c("x", "y")])
        hull_pts <- concaveman(pts_matrix, concavity = concavity, length_threshold = 0) %>%
          as.data.frame()
        colnames(hull_pts) <- c("x", "y")
        hull_pts[[ROI_colnm]] <- cl
      }
      hull_list[[as.character(cl)]] <- hull_pts
    }
  }
  if (length(hull_list) > 0) {
    hull_data <- bind_rows(hull_list)
  } else {
    clog_warn(paste0("No cluster has enough points to compute ",region_detect_method," hull."))
    return(NULL)
  }

  p1 <- ggplot() +
    geom_polygon(
      data = hull_data,
      aes(x = x, y = y, group = !!sym(ROI_colnm), fill = !!sym(ROI_colnm)),
      alpha = 0.25,
      color = "black",
      linewidth = 0.6
    ) +
    geom_point(
      data = meta_data,
      aes(x = x, y = y, color = !!sym(ROI_colnm)),
      size = 1
    ) +
    labs(x = "Column", y = "Row", title = samp_id) +
    theme_bw() +
    coord_fixed(ratio = 1) +
    scale_y_reverse() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_rect(fill="white",color = "white")
    )
  ggsave(filename = paste0(photo_dir,"/",samp_id,"_DBSCAN_",region_detect_method,"_(NicheDetect_STS).pdf"),
         plot = p1, width =8, height = 8)

  return(hull_data)
}

#' Determine which spots fall within ROI hull polygons
#'
#' Uses point-in-polygon algorithm to identify spots belonging to each ROI.
#' Handles spots that may belong to multiple ROIs.
#'
#' @param meta_data Data frame containing all spots
#' @param hull_data Data frame containing hull polygon coordinates
#' @param ROI_colnm Character, column name containing ROI labels in hull data
#' @param samp_id Character, sample identifier for plot titles
#' @param photo_dir Character, directory path for saving plots
#'
#' @return Filtered data frame containing spots that belong to ROIs, with
#'         added ROI_label and ROI_count columns
#'
#' @importFrom sp point.in.polygon
#'
#' @keywords internal
#'
#' @noRd
.point_in_hull <- function(meta_data = NULL,
                           hull_data = NULL,
                           ROI_colnm = NULL,
                           samp_id = NULL,photo_dir = NULL) {
  clog_normal(paste0("Point in hull detection ..."))

  # >
  spot_name <- rownames(meta_data)
  x_all <- meta_data$x
  y_all <- meta_data$y
  polygon_list <- split(hull_data, hull_data[ROI_colnm])
  list_inside <- lapply(polygon_list, function(poly) {
    px <- poly$x # 提取多边形 x, y 坐标
    py <- poly$y
    inside <- point.in.polygon(x_all, y_all, px, py) # 0=outside, 1,2=inside/border,3=vertex
    return(inside)
  })
  logic_inside <- lapply(list_inside, function(x) {
    return(x > 0)
  })
  logic_edge <- lapply(list_inside, function(x) {
    return(x > 1) # TRUE if inside，根据逻辑值取点
  })
  logic_inside_matrix <- bind_cols(logic_inside)
  logic_edge_matrix <- bind_cols(logic_edge)
  colnm <- colnames(logic_inside_matrix) # 多边形名字
  ROI_label <- apply(logic_inside_matrix, 1, function(row) {
    included_rois <- colnm[which(row)]
    if (length(included_rois) == 0) {
      return(NA_character_)
    } else {
      return(paste(included_rois, collapse = "; "))
    }
  })
  ROI_count <- rowSums(logic_inside_matrix)
  ROI_edge <- rowSums(logic_edge_matrix) > 0 # 这个边界是多边形的边界

  # >
  meta_data <- meta_data %>%
    # cbind(data.frame(ROI_label,ROI_count,ROI_edge)) %>%
    cbind(data.frame(ROI_label,ROI_count)) %>%
    filter(!is.na(ROI_label)) # is.na is true, not neg
  clog_normal(paste0("The Results of Point in hull detection:"))
  print(
    table(meta_data["ROI_label"])
  )
  if(max(meta_data$ROI_count) >1){
    clog_warn("Detected points belong to multiple ROIs, these points will be removed in final results.")
  }
  meta_data_fil <- meta_data %>%
    filter(ROI_count <2)

  p1 <- ggplot(data = meta_data,
               aes(x = x, y = y, color = ROI_label)) +
    geom_polygon(
      data = hull_data, aes(x = x, y = y, group = !!sym(ROI_colnm)),
      alpha = 0.25,color = "black",fill = NA,
      linewidth = 0.25
    ) +
    geom_point(size = 0.25,shape = 16,stroke = 0
    ) +
    labs(x = "Column", y = "Row", title = samp_id) +
    theme_bw() +
    coord_fixed(ratio = 1) +
    scale_y_reverse() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_rect(fill="white",color = "white")
    )
  p2 <- ggplot(data = meta_data_fil,
               aes(x = x, y = y, color = ROI_label)) +
    geom_polygon(
      data = hull_data, aes(x = x, y = y, group = !!sym(ROI_colnm)),
      alpha = 0.25,color = "black",fill = NA,
      linewidth = 0.25
    ) +
    geom_point(size = 0.25,shape = 16,stroke = 0
    ) +
    labs(x = "Column", y = "Row", title = samp_id) +
    theme_bw() +
    coord_fixed(ratio = 1) +
    scale_y_reverse() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_rect(fill="white",color = "white")
    )
  pdf(file = paste0(photo_dir,"/",samp_id,"_DBSCAN_PointInHull_(NicheDetect_STS).pdf"),
      width = 12,height = 6)
  print(p1 + p2)
  dev.off()
  return(meta_data_fil)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# STS_ROI: Spot_level
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Convert spatial coordinates to adjacency matrix
#'
#' Constructs a sparse adjacency matrix from spatial coordinates using k-nearest
#' neighbors and distance threshold.
#'
#' @param coords Data frame containing 'x' and 'y' columns with spatial coordinates
#' @param k_neighbors Integer, number of nearest neighbors to consider
#' @param dist_thres Numeric, distance threshold for adjacency (default: uses
#'        median of first neighbor distances)
#'
#' @return A sparse adjacency matrix (dgCMatrix) where TRUE indicates neighboring spots
#'
#' @import FNN Matrix
#'
#' @keywords internal
#' @noRd
.as_adjMatrix <- function(coords = NULL, k_neighbors = NULL, dist_thres = NULL ) {
  if(ncol(coords) < 2 || !all(c("x", "y") %in% colnames(coords))){
    clog_error("Coordinate data must contain 'x' and 'y' columns.")
  }
  coords <- coords[, c("x", "y"), drop = FALSE]
  n_cells <- nrow(coords)
  coords_mat <- as.matrix(coords)
  knn_result <- get.knn(coords_mat, k = k_neighbors)
  if(is.null(dist_thres)){
    dist_thres <- quantile(knn_result$nn.dist[, 1], 0.5)
  }
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
  adj_matrix <- adj_matrix | t(adj_matrix)
  diag(adj_matrix) <- FALSE # 对称化矩阵，去掉自环
  adj_matrix <- as(adj_matrix, "dgCMatrix")
  return(adj_matrix)
}


#' Calculate clusters from adjacency matrix
#'
#' Identifies connected components (clusters/ROIs) from an adjacency matrix
#' and filters by minimum cluster size.
#'
#' @param adj_matrix Sparse adjacency matrix (dgCMatrix) from `.as_adjMatrix`
#' @param pos_meta_data Data frame containing spatial metadata with rownames as cell IDs
#' @param ROI_size Integer, minimum number of spots required for a valid ROI
#'
#' @return Filtered metadata with added 'ROI_label' and 'ROI_count' columns,
#'         only retaining spots in clusters meeting the size threshold
#'
#' @import dplyr tibble
#' @importFrom igraph graph_from_adjacency_matrix
#'
#' @keywords internal
#' @noRd
.CalAdjMatrixCluster <- function(adj_matrix = NULL,pos_meta_data = NULL,ROI_size = NULL){
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected") # !!! time-consuming step
  comps <- igraph::components(g)$membership
  pos_meta_data$ROI_label <- paste0("ROI_", comps)
  pos_meta_data <- pos_meta_data %>%
    rownames_to_column("cell_id") %>%
    group_by(ROI_label) %>%
    mutate(ROI_count = n()) %>%
    ungroup() %>%
    filter(ROI_count >= ROI_size) %>%
    column_to_rownames("cell_id")

  return(pos_meta_data)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  CompareNiche
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Compare two niche detection results
#'
#' Compares two sets of niche/ROI detection results using Venn diagrams and
#' distance analysis. Identifies overlapping spots and calculates distances
#' between niche centers.
#'
#' @param STID_obj An STID object containing niche detection results
#' @param meta_key1 Character, first metadata key for comparison
#' @param meta_key2 Character, second metadata key for comparison
#' @param bins Integer, number of bins for distance distribution plots (default: 15)
#'
#' @return NULL (invisible), generates plots and prints comparison statistics
#'
#' @import ggplot2
#' @import FNN
#' @import ggvenn ggvenn
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare two niche detection methods
#' CompareNiche(
#'   STID_obj = ist_object,
#'   meta_key1 = "M2_NicheDetect_Lasso_20240101",
#'   meta_key2 = "M2_NicheDetect_STS_20240101"
#' )
#' }
CompareNiche <- function(STID_obj = NULL,
                         meta_key1 = NULL,
                         meta_key2 = NULL,
                         bins = 15
                         ) {

  # >>> Start pipeline
  clog_start()
  if (!inherits(STID_obj, "STID")) {
    clog_error("Input object is not an STID object.")
  }
  .check_null_args(meta_key1, meta_key2)
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info", sub_key = "samp_colnm")[[1]]
  x_colnm <- GetInfo(STID_obj, info_key = "data_info", sub_key = "x_colnm")[[1]]
  y_colnm <- GetInfo(STID_obj, info_key = "data_info", sub_key = "y_colnm")[[1]]

  # >>> Start main pipeline
  meta_data1 <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key1)[[1]]
  meta_data2 <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key2)[[1]]
  meta_data1_Niche <- meta_data1 %>% filter(ROI_label != "neg")
  meta_data2_Niche <- meta_data2 %>% filter(ROI_label != "neg")
  meta_data1_Niche_coord <- meta_data1_Niche %>% dplyr::select(x,y)
  meta_data2_Niche_coord <- meta_data2_Niche %>% dplyr::select(x,y)
  meta_data1_cells <- rownames(meta_data1_Niche)
  meta_data2_cells <- rownames(meta_data2_Niche)

  #> ggVennDiagram
  clog_normal("Plotting Venn diagram ...")
  venn_list <- list(Niche1 = meta_data1_cells, Niche2 = meta_data2_cells)
  p1 <- ggvenn(
    data = venn_list,
    show_elements = F,
    label_sep = "\n",
    show_percentage = T,
    digits = 2,
    fill_color = c("#E41A1C", "#1E90FF"),
    fill_alpha = 0.4,
    stroke_color = "white",
    stroke_alpha = 0.5,
    stroke_size = 0.75,
    stroke_linetype = "solid",
    set_name_color = "black",
    set_name_size = 6,
    text_color = "black",
    text_size = 5
  )
  print(p1)

  #> knnx
  clog_normal("Calculating distance to nearest niche center ...")
  meta_data1_center <- meta_data1_Niche %>% filter(ROI_center) %>%
    dplyr::select(x,y)
  meta_data2_center <- meta_data2_Niche %>% filter(ROI_center) %>%
    dplyr::select(x,y)
  meta_data1_edge <- meta_data1_Niche %>% filter(ROI_edge) %>%
    pull(All_Dist2ROIcenter) %>%
    median(na.rm = TRUE)
  meta_data2_edge <- meta_data2_Niche %>% filter(ROI_edge) %>%
    pull(All_Dist2ROIcenter) %>%
    median(na.rm = TRUE)
  res_knnx1 <- get.knnx(data = meta_data1_center,
                        query = meta_data2_Niche_coord, # 返回行数等于coord_all，有k列，值是coord_pos
                        k = 1,algorithm = "kd_tree"
  )
  res_knnx2 <- get.knnx(data = meta_data2_center,
                        query = meta_data1_Niche_coord, # 返回行数等于coord_all，有k列，值是coord_pos
                        k = 1,algorithm = "kd_tree"
  )
  knnx1_dist <- data.frame(dist = res_knnx1$nn.dist[,1], grp = "Niche2_to_Niche1_center")
  knnx2_dist <- data.frame(dist = res_knnx2$nn.dist[,1], grp = "Niche1_to_Niche2_center")

  #>
  clog_normal("Plotting distance distribution ...")
  plot_data <- bind_rows(knnx1_dist, knnx2_dist)
  plot_data_edge <- data.frame(
    dist = c(meta_data1_edge, meta_data2_edge),
    grp = c("Niche2_to_Niche1_center", "Niche1_to_Niche2_center")
  )
  clog_normal(paste0("Niche1 edge distance: ", round(meta_data1_edge,2)))
  clog_normal(paste0("Niche2 edge distance: ", round(meta_data2_edge,2)))
  p1 <- ggplot(plot_data, aes( x = dist, fill = grp,color = grp)) +
    # geom_density( alpha = 0.5,bw = 10,trim = T)
    geom_histogram(aes(y = after_stat(count)), bins = bins, alpha = 0.8,color = "grey98",linewidth = 0.25) +
    geom_text(aes(label=ifelse(after_stat(count) == 0, "", after_stat(count)),
                  y=after_stat(count)), stat='bin', bins = bins, vjust=-0.5, size=3, color='black') +
    geom_vline(data = plot_data_edge, aes(xintercept = dist),
               color = "red",
               linetype = "dashed", linewidth = 0.75) +
    facet_wrap(~grp, ncol = 1,scales = "free") +
    scale_fill_manual(values = c("#99CCFF","#FF9999")) +
    labs(
      title = "Distance to Nearest Niche Center",
      x = "Distance", y = "Count",
      fill = NULL, color = NULL
    ) +
    theme_test() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_rect(fill="grey98"))
  print(p1)

  # >>> Final
  clog_end()
  return(invisible(NULL))
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NicheDetect_Spot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Spot-level niche detection
#'
#' Performs niche detection at the individual spot level, treating each positive
#' spot as its own ROI. Calculates distances from all spots to the nearest
#' positive spot.
#'
#' @param STID_obj An STID object containing spatial transcriptomics data
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param pos_colnm Character, column name containing positive spot labels
#' @param neg_value Character, value indicating negative spots (default: "neg")
#' @param meta_key Character, metadata key containing positive spot information
#' @param description Character, description of the analysis (default: NULL)
#' @param grp_nm Character, group name for output organization (default: NULL, uses timestamp)
#' @param dir_nm Character, directory name for output (default: "M2_NicheDetect_Spot")
#'
#' @return Returns the modified STID object with added metadata containing
#'         ROI labels and distances to nearest positive spots
#'
#' @import Seurat
#' @import ggdensity
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform spot-level niche detection
#' STID_obj <- NicheDetect_Spot(
#'   STID_obj = ist_object,
#'   meta_key = "M1_SpotDetect_Gene_20240101",
#'   pos_colnm = "Label_geneA"
#' )
#' }
NicheDetect_Spot <- function(STID_obj = NULL, loop_id = "LoopAllSamp",
                             pos_colnm = NULL, neg_value = "neg",
                             meta_key = NULL, description = NULL,
                             grp_nm = NULL, dir_nm = "M2_NicheDetect_Spot"
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
  .check_null_args(meta_key, pos_colnm)
  loop_single <- .check_loop_single(STID_obj = STID_obj,loop_id = loop_id, mode = 1)
  all_single <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "samp_id")[[1]]
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  meta_data <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]]
  .check_column_exist(meta_data, pos_colnm)
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  x_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "x_colnm")[[1]]
  y_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "y_colnm")[[1]]
  data_format <- GetInfo(STID_obj, info_key = "data_info",sub_key = "data_format")[[1]]
  data_platform <- GetInfo(STID_obj, info_key = "data_info",sub_key = "data_platform")[[1]]
  binsize <- GetInfo(STID_obj, info_key = "data_info",sub_key = "binsize")[[1]]
  meta2STID_list <- list()
  n <- 1
  for(i in seq_along(all_single)){
    i_single <- all_single[i]
    i_meta_data <- meta_data %>% filter(!!sym(samp_colnm) == i_single)
    i_meta_data <- i_meta_data[, c(samp_colnm, x_colnm, y_colnm, pos_colnm)]
    if(!i_single %in% loop_single){
      meta2STID_list[[i_single]] <- i_meta_data
      next
    }
    clog_loop(paste0("Processing samp_id: ", i_single, " (", n, "/", length(loop_single), ")"))
    n <- n + 1

    #> pos_meta_data
    pos_meta_data <- i_meta_data %>%
      filter(.,!!sym(pos_colnm) != neg_value)
    len_pos_value <- nrow(pos_meta_data)
    if(len_pos_value == 0){
      clog_warn(paste0("No positive cells found in samp_id: ", i_single, ", skiping..."))
      meta2STID_list[[i_single]] <- i_meta_data
      next
    }
    i_meta_data <- i_meta_data %>%
      mutate(ROI_label = ifelse(!!sym(pos_colnm) != neg_value, rownames(.), "neg"),
             ROI_edge = ifelse(ROI_label == "neg", FALSE, TRUE),
             ROI_center = ROI_edge)
    i_meta_data <- .CalDist2Center(meta_data = i_meta_data,
                                   ROI_label_colnm = "ROI_label",
                                   ROI_center_colnm = "ROI_center",
                                   CalMode = "Spot")
    meta2STID_list[[i_single]] <- i_meta_data
  }

  #> AddMetaData
  meta2STID <- bind_rows(meta2STID_list) %>%
    as.data.frame()
  meta2STID <- meta2STID[rownames(meta_data), , drop = FALSE] # keep order
  new_meta_key <- paste0(dir_nm,"_",grp_nm)
  STID_obj <- AddMetaData(STID_obj = STID_obj,
                         meta_key = new_meta_key,
                         add_data = meta2STID,
                         dir_nm = dir_nm,
                         grp_nm = grp_nm,
                         asso_key = meta_key,
                         description = description)

  # >>> Final
  .save_function_params("NicheDetect_Spot", envir = environment(), file = paste0(output_dir,"Log_function_params_(NicheDetect_Lasso).log") )
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(NicheDetect_Spot).log")) %>% invisible()
  return(STID_obj)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NicheExpand
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Expand niche regions based on distance threshold
#'
#' Expands existing niche regions by including neighboring spots within a
#' specified distance threshold. Useful for defining broader influence zones
#' around detected niches.
#'
#' @param STID_obj An STID object containing niche detection results
#' @param meta_key Character, metadata key containing niche information
#' @param loop_id Character, sample grouping identifier (default: "LoopAllSamp")
#' @param pos_colnm Character, column name containing positive spot labels
#' @param neg_value Character, value indicating negative spots (default: "neg")
#' @param center_colnm Character, column name indicating niche center spots
#'        (default: NULL)
#' @param expand_dist Numeric, distance threshold for expansion (default: 1)
#' @param description Character, description of the analysis
#' @param grp_nm Character, group name for output organization (default: NULL, uses timestamp)
#' @param dir_nm Character, directory name for output (default: "M2_NicheExpand")
#'
#' @return Returns the modified STID object with expanded ROI labels and
#'         updated distance calculations
#'
#' @import Seurat
#' @import FNN
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Expand niches by including spots within distance 5
#' STID_obj <- NicheExpand(
#'   STID_obj = ist_object,
#'   meta_key = "M2_NicheDetect_STS_20240101",
#'   pos_colnm = "ROI_label",
#'   expand_dist = 5
#' )
#' }
NicheExpand <- function(STID_obj = NULL, meta_key = NULL, loop_id = "LoopAllSamp",
                        pos_colnm = NULL,neg_value = "neg", center_colnm = NULL, expand_dist = 1,
                        description = NULL, grp_nm = NULL, dir_nm = "M2_NicheExpand"
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
  .check_null_args(meta_key, pos_colnm)
  loop_single <- .check_loop_single(STID_obj = STID_obj,loop_id = loop_id, mode = 1)
  all_single <- GetInfo(STID_obj, info_key = "samp_info",sub_key = "samp_id")[[1]]
  # >>> End check

  # >>> dir
  dir_list <- .create_directory(grp_nm,dir_nm)
  output_dir <- dir_list$output_dir
  photo_dir <- dir_list$photo_dir
  grp_nm <- dir_list$grp_nm

  # >>> Start main pipeline
  meta_data <- GetMetaData(STID_obj = STID_obj, meta_key = meta_key)[[1]]
  .check_column_exist(meta_data, pos_colnm)
  samp_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "samp_colnm")[[1]]
  x_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "x_colnm")[[1]]
  y_colnm <- GetInfo(STID_obj, info_key = "data_info",sub_key = "y_colnm")[[1]]
  data_format <- GetInfo(STID_obj, info_key = "data_info",sub_key = "data_format")[[1]]
  if(data_format  == "square_grid"){
    k <- 8; edge_thres <- 1
  }else if(data_format == "hex_grid" & data_platform == "Visium"){
    k <- 6; edge_thres <- 1
  }
  meta2STID_list <- list()
  n <- 1
  for(i in seq_along(all_single)){
    i_single <- all_single[i]
    i_meta_data <- meta_data %>% filter(!!sym(samp_colnm) == i_single)
    i_meta_data <- i_meta_data[,c(samp_colnm,x_colnm,y_colnm,center_colnm,pos_colnm)]
    if(!i_single %in% loop_single){
      meta2STID_list[[i_single]] <- i_meta_data
      next
    }
    clog_loop(paste0("Processing samp_id: ", i_single, " (", n, "/", length(loop_single), ")"))
    n <- n + 1


    #>
    pos_meta_data <- i_meta_data %>%
      filter(.,!!sym(pos_colnm) != neg_value)
    len_pos_value <- nrow(pos_meta_data)
    if(len_pos_value == 0){
      clog_warn(paste0("No positive cells found in samp_id: ", i_single, ", skiping..."))
      meta2STID_list[[i_single]] <- i_meta_data
      next
    }
    clog_normal(paste0("Calculating distances to positive ROIs, it may take a while..."))
    coord_all <- i_meta_data[,c(x_colnm,y_colnm)]
    coord_pos <- pos_meta_data[,c(x_colnm,y_colnm)]
    result_knnx <- get.knnx(data = coord_pos, query = coord_all, k = 1)
    min_dist_idx <- result_knnx$nn.index[,1]
    min_dist_values <- result_knnx$nn.dist[,1]
    i_meta_data <- i_meta_data %>%
      mutate(
        minspot = rownames(coord_pos)[min_dist_idx],
        mindist = min_dist_values,
        posROI = ifelse(mindist <= expand_dist, "pos","neg"),
        ROI_label = if_else(posROI == "pos",.[minspot,pos_colnm],!!sym(pos_colnm))) %>% # replace ROI_label
      dplyr::select(-minspot)
    if(!is.null(center_colnm)){
      i_meta_data <- .CalROIedge(meta_data = i_meta_data,
                                 ROI_label_colnm = "ROI_label",
                                 ROI_edge_colnm = "ROI_edge",
                                 k = k,edge_thres = edge_thres)
      i_meta_data <- .CalDist2Center(meta_data = i_meta_data,
                                     ROI_label_colnm = "ROI_label",
                                     ROI_center_colnm = "ROI_center",
                                     CalMode = "ROI")
    }
    meta2STID_list[[i_single]] <- i_meta_data
  }

  #> AddMetaData
  meta2STID <- bind_rows(meta2STID_list) %>%
    as.data.frame()
  meta2STID <- meta2STID[rownames(meta_data), , drop = FALSE] # keep order
  new_meta_key <- paste0(dir_nm,"_",grp_nm)
  STID_obj <- AddMetaData(STID_obj = STID_obj,
                         meta_key = new_meta_key,
                         add_data = meta2STID,
                         dir_nm = dir_nm,
                         grp_nm = grp_nm,
                         asso_key = meta_key,
                         description = description)

  # >>> Final
  .save_function_params("NicheExpand", envir = environment(), file = paste0(output_dir,"Log_function_params_(NicheDetect_Lasso).log") )
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(NicheExpand).log")) %>% invisible()
  return(STID_obj)
}


