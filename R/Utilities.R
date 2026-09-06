

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# clog
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Logging functions for pipeline execution
#'
#' A comprehensive set of logging functions for tracking pipeline execution
#' with different log levels and message types. Supports timestamped messages
#' with various prefixes for different purposes.
#'
#' @param message Character, the log message to display
#' @param char Character, prefix character(s) for the log message (default: ">>>")
#' @param level Character, log level - one of "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"
#'
#' @return NULL (invisible), messages are printed to console
#'
#' @rdname clog
#'
#' @export
#'
#' @examples
#' \dontrun{
#' clog("This is an info message")
#' clog_start()
#' clog_step("Starting analysis step")
#' clog_loop("Processing sample 1 of 10")
#' clog_normal("Normal operation message")
#' clog_warn("This is a warning")
#' clog_end()
#' }
clog <- function(message = "", char = ">>>",level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] [%s] %s %s\n",
                       timestamp,
                       level,
                       char,
                       message)
  cat(log_entry)
}

#' @rdname clog
#' @export
clog_start  <- function(){
  clog(message = "",char = "[START]", level = "INFO")
  cat("\n")
}

#' @rdname clog
#' @export
clog_end  <- function(){
  cat("\n")
  clog(message = "",char = "[END]", level = "INFO")
}

#' @rdname clog
#' @export
clog_check <- function(){
  clog(message = "",char = "[check]", level = "INFO")
}

#' @rdname clog
#' @export
clog_SS <- function(message){
  clog(message = message,char = "[SS]", level = "INFO")
}

#' @rdname clog
#' @export
clog_MS <- function(message){
  clog(message = message,char = "[MS]", level = "INFO")
}

#' @rdname clog
#' @export
clog_step <- function(message){
  cat("\n")
  clog(message,char = "[STEP]", level = "INFO")
}

#' @rdname clog
#' @export
clog_loop <- function(message){
  cat("\n")
  clog(message,char = "[loop]", level = "INFO")
}


#' @rdname clog
#' @export
clog_normal <- function(message){clog(message,char = " ---- ", level = "INFO")}

#' @rdname clog
#' @export
clog_title <- function(message){clog(message,char = " >>>> ", level = "INFO")}

#' @rdname clog
#' @export
clog_warn <- function(message){clog(message,char = " !!!! ", level = "WARN")}

#' @rdname clog
#' @param call. Logical, whether to include call information in error messages (default: TRUE)
#' @export
clog_error <- function(message, call. = TRUE){
  clog(message,char = "!!!", level = "ERROR")
  stop(message,call. =  call.)
  # quit('no', -1)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# other utilities：.check and null processing
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check for NULL arguments
#'
#' Validates that none of the provided arguments are NULL. If any are NULL,
#' throws an error listing the NULL arguments.
#'
#' @param ... Arguments to check for NULL values
#'
#' @return NULL (invisible), throws error if any argument is NULL
#'
#' @keywords internal
#'
#' @noRd
.check_not_any_null <- function(...) {
  dots <- list(...)
  names <- as.character(substitute(list(...)))[-1]
  null_names <- names[vapply(seq_along(dots), function(i) is.null(dots[[i]]), TRUE)]
  if (length(null_names) > 0) {
    msg <- paste("These arguments are NULL. Please provide:", paste(sQuote(null_names), collapse = ", "))
    clog_error(msg)
  }
}

#' Check that at least one argument is not NULL
#'
#' Validates that at least one of the provided arguments is non-NULL.
#' Logs information about NULL and non-NULL arguments.
#'
#' @param ... Arguments to check
#'
#' @return NULL (invisible)
#'
#' @keywords internal
#'
#' @noRd
.check_not_all_null <- function(...) {
  dots <- list(...)
  names <- as.character(substitute(list(...)))[-1]
  null_names <- names[vapply(seq_along(dots), function(i) is.null(dots[[i]]), TRUE)]
  not_null_names <- setdiff(names, null_names)
  len_null <- length(null_names)
  len_not_null <- length(not_null_names)
  if (len_null == length(dots)) {
    clog_error("All arguments are NULL. At least one must be non-NULL.")
  } else {
    if (len_null > 0) {
      clog_normal(paste("The following arguments are NULL:",
                         paste(sQuote(null_names), collapse = ", ")))
    }
    if (len_not_null > 0) {
      clog_normal(paste("The following arguments are non-NULL:",
                        paste(sQuote(not_null_names), collapse = ", ")))
    }
  }
}

#' Check that argument has length 1
#'
#' Validates that the provided argument is of length 1.
#'
#' @param arg Argument to check
#'
#' @return NULL (invisible), throws error if argument length != 1
#'
#' @keywords internal
#'
#' @noRd
.check_one_arg <- function(arg) {
  arg_name <- deparse(substitute(arg))
  if (length(arg) != 1) {
    clog_error(paste("Argument", sQuote(arg_name), "must be of length 1."))
  }
}

#' Check if all arguments are NULL
#'
#' @param ... Arguments to check
#'
#' @return Logical, TRUE if all arguments are NULL, FALSE otherwise
#'
#' @keywords internal
#'
#' @noRd
.all_null <- function(...) {
  if (all(vapply(list(...), is.null, logical(1)))) {
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#' Check if any arguments are NULL
#'
#' @param ... Arguments to check
#'
#' @return Logical, TRUE if any argument is NULL, FALSE otherwise
#'
#' @keywords internal
#'
#' @noRd
.any_null <- function(...) {
  if (any(vapply(list(...), is.null, logical(1)))) {
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#' Convert NULL to NA
#'
#' Replaces NULL values with appropriate NA type based on specified type.
#'
#' @param x Object to check
#' @param type Character, target type - "character", "numeric", "integer", or "logical"
#'
#' @return Original x if not NULL, otherwise NA of specified type
#'
#' @keywords internal
#'
#' @noRd
.null_to_na <- function(x, type = "character") {
  if (is.null(x)) {
    switch(type,
           character = NA_character_,
           numeric   = NA_real_,
           integer   = NA_integer_,
           logical   = NA
    )
  } else x
}

#' Convert NULL to 0
#'
#' Replaces NULL values with 0 or appropriate zero-equivalent.
#'
#' @param x Object to check
#' @param type Character, target type - "numeric", "integer", or "logical"
#'
#' @return Original x if not NULL, otherwise 0 of specified type
#'
#' @keywords internal
#'
#' @noRd
.null_to_0 <- function(x, type = "numeric") {
  if (is.null(x)) {
    switch(type,
           numeric = 0,
           integer = 0L,
           logical = FALSE
    )
  } else x
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# other utilities：check_column/features
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check if columns exist in metadata
#'
#' Validates that specified columns exist in a metadata data frame.
#'
#' @param metadata A data frame to check columns in
#' @param ... Column names to check for existence
#'
#' @return TRUE if all columns exist, throws error otherwise
#'
#' @keywords internal
#'
#' @noRd
.check_column_exist <- function(metadata, ...) {
  if(!inherits(metadata, "data.frame")){
    clog_error("'metadata' must be a data.frame")
  }
  data_name <- as_label(substitute(metadata))
  param_names <- substitute(list(...)) %>% as.character() %>% .[-1]
  ref <- colnames(metadata)
  args <- list(...)
  if(length(args) == 0){
    clog_error(paste0("No columns specified to check in '", data_name, "'. Please provide column names as arguments."))
  }
  if(is.null(args)){
    clog_warn("You have provided NULL for column names.")
  }
  for (i in seq_along(args)) {
    param_name <- param_names[i]
    cols <- args[[i]]
    missing <- setdiff(cols, ref)
    if (length(missing) > 0) {
      clog_error(paste0(
        "Column(s) not found in '", data_name, "' for parameter '", param_name, "': ",
        paste(sQuote(missing), collapse = ", ")
      ))
    }
  }
  invisible(TRUE)
}

#' Validate features exist in STID object
#'
#' Checks which specified features exist in the STID object and returns valid ones.
#' Logs warnings for missing features.
#'
#' @param STID_obj An STID or Seurat object
#' @param features Character vector of feature/gene names to validate
#'
#' @return Character vector of valid features that exist in the object
#'
#' @keywords internal
#'
#' @noRd
.check_features_exist <- function(STID_obj, features) {
  if(!is.null(features)){
    features <- features[nzchar(features) & !is.na(features)]
    valid_features <- features[features %in% rownames(STID_obj)]
    missing_features <- features[!features %in% rownames(STID_obj)]

    if (length(missing_features) > 0) {
      clog_warn(paste("The following features not found:",
                      paste(missing_features, collapse = ", ")))
    }
    if (length(valid_features) == 0) {
      clog_error("No valid features provided in 'feature' argument.")
    }
    return(valid_features)
  }else{
    clog_error("features is NULL")
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# other utilities：directory
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Create output and photo directories for analysis
#'
#' Creates outputdata and photo directories with automatic timestamp naming
#' and interactive handling of existing directories.
#'
#' @param grp_nm Character, group name (uses timestamp if NULL)
#' @param dir_nm Character, directory name (required)
#'
#' @return List with output_dir and photo_dir paths
#'
#' @keywords internal
#' @noRd
.create_directory <- function(grp_nm = NULL,dir_nm = NULL) {
  .check_not_any_null(dir_nm)
  old_opt <- getOption("error")
  on.exit(options(error = old_opt), add = TRUE)

  now_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  if(is.null(grp_nm)){
    grp_nm <- now_time
  }

  clog_normal("Create outputdata and photo directories")
  output_dir <- paste0("./outputdata/",dir_nm,"/",grp_nm,"/")
  photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm,"/")
  clog_normal(paste0("Output directory: ", output_dir))
  clog_normal(paste0("Photo directory: ", photo_dir))
  logic <- .create_directory_interactive(output_dir)
  if(!logic) {options(error = NULL);clog_error("User exited without taking any action")}
  logic <- .create_directory_interactive(photo_dir)
  if(!logic) {options(error = NULL);clog_error("User exited without taking any action")}

  return(list(output_dir = output_dir, photo_dir = photo_dir,
              grp_nm = grp_nm, dir_nm = dir_nm))
}


#' Interactive directory creation
#'
#' Creates a directory with interactive handling of existing directories.
#' Provides options to delete and recreate, keep existing, or exit.
#'
#' @param path Character, directory path to create
#'
#' @return Logical, TRUE if directory is ready for use, FALSE if user chose to exit
#'
#' @keywords internal
#'
#' @noRd
.create_directory_interactive <- function(path) {
  if (dir.exists(path)) {
    choice <- menu(
      choices = c("Delete existing directory and create new one",
                  "Keep existing directory and continue",
                  "Exit without doing anything"),
      title = clog_normal(paste("Directory '", path, "' already exists. Please choose an action:", sep = ""))
    )
    if (choice == 1) {
      unlink(path, recursive = TRUE)  # Remove directory and all contents
      dir.create(path, recursive = TRUE, showWarnings =F)
      clog_normal(paste0("Deleted '", path, "' and created a new directory."))
      return(TRUE)
    } else if (choice == 2) {
      clog_normal(paste0("Kept existing directory '", path, "' and continuing execution."))
      return(TRUE)
    } else if (choice == 3) {
      # clog_error("User chose to exit. No directory created.", call. = F)
      # return("User chose to exit. No directory created.")
      return(FALSE)
    } else {
      clog_error("Invalid choice. No operation performed.", call. = F)
    }
  } else {
    dir.create(path, recursive = TRUE, showWarnings =F)
    clog_normal(paste0("Created new directory: ", path))
    return(TRUE)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# other utilities
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Save function parameters to log file
#'
#' Captures and saves function parameters (both defaults and passed values)
#' to a log file for reproducibility.
#'
#' @param func_name Character, name of the function
#' @param envir Environment to evaluate parameters in (default: parent.frame())
#' @param file Character, output file path
#'
#' @return NULL (invisible), writes to log file
#'
#' @keywords internal
#'
#' @noRd
.save_function_params <- function(func_name, envir = parent.frame(), file = "function_params.log") {
  if (!exists(func_name)) {
    warning("Function ", sQuote(func_name), " not found.")
    return(invisible(NULL))
  }
  func_def <- get(func_name, mode = "function")

  # default values
  func_formals <- formals(func_def)
  all_param_names <- names(func_formals)

  # parameters passed in the call
  call_expr <- match.call(
    definition = func_def,
    call = sys.call(sys.parent()),
    expand.dots = TRUE
  )
  called_params <- as.list(call_expr)
  called_params <- called_params[-1]  # remove function name


  # merge, prioritize called params
  final_values <- func_formals
  final_values <- final_values[names(final_values) != "..."]
  if (length(called_params) > 0) {
    for (nm in names(called_params)) {
      if (nm %in% all_param_names) {
        # final_values[[nm]] <- eval(called_params[[nm]], envir = envir) #  don't eval
        final_values[[nm]] <- called_params[[nm]]
      }
    }
  }
  param_str <- sapply(all_param_names, function(nm) {
    val <- final_values[[nm]]
    return(deparse(val)) # must need deparse
  }, USE.NAMES = TRUE)

  # >>>
  conn <- file(file, open = "a")
  on.exit(close(conn), add = TRUE)
  writeLines(c(
    "=== Function Parameters ===",
    paste0(names(param_str), ": ", param_str),
    "",
    "=== Timestamp ===",
    paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    ""
  ), con = conn)

  return(invisible(NULL))
}

#' Remove all sink connections
#'
#' Closes all open sink connections to clean up file redirections.
#'
#' @return NULL (invisible)
#'
#' @keywords internal
#'
#' @noRd
.remove_sinks <- function() {
  while (sink.number() > 0) {
    sink()
  }
}

#' Check if numeric values are integers
#'
#' Determines which numeric values are effectively integers within
#' machine precision.
#'
#' @param x Numeric vector
#'
#' @return Logical vector indicating which values are integers
#'
#' @keywords internal
#'
#' @noRd
.is_integer <- function(x) {
  if (!is.numeric(x)) return(rep(FALSE, length(x)))
  is.finite(x) & (abs(x - round(x)) < .Machine$double.eps^0.5)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# theme
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Spatial plot theme
#'
#' A ggplot2 theme optimized for spatial transcriptomics plots with
#' minimal grid lines and clean background.
#'
#' @param base_size Numeric, base font size (default: 11)
#'
#' @return A ggplot2 theme object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ggplot(spatial_data, aes(x, y, color = expression)) +
#'   geom_point() +
#'   theme_spatial()
#' }
theme_spatial <- function(base_size = 11) {
  theme_bw(base_size = base_size) + # %+replace% 是替换，而非添加 +
    theme(
      # legend.position = "none",
      panel.grid = element_blank(),
      # plot.margin = margin(0,0,0,0,'cm'),


      # plot.title = element_text(face = 'bold',hjust = 0.5),
      # axis.title = element_text(face = 'bold',hjust = 0.5),
      # axis.text.x = element_text(face = 'plain',hjust = 0.5,angle = 0),
      # axis.text.y  = element_text(face = 'bold',hjust = 1),
      # legend.title = element_text(face = 'bold',hjust = 0),
      # legend.text = element_text(face = 'bold',hjust = 0),
      # legend.key.size = unit(8, "pt"),
      # legend.margin = margin(0,0,0,0,'cm'),
      # legend.box.spacing = unit(0.2, "cm"),
      # legend.justification = "left"
    )
}

#' Common analysis theme
#'
#' A versatile ggplot2 theme for general analysis plots with options
#' for continuous or discrete x-axis.
#'
#' @param base_size Numeric, base font size (default: 11)
#' @param x_type Character, x-axis type - "continuous" or "discrete"
#'        (default: "continuous")
#'
#' @return A ggplot2 theme object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # For continuous x-axis
#' ggplot(continuous_data, aes(x, y)) +
#'   geom_point() +
#'   theme_common(x_type = "continuous")
#'
#' # For discrete x-axis with rotated labels
#' ggplot(discrete_data, aes(category, value)) +
#'   geom_bar() +
#'   theme_common(x_type = "discrete")
#' }
theme_common <- function(base_size = 11,x_type = "continuous") {
  res_theme <- theme_bw(base_size = base_size) + # %+replace% 是替换，而非添加 +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = 'bold',hjust = 0.5),
      axis.title = element_text(face = 'bold',hjust = 0.5),
      # axis.text.x = element_text(face = 'plain',hjust = 1,angle = 45),
      axis.text.y  = element_text(face = 'plain',hjust = 1),
      legend.title = element_text(face = 'bold',hjust = 0),
      legend.text = element_text(face = 'plain',hjust = 0),
      # legend.key.size = unit(8, "pt"),
      # legend.margin = margin(0,0,0,0,'cm'),
      # legend.box.spacing = unit(0.2, "cm"),
      legend.justification = "left",
      strip.background = element_rect(fill="grey97")
    )
  if(x_type == "discrete"){
    res_theme <- res_theme + theme(
      axis.text.x = element_text(face = 'plain',hjust = 1,angle = 45)
    )
  } else if(x_type == "continuous"){
    res_theme <- res_theme + theme(
      axis.text.x = element_text(face = 'plain',hjust = 0.5,angle = 0)
    )
  }
  return(res_theme)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# color process
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Adjust hex color values
#'
#' Modifies hex colors by adjusting hue, saturation, and value components.
#'
#' @param hex_vec Character vector of hex color codes
#' @param add_h Numeric, adjustment to hue (-1 to 1, default: 0)
#' @param add_s Numeric, adjustment to saturation (-1 to 1, default: 0)
#' @param add_v Numeric, adjustment to value/brightness (-1 to 1, default: 0)
#'
#' @return Character vector of adjusted hex color codes
#'
#' @import grDevices
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' .hex_add("#FF0000", add_s = -0.3, add_v = 0.2)
#' }
.hex_add <- function(
    hex_vec,
    add_h=0, # 色相
    add_s=0, # 饱和度
    add_v=0  # 明度
) {
  res_hex_vec <- c()
  for(i in 1:length(hex_vec)) {
    hex_color <- hex_vec[i]
    res_rgb <- col2rgb(hex_color)
    r <- res_rgb[1]
    g <- res_rgb[2]
    b <- res_rgb[3]

    res_hsv <- rgb2hsv(r,g,b,maxColorValue = 255)
    h <- res_hsv[1]
    s <- res_hsv[2]
    v <- res_hsv[3]
    # print(paste("The h,s,v value of",hex_color,"is",h,s,v,sep=" "))

    if(h+add_h>1) {
      h <- 1
      # print(paste0("The h value of ",hex_color," is greater than 1"))
    }else if(h+add_h<0) {
      h <- 0
      # print(paste0("The h value of ",hex_color," is less than 0"))
    } else {
      h <- h+add_h
    }
    if(s+add_s>1) {
      s <- 1
      # print(paste0("The s value of ",hex_color," is greater than 1"))
    }else if(s+add_s<0) {
      s <- 0
      # print(paste0("The s value of ",hex_color," is less than 0"))
    } else {
      s <- s+add_s
    }
    if(v+add_v>1) {
      v <- 1
      # print(paste0("The v value of ",hex_color," is greater than 1"))
    }else if(v+add_v<0) {
      v <- 0
      # print(paste0("The v value of ",hex_color," is less than 0"))
    } else {
      v <- v+add_v
    }

    res_hex <- hsv(h=h, s=s, v=v)
    res_hex_vec[i] <- res_hex
  }
  return(res_hex_vec)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# color
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Predefined color palettes
#'
#' A collection of color palettes for consistent visualization across the package.
#' Includes palettes for discrete and continuous scales, with options for
#' different backgrounds.
#'
#' @format
#' \describe{
#'   \item{COLOR_DIS_CON}{List with discrete (2 colors) and continuous (5 colors) palettes}
#'   \item{COLOR_LIST}{Named list of various color palettes:
#'     \itemize{
#'       \item{PALETTE_WHITE_BG}{Extended palette for white background}
#'       \item{PALETTE_BLACK_BG}{Extended palette for black background}
#'       \item{PALETTE_7_VECTOR}{7-color palette with variations}
#'       \item{PALETTE_7_CLASSIC}{Classic 7-color palette}
#'       \item{PALETTE_7_LIGH}{Light 7-color palette}
#'       \item{PALETTE_9_CLASSIC}{Classic 9-color palette}
#'       \item{PALETTE_9_VIRIDIS}{Viridis-based 9-color palette}
#'     }
#'   }
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Use discrete palette
#' ggplot(data, aes(x, y, color = group)) +
#'   geom_point() +
#'   scale_color_manual(values = COLOR_DIS_CON$dis)
#'
#' # Use continuous palette
#' ggplot(data, aes(x, y, color = value)) +
#'   geom_point() +
#'   scale_color_gradientn(colors = COLOR_DIS_CON$con)
#'
#' # Use white background palette
#' ggplot(data, aes(x, y, color = category)) +
#'   geom_point() +
#'   scale_color_manual(values = COLOR_LIST$PALETTE_WHITE_BG)
#' }
COLOR_DIS_CON <- list(
  dis = c("grey95","#E34D4A"),
  con = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
)

#' @rdname COLOR_DIS_CON
#' @export
COLOR_LIST <- list(
  PALETTE_WHITE_BG  = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#66C2A5",
    "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3CC", "#8DD3C7", "#BEBADA", "#FB8072",
    "#FDB462", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#FB9A99",
    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFC0CB", "#B15928", "#B3DE69", "#FBB4AE", "#CCC49C",
    "#F768A1", "#AE017E", "#7A0177", "#49006A", "#E31A1C", "#FF7F00", "#FFFF33", "#A6CEE3", "#1F78B4", "#B2DF8A",
    "#33A02C", "#FB9A99", "#FDBF6F"
  ),
  PALETTE_BLACK_BG = c(
    "#FF8C00", "#FFE4E1", "#0000A6", "#1CE6FF", "#F6CAE5", "#8983BF", "#808000", "#FF4A46",
    "#F5F5DC", "#927d85", "#A30059", "#FF90C9", "#8FBC8F", "#004EB0", "#f5cac3", "#BA0900",
    "#A52A2A", "#8FB0FF", "#e14aec", "#63FFAC", "#F0988C", "#e8612C", "#669bbc", "#FFD39B",
    "#FFBBFF", "#A1A9D0", "#84a59d", "#FFFF00", "#71c33a", "#00868B", "#edc092", "#C4A5DE",
    "#faf3dd", "#A98BC6", "#bc6c25", "#2F7FC1", "#faa307", "#F8AC8C", "#9395E7", "#B1CE46",
    "#63E398", "#06d6a0", "#FFD700", "#7FFF00", "#00FFFF", "#FA8072"
  ),
  PALETTE_7_VECTOR = c("#F81B02FF" ,"#FC7715FF" ,"#FFCE33" ,"#AFBF41FF" ,"#50C49FFF" ,"#3B95C4FF" ,"#B560D4FF",
                      "#F8422D" ,"#FC9041" ,"#FFD960" ,"#B3BF62" ,"#72C4AA" ,"#5DA1C4" ,"#BF85D4",
                      "#F86959" ,"#FCAA6D" ,"#FFE38C" ,"#B7BF84" ,"#95C4B5" ,"#80ADC4" ,"#C9AAD4",
                      "#F88A7E" ,"#FCC093" ,"#FFEDB3" ,"#BBBFA1" ,"#B2C4BE" ,"#9DB7C4" ,"#D1CAD4"),
  PALETTE_7_CLASSIC = c("#F81B02FF" ,"#FC7715FF" ,"#FCB11C" ,"#AFBF41FF" ,"#50C49FFF" ,"#3B95C4FF" ,"#B560D4FF") ,
  PALETTE_7_LIGH = c("#F88A7E", "#FCC093", "#FCDB9A", "#BBBFA1", "#B2C4BE", "#9DB7C4", "#D1CAD4"),
  PALETTE_9_CLASSIC = c("#F81A01", "#FB660F", "#FC9518", "#E9B527", "#AEBF40", "#73C289", "#4BACB2", "#678AC8", "#B45FD3"),
  PALETTE_9_VIRIDIS = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF", "#5DC863FF", "#AADC32FF", "#FDE725FF")
)





