

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# .STID_options
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get default parallelization plan based on OS

#' Returns the default parallelization plan to use based on the operating system.
# For Windows, it returns "multisession", and for other OS, it returns "multicore".
#' @return Character, "multisession" for Windows, "multicore" for other OS
#'
#' @keywords internal
#' @noRd
.get_default_plan <- function() {
  if (.Platform$OS.type == "windows") {
    "multisession"
  } else {
    "multicore"
  }
}

#' Get default number of parallel workers
#'
#' Returns the default number of parallel workers to use, which is the minimum of 4 or the number of available cores minus 2.
#' @return Integer, number of parallel workers
#' @keywords internal
#' @noRd
.get_parallel_workers <- function() {
  min(4, parallel::detectCores() - 2)
}

#' Internal environment for STID options
#'
#' @keywords internal
#' @noRd
.STID_options <- new.env(parent = emptyenv())
.STID_options$parallel_workers <- .get_parallel_workers()
.STID_options$default_plan <- .get_default_plan()
.STID_options$future.globals.maxSize <- 1 * 1024^3 # 1GB


#' Set STID options
#'
#' Internal function to set STID options if not already set.
#'
#' @keywords internal
#' @noRd
.set_STID_options <- function() {
  for (i in seq_along(.STID_options)) { # .STID_options
    op <- names(.STID_options)[i]
    opt_list <- as.list(.STID_options)
    if (is.null(getOption(op))) {
      do.call(options, stats::setNames(list(opt_list[[i]]), op))
    }
  }
}
.set_STID_options()


#' Get STID options
#'
#' Returns a list of current STID options.
#'
#' @return A list of STID options
#'
#' @export
#'
#' @examples
#' \dontrun{
#' opts <- get_STID_options()
#' print(opts$parallel_workers)
#' }
get_STID_options <- function() {
  as.list(.STID_options)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# .STID_options and STID_globals
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' STID package global variables
#'
#' Internal global variables and constants used throughout the STID package.
#'
#' @keywords internal
#' @noRd
.STID_globals <- new.env(parent = emptyenv())
.STID_globals$supported_types <- c("scRNA", "stRNA")
.STID_globals$supported_formats <- c("square_grid", "hex_grid")
.STID_globals$supported_platforms <- c("StereoSeq", "Visium","VisiumHD","SlideSeq","unknown")
.STID_globals$supported_hosts <- c("human", "mouse", "unknown")
.STID_globals$supported_pathogens <- c("virus", "bacteria", "parasite", "unknown")
.STID_globals$score_methods <- c("AddModuleScore", "AUCell", "UCell", "MeanExp", "SumExp")
.STID_globals$blur_methods <- c("isoblur", "medianblur")
.STID_globals$plot_methods <- c("merge", "single")
.STID_globals$samp_modes <- c("SS", "MS")


#' Get STID global constants
#'
#' Returns a list of global constants used throughout the package.
#'
#' @return A list containing supported formats, hosts, pathogens, etc.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' globals <- get_STID_globals()
#' print(globals$supported_formats)
#' }
get_STID_globals <- function() {
  as.list(.STID_globals)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package onLoad, onAttach and onUnload
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Package load hook
#'
#' @param libname Library name
#' @param pkgname Package name
#'
#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {

  # Set options if not already set
  .set_STID_options()
}


#' Package attach hook
#'
#' @param libname Library name
#' @param pkgname Package name
#'
#' @keywords internal
#' @noRd
.onAttach <- function(libname, pkgname) {
  # Get package version
  pkgVersion <- tryCatch(
    utils::packageVersion(pkgname),
    error = function(e) "unknown"
  )

  # Startup message
  packageStartupMessage(
    paste0(
      "\n",
      "==================================================\n",
      "  STID (Spatial Transcriptomics toolkit for Infectious Diseases)\n",
      "  Version ", pkgVersion, "\n",
      "==================================================\n",
      "\n",
      "Welcome to STID package for analyzing \n",
      "Spatial Transcriptomics of Infectious Disease.\n",
      "\n",
      "Key functionalities:\n",
      "  - Data conversion (h5ad2rds, rds2h5ad)\n",
      "  - Preprocessig (Seurat_pipeline, anno_SingleR)\n",
      "  - Background correction (CorrectBackgroud)\n",
      "  - Spot detection (SpotDetect_Gene, SpotDetect_Geneset)\n",
      "  - Niche identification (NicheDetect_Lasso, NicheDetect_STS, ...)\n",
      "  - Single-Sample analysis (CalSampComp, CalSampDEGs, ...)\n",
      "  - Multi-Sample analysis (CalSampPathoTrack, CalSampGeneTrend, ...)\n",
      "\n",
      "For more information, use ?STID see the package documentation.\n",
      "Github URL: https://github.com/YulongQin/STID \n",
      "Repo URL: https://yulongqin.github.io/STID \n",
      "==================================================\n"
    )
  )
}


#' Package unload hook
#'
#' Clean up parallel backend and close open connections.
#'
#' @param libpath Library path
#'
#' @keywords internal
#' @noRd
.onUnload <- function(libpath) {
  # Close any open sink connections
  while (sink.number() > 0) {
    sink()
  }

  # Clean up future parallel backend if running
  if (requireNamespace("future", quietly = TRUE)) {
    try(future::plan("sequential"), silent = TRUE)
  }

  invisible()
}



