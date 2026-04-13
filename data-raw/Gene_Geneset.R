## code to prepare `Gene_Geneset` dataset goes here

library(vroom)
read_txt_as_df_vroom <- function(path, return_filenames_only = FALSE, ...) {
  entries <- list.files(
    path,
    full.names = TRUE,
    include.dirs = TRUE,
    recursive = FALSE,
    no.. = TRUE  # 排除 . 和 ..
  )

  if (length(entries) == 0) return(list())

  names_entries <- basename(entries)
  result <- list()

  for (i in seq_along(entries)) {
    entry <- entries[i]
    name <- names_entries[i]

    if (dir.exists(entry)) {
      result[[name]] <- read_txt_as_df_vroom(
        entry,
        return_filenames_only = return_filenames_only,
        ...
      )
    } else if (grepl("\\.txt$", entry, ignore.case = TRUE)) {
      clean_name <- tools::file_path_sans_ext(name)

      if (return_filenames_only) {
        result[[clean_name]] <- name
      } else {
        df <- tryCatch({
          vroom(entry, ...)
        }, error = function(e) {
          warning("Failed to read '", entry, "': ", conditionMessage(e))
          NULL
        })
        result[[clean_name]] <- df
      }
    }
  }

  return(result)
}

Gene_Geneset <- read_txt_as_df_vroom(
  path = "./inputdata/Gene_Geneset/",
  delim = "\t"
)
saveRDS(Gene_Geneset,file = "./rds/Gene_Geneset.rds")

file_structure <- read_txt_as_df_vroom(
  path = "./inputdata/Gene_Geneset/",
  delim = "\t",
  return_filenames_only = TRUE
)
str(file_structure)









usethis::use_data(Gene_Geneset, overwrite = TRUE)
