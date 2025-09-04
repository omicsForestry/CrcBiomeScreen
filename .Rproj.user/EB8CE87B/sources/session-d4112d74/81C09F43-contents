#' Load all functions for testing
#' 
#' This function sources all R files in the current directory
#' for testing purposes. This is not meant for the final package.
#'
load_functions <- function() {
  # Get all R files in the current directory (except main.R itself)
  r_files <- list.files(path = "R", pattern = "\\.R$", full.names = TRUE)
  r_files <- r_files[!grepl("main\\.R$", r_files)]
  
  # Load all files
  for(file in r_files) {
    message("Loading: ", basename(file))
    source(file, local = parent.frame())
  }
  
  message("All functions loaded successfully!")
}
load_functions()
