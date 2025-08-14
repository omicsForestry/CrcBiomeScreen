# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("Welcome to CrcBiomeScreen!")
#
#   # C heck for critical packages
#   critical_packages <- c("phyloseq", "microbiome", "dplyr", "ggtree")
#   missing_critical <- critical_packages[!sapply(critical_packages, requireNamespace, quietly = TRUE)]
#
#   if (length(missing_critical) > 0) {
#     packageStartupMessage("Warning: Missing critical packages: ",
#                           paste(missing_critical, collapse = ", "))
#     packageStartupMessage("Run CrcBiomeScreen::install_dependencies() to install all required packages.")
#   }
# }

.onLoad <- function(libname, pkgname) {
  deps <- c(
    "caret", "curatedMetagenomicData", "ggplot2", "ggplotify", "ggpubr",
    "ggrepel", "ggtree", "glmnet", "GUniFrac", "Matrix", "microbiome",
    "microbiomeMarker", "phyloseq", "ranger", "vegan", "xgboost", "dplyr",
    "tibble", "tidyr", "foreach", "pROC", "progress", "doParallel",
    "doFuture", "future", "future.apply", "progressr", "ape"
  )
  missing <- deps[!deps %in% installed.packages()[, "Package"]]
  if (length(missing) > 0) {
    message("Installing missing dependencies: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  }
}
