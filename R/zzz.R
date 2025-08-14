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
#
# .onLoad <- function(libname, pkgname) {
#   deps <- c(
#     "caret", "curatedMetagenomicData", "ggplot2", "ggplotify", "ggpubr",
#     "ggrepel", "ggtree", "glmnet", "GUniFrac", "Matrix", "microbiome",
#     "microbiomeMarker", "phyloseq", "ranger", "vegan", "xgboost", "dplyr",
#     "tibble", "tidyr", "foreach", "pROC", "progress", "doParallel",
#     "doFuture", "future", "future.apply", "progressr", "ape"
#   )
#   missing <- deps[!deps %in% installed.packages()[, "Package"]]
#   if (length(missing) > 0) {
#     message("Installing missing dependencies: ", paste(missing, collapse = ", "))
#     install.packages(missing, dependencies = TRUE)
#   }
# }

.onLoad <- function(libname, pkgname) {
  pkgs <- list(
    Matrix = "1.6-5",
    MASS = "7.3-60.0.1",
    ggplot2 = "3.5.2",
    mgcv = "1.9-1",
    rstatix = "0.7.2"
  )
  for (p in names(pkgs)) {
    if (!requireNamespace(p, quietly = TRUE)) {
      remotes::install_version(p, version = pkgs[[p]], dependencies = TRUE)
    }
  }
}

