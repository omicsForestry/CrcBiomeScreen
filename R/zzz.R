.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to CrcBiomeScreen!")

  # C heck for critical packages
  critical_packages <- c("phyloseq", "microbiome", "dplyr", "ggtree")
  missing_critical <- critical_packages[!sapply(critical_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_critical) > 0) {
    packageStartupMessage("Warning: Missing critical packages: ",
                          paste(missing_critical, collapse = ", "))
    packageStartupMessage("Run CrcBiomeScreen::install_dependencies() to install all required packages.")
  }
}
