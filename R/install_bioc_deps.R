install_bioc_deps <- function() {
  # List of Bioconductor packages you need
  bioc_pkgs <- c(
    "phyloseq",
    "microbiomeMarker",
    "microbiome",
    "ggtree",
    "curatedMetagenomicData"
  )

  # Check which are missing
  missing_pkgs <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]

  if (length(missing_pkgs) > 0) {
    message("Installing missing Bioconductor packages: ", paste(missing_pkgs, collapse = ", "))

    # Install BiocManager if needed
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
  } else {
    message("All Bioconductor packages are already installed.")
  }
}
