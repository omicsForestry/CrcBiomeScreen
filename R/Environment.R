# The environment setup script for the project
# This script checks for the required R packages and installs any that are missing.
# It also loads the necessary libraries for the analysis. 
# Define required R packages
packages <- c("stats","caret", "ranger", "progress", "doParallel", "future.apply", 
              "doFuture", "future", "foreach", "progressr", "dplyr", 
              "devtools", "glmnet", "pROC", "GUniFrac", "Matrix","microbiomeMarker",
              "curatedMetagenomicData","tidyr","dplyr","tibble","doParallel","foreach",
              "ggplot2","phyloseq","ggpubr","ggrepel","ggplotify","ggtree",
              "ape","BiocManager","microbiome","vegan","xgboost")

# Function to check and install missing packages
# This function checks if a package is installed and installs it if not.
# It also handles Bioconductor packages separately.
# The function takes the package name as an argument and checks if it is installed.
# If the package is not installed, it installs it using install.packages() or BiocManager::install() for Bioconductor packages.
# The function also includes an option to install dependencies.
options("repos" = c(CRAN="https://cran.ma.imperial.ac.uk/"))
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing missing package:", pkg))
    if (pkg %in% c("microbiomeMarker","GUniFrac")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# Install missing packages
# This line applies the install_if_missing function to each package in the packages vector.
invisible(lapply(packages, install_if_missing))

# Load required libraries
# This line loads each package in the packages vector using the library() function.
# The character.only = TRUE argument allows for the package names to be passed as strings.
# This is useful for dynamically loading packages based on the names in the vector. 
lapply(packages, library, character.only = TRUE)


