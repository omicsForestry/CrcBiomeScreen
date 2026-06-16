# CrcBiomeScreen

![Pipeline](https://github.com/user-attachments/assets/3bbae590-e68c-4d25-8f35-8b64543ab7a2)

An R package for colorectal cancer screening and microbiome analysis.

## Installation

### Option 1 — Install directly from GitHub
You can install **CrcBiomeScreen** directly in R:

```r
install.packages("remotes")
remotes::install_github("omicsForestry/CrcBiomeScreen", force = TRUE)
```

### Option 2 — Use the provided conda environment

This will set up all system and R dependencies.
After activating the environment, you still need to install the package itself from GitHub.

```bash
# Create and activate environment
conda env create -f environment.yml
conda activate CrcBiomeScreen
```

Then, inside R:

```R
install.packages("remotes")
remotes::install_github("omicsForestry/CrcBiomeScreen", force = TRUE)
```

## Usage

```r
library(CrcBiomeScreen)
```

### Run Vignette
```r
# Install BiocStyle for vignette rendering
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BiocStyle")
library(BiocStyle)
remotes::install_github("omicsForestry/CrcBiomeScreen", build_vignettes = TRUE, force = TRUE)
vignette("CrcBiomeScreen")
```

## Package Structure

- **R/**: Core functions
- **data/**: Example datasets  
- **vignettes/**: Usage examples
- **tests/**: Testing functions
- **man/** : Documentation files

## Roadmap

- CompareModels()
- SelectImportanceFeatures()  
- SHAP values analysis

This package is licensed under the [MIT License](LICENSE.md).
