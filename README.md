# CrcBiomeScreen

![Pipeline](https://github.com/user-attachments/assets/3bbae590-e68c-4d25-8f35-8b64543ab7a2)

An R package for colorectal cancer screening and microbiome analysis.

## Installation

```r
# Install from GitHub
install.packages(c("remotes"))
remotes::install_github("omicsForestry/CrcBiomeScreen", force = TRUE)
```

Or use conda environment:
```bash
conda env create -f environment.yml
conda activate CrcBiomeScreen
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
