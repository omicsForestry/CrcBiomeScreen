# CrcBiomeScreen

![Pipeline](https://github.com/user-attachments/assets/5434dac9-5392-4825-884e-06a56d232a1e)

An R package for colorectal cancer screening and microbiome analysis.

## Installation

```r
# Install from GitHub
devtools::install_github("iChronostasis/CrcBiomeScreen", build_vignettes = TRUE)
```

Or use conda environment:
```bash
conda env create -f environment.yml
conda activate CrcBiomeScreen
```

## Usage

```r
library(CrcBiomeScreen)
vignette("CrcBiomeScreen")
```

## Package Structure

- **R/**: Core functions
- **data/**: Example datasets  
- **vignettes/**: Usage examples
- **tests/**: Testing functions

## Roadmap

- CompareModels()
- SelectImportanceFeatures()  
- SHAP values analysis
