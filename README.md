# CrcBiomeScreen

![Pipeline](https://github.com/user-attachments/assets/5434dac9-5392-4825-884e-06a56d232a1e)

## Description
* **Vignette.R**: Toy sample workflow
* **Environment.R**: Environment setup
* **R/**: Core functions
* **dev/**: Development functions
* **Dataset/**: Toy datasets
* **tests/**: Testing functions

## Installation
### Install the R packages repository
```R
library(devtools)
devtools::install_github("iChronostasis/CrcBiomeScreen", force = TRUE)
```

Or you could choose construct the environment by using the conda:
```bash
conda env create -f environment.yml
conda activate CrcBiomeScreen
```
Then download the R package from github.

## Usage
 * Run the vignette to see how to use this workflow.
```R
library(CrcBiomeScreen)
vignette("vignettes/CrcBiomeScreen.Rmd")
```

## Plan
 * In the future, more functions(modules) will put in this workflow.➡️
   * CompareModels()
   * SelectImportanceFeatures()
   * SHAP values...
