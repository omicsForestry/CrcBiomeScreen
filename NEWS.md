# CrcBiomeScreen 0.99.0 (2025-10-21)

## Overview
**CrcBiomeScreen** is an R package designed to streamline microbiome-based colorectal cancer (CRC) screening workflows.  
It provides standardized functions for preprocessing, taxonomic data handling, machine learning model training,  
and cross-cohort validation — supporting reproducible and interpretable microbiome analysis for biomarker discovery.

This version marks the **first public release** of the package, submitted to **Bioconductor**.

---

## New Features

### Data Processing
- **`SplitTaxas()`** now automatically detects taxonomy separators  
  (supports `|`, `.`, `_`, and `;`) to handle common formats from MetaPhlAn, QIIME, etc.  
  Adds:
  - `OriginalTaxa` column to retain the raw taxonomy string.
  - Special handling for *uncultured* and *unclassified* taxa  
    (e.g. converts `f__Rikenellaceae|g__unclassified` → `Rikenellaceae_unclassified`).

- **`KeepTaxonomicLevel()`** filters data at a user-defined rank (e.g., `Genus`, `Family`, `Species`)  
  and automatically collapses lower-level abundances.  
  Handles multiple nested unclassified levels (e.g., `D_2__Clostridia.D_3__uncultured.D_4__uncultured`).

---

### Machine Learning
- **`TrainModel()`** serves as a unified interface for training Random Forest (RF) and XGBoost classifiers.  
  - Integrates with internal preprocessing and class-weight balancing.  
  - Uses `withr::with_seed()` for local reproducibility instead of `set.seed()`.  
  - Adds support for both *training* and *external validation* workflows.

- **`EvaluateRF()`** and **`EvaluateXGBoost()`** now return standardized performance metrics  
  (AUC, accuracy, recall, F1) and store model outputs in the main object structure.

- **`ValidateModelOnData()`** supports model evaluation across independent datasets.

---

### Data Object Structure
- Introduced the **`CrcBiomeScreenObject`** class to store:
  - Absolute and relative abundance tables  
  - Taxonomic annotations  
  - Model training results  
  - Validation data and predictions  
  - Visual outputs (e.g., ROC, variable importance)

This ensures data provenance and reproducibility across the full workflow.

---

### Utility and Visualization
- Added built-in plotting functions for:
  - Model ROC curves (`PlotAUC`)
---

## Documentation and Testing
- Added **comprehensive vignette**:
  - Demonstrates preprocessing, genus-level filtering, model training, and validation.
  - Includes examples for multiple taxonomy formats and classifier comparisons.

- Implemented **unit tests** under `tests/testthat/` for key components:
  - Taxa splitting
  - Level filtering
  - Model reproducibility

---

## Technical and Compliance Updates
- Adopted **MIT license**.
- Added `BugReports` and `URL` fields in `DESCRIPTION`.
- Removed unnecessary system files (`.Rproj`, `.DS_Store`).
- Replaced all instances of `set.seed()` with `withr::with_seed()` for Bioconductor compliance.
- Reduced hard-coded dependencies and improved optional imports.

---

## Future Plans
- Add support for additional classifiers.
- Integrate feature interpretation.
- Provide reproducible benchmarking across public CRC datasets.

---

**Maintainer:** Li Chengxin (University of Leeds)  
**Date:** 2025-10-21  
**License:** MIT  
**Repository:** [https://github.com/iChronostasis/CrcBiomeScreen](https://github.com/iChronostasis/CrcBiomeScreen)
