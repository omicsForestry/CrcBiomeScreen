# CrcBiomeScreen Vignette
# ------------------------------------------------------------------------------
rm(list = ls())
library(devtools)
devtools::install_github("iChronostasis/CrcBiomeScreen",force = TRUE)
setwd("/home/CRCscreening/CRCscreening-Workflow/")
# ------------------------------------------------------------------------------
# Start the CrcBiomeScreening workflow
# Load required libraries
library(CrcBiomeScreen)
library(ggplot2)
library(dplyr)
# ------------------------------------------------------------------------------
## Get the toy data from curatedMetagenomicData
library(curatedMetagenomicData)
toydata <- curatedMetagenomicData(
            "ThomasAM_2018a.relative_abundance"
            , dryrun = FALSE, rownames = "short")

## Create the CrcBiomeScreenObject
CrcBiomeScreenObject <- CreateCrcBiomeScreenObject(RelativeAbundance = toydata[[1]]@assays@data@listData$relative_abundance,
                                                  TaxaData = toydata[[1]]@rowLinks$nodeLab,
                                                  SampleData = toydata[[1]]@colData)
## Split the taxa into multiple levels
CrcBiomeScreenObject <- SplitTaxas(CrcBiomeScreenObject)

## Keep only the genus level data
CrcBiomeScreenObject <- KeepGenusLevel(CrcBiomeScreenObject)

## Normalize the data using GMPR/TSS
CrcBiomeScreenObject <- NormalizeData(CrcBiomeScreenObject, method = "GMPR")

CrcBiomeScreenObject <- NormalizeData(CrcBiomeScreenObject, method = "TSS")

# ------------------------------------------------------------------------------
# Prepare the validation data

## Load the validation data from curatedMetagenomicData
ValidationData_curated <- curatedMetagenomicData(
            paste0("ZellerG_2014",".","relative_abundance")
            , dryrun = FALSE, rownames = "short")

## Create the CrcBiomeScreenObject for validation data
ValidationData <- CreateCrcBiomeScreenObject(RelativeAbundance = ValidationData_curated[[1]]@assays@data@listData$relative_abundance,
                                                  TaxaData = ValidationData_curated[[1]]@rowLinks$nodeLab,
                                                  SampleData = ValidationData_curated[[1]]@colData)

## Split the taxa into multiple levels
ValidationData <- SplitTaxas(ValidationData)

## Keep only the genus level data
ValidationData <- KeepGenusLevel(ValidationData)

## Normalize the validation data using GMPR/TSS
ValidationData <- NormalizeData(ValidationData, method = "GMPR")

## Keep only the taxa that are present in both training and validation data
ValidationData$NormalizedData <- ValidationData$NormalizedData[, colnames(ValidationData$NormalizedData) %in%
                                            colnames(CrcBiomeScreenObject$NormalizedData)]

CrcBiomeScreenObject$NormalizedData <- CrcBiomeScreenObject$NormalizedData[, colnames(CrcBiomeScreenObject$NormalizedData) %in%
                                            colnames(ValidationData$NormalizedData)]

## Check the labels of the study condition in sample data
table(CrcBiomeScreenObject$SampleData$study_condition)

# ------------------------------------------------------------------------------
# Modeling and Evaluation
CrcBiomeScreenObject <- FilterDataSet(CrcBiomeScreenObject,
                                         label = c("CRC","control"),
                                         condition_col = "study_condition")
table(CrcBiomeScreenObject$SampleData$study_condition)
## Split the data into training and test sets by using the labels and set the partition
CrcBiomeScreenObject <- SplitDataSet(CrcBiomeScreenObject, label = c("control","CRC"), partition = 0.7)

## Optinal: quality control by cmdscale
CrcBiomeScreenObject <- qcByCmdscale(CrcBiomeScreenObject,
                                    TaskName = "Normalize_ToyData_filtered_qc",
                                    normalize_method = "GMPR")

## Example: check balance in your classification labels
checkClassBalance(CrcBiomeScreenObject$ModelData$TrainLabel)


## Train the models using Random Forest and XGBoost
CrcBiomeScreenObject <- TrainModels(CrcBiomeScreenObject,
                                    model_type = "RF",
                                    TaskName = "ToyData_RF",
                                    ClassBalance = TRUE,
                                    TrueLabel = "CRC",
                                    num_cores = 10)

CrcBiomeScreenObject <- TrainModels(CrcBiomeScreenObject,
                                    model_type = "XGBoost",
                                    TaskName = "ToyData_XGBoost",
                                    ClassBalance = TRUE,
                                    TrueLabel = "CRC",
                                    num_cores = 10)

## Evaluate the models on the test set
CrcBiomeScreenObject <- EvaluateModel(CrcBiomeScreenObject,
                                       model_type = "RF",
                                       TaskName = "ToyData_RF_Test",
                                       TrueLabel = "CRC",
                                       PlotAUC = TRUE)

CrcBiomeScreenObject <- EvaluateModel(CrcBiomeScreenObject,
                                       model_type = "XGBoost",
                                       TaskName = "ToyData_XGBoost_Test",
                                       TrueLabel = "CRC",
                                       PlotAUC = TRUE)

# ------------------------------------------------------------------------------
# Validate the models on the validation data

## Filter the validation data to keep only the samples with labels "CRC" and "control"
ValidationData_filtered <- FilterDataSet(ValidationData,
                                         label = c("CRC","control"),
                                         condition_col = "study_condition")

## Quality control the validation data by cmdscale
ValidationData_filtered_qc <- qcByCmdscale(ValidationData_filtered,
                                           TaskName = "Normalize_ValidationData_filtered_qc",
                                           normalize_method = "GMPR")

## Validate the models on the validation data using Random Forest and XGBoost
CrcBiomeScreenObject <- ValidateModelOnData(CrcBiomeScreenObject,
                                            model_type = "RF",
                                            ValidationData = ValidationData_filtered_qc,
                                            TaskName = "ValidationData_RF_Validation",
                                            TrueLabel = "CRC",
                                            PlotAUC = TRUE)

CrcBiomeScreenObject <- ValidateModelOnData(CrcBiomeScreenObject,
                                            model_type = "XGBoost",
                                            ValidationData = ValidationData_filtered_qc,
                                            TaskName = "ValidationData_XGBoost_Validation",
                                            TrueLabel = "CRC",
                                            PlotAUC = TRUE)
# ------------------------------------------------------------------------------
# Run the screening workflow by using one function
CrcBiomeScreenObject <- RunScreening(CrcBiomeScreenObject,
                                    normalize_method = "GMPR",
                                    model = "RF",
                                    partition = 0.7,
                                    split.requirement =
                                    c(label = c("control","CRC"),
                                      condition_col = "study_condition"),
                                    ClassBalance = TRUE,
                                    n_cv = 10,
                                    TaskName = "GMPR_toydata",
                                    num_cores = 10,
                                    ValidationData = ValidationData_filtered_qc,
                                    TrueLabel = "CRC")







