rm(list = ls())
setwd("/home/CRCscreening/CRCscreening-Workflow/")
setwd("/mnt/scratch/ngzh5554/CRCscreening-Workflow")
# Load required libraries
# source("R/Environment.R")
# source("R/main.R")

# Load the data_relative abundance
# toydata <- readRDS("toy_data.rds")
# ThomasaAM_a
# toydata <- curatedMetagenomicData(
#             "ZellerG_2014.relative_abundance"
#             , dryrun = FALSE, rownames = "short")

toydata <- curatedMetagenomicData(
            "ThomasAM_2018a.relative_abundance"
            , dryrun = FALSE, rownames = "short")

CrcBiomeScreenObject <- CreateCrcBiomeScreenObject(RelativeAbundance = toydata[[1]]@assays@data@listData$relative_abundance,
                                                  TaxaData = toydata[[1]]@rowLinks$nodeLab,
                                                  SampleData = toydata[[1]]@colData)
CrcBiomeScreenObject <- SplitTaxas(CrcBiomeScreenObject)
CrcBiomeScreenObject <- KeepGenusLevel(CrcBiomeScreenObject)

CrcBiomeScreenObject <- NormalizeData(CrcBiomeScreenObject, method = "GMPR")
# CrcBiomeScreenObject <- qcByCmdscale(CrcBiomeScreenObject,
#                                               TaskName = "ToyData",
#                                               normalize_method = "GMPR")

# ------------------------------------------------------------------------------
ValidationData_curated <- curatedMetagenomicData(
            paste0("ZellerG_2014",".","relative_abundance")
            , dryrun = FALSE, rownames = "short")

ValidationData <- CreateCrcBiomeScreenObject(RelativeAbundance = ValidationData_curated[[1]]@assays@data@listData$relative_abundance,
                                                  TaxaData = ValidationData_curated[[1]]@rowLinks$nodeLab,
                                                  SampleData = ValidationData_curated[[1]]@colData)
ValidationData <- SplitTaxas(ValidationData)
ValidationData <- KeepGenusLevel(ValidationData)
ValidationData <- NormalizeData(ValidationData, method = "GMPR")

ValidationData$NormalizedData <- ValidationData$NormalizedData[, colnames(ValidationData$NormalizedData) %in%
                                            colnames(CrcBiomeScreenObject$NormalizedData)]

CrcBiomeScreenObject$NormalizedData <- CrcBiomeScreenObject$NormalizedData[, colnames(CrcBiomeScreenObject$NormalizedData) %in%
                                            colnames(ValidationData$NormalizedData)]

table(CrcBiomeScreenObject$SampleData$study_condition)

# ------------------------------------------------------------------------------
CrcBiomeScreenObject <- SplitDataSet(CrcBiomeScreenObject, label = c("control","CRC"), partition = 0.7)

# Optinal: quality control by cmdscale
CrcBiomeScreenObject <- qcByCmdscale(CrcBiomeScreenObject,
                                    TaskName = "Normalize_ToyData_filtered_qc",
                                    normalize_method = "GMPR")
# Example: check balance in your classification labels
checkClassBalance(CrcBiomeScreenObject$ModelData$TrainLabel)

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

CrcBiomeScreenObject <- RunScreening(CrcBiomeScreenObject,
                                    normalize_method = "GMPR",
                                    model = "RF",
                                    split.requirement =
                                    c(label = c("control","CRC"),
                                      condition_col = "study_condition"),
                                    TaskName = "GMPR_toydata",
                                    num_cores = 10,
                                    ValidationData = ValidationData,
                                    TrueLabel = "Cancer")



# ValidationData_curated <- curatedMetagenomicData(
#             paste0("ZellerG_2014",".","relative_abundance")
#             , dryrun = FALSE, rownames = "short")
# # saveRDS(ValidationData_curated, "ValidationData.rds")
# ValidationData <- CreateCrcBiomeScreenObject(RelativeAbundance = ValidationData_curated[[1]]@assays@data@listData$relative_abundance,
#                                                   TaxaData = ValidationData_curated[[1]]@rowLinks$nodeLab,
#                                                   SampleData = ValidationData_curated[[1]]@colData)
# ValidationData <- SplitTaxas(ValidationData)
# ValidationData <- KeepGenusLevel(ValidationData)
# ValidationData <- NormalizeData(ValidationData, method = "GMPR",TaskName = "Normalize_ValidationData")

ValidationData_filtered <- FilterDataSet(ValidationData,
                                 label = c("CRC","control"),
                                 condition_col = "study_condition")

ValidationData_filtered_qc <- qcByCmdscale(ValidationData_filtered,
                                              TaskName = "Normalize_ValidationData_filtered_qc",
                                              normalize_method = "GMPR")

CrcBiomeScreenObject <- ValidateModelOnData(CrcBiomeScreenObject,
                                       model_type = "RF",
                                       ValidationData = ValidationData_filtered_qc,
                                       TaskName = "ValidationData_RF_Validation",
                                       TrueLabel = "CRC",
                                       PlotAUC = TRUE)

CrcBiomeScreenObject <- ValidateModelOnData(CrcBiomeScreenObject,
                                       model_type = "XGBoost",
                                       ValidationData = ValidationData,
                                       TaskName = "ValidationData_XGBoost_Validation",
                                       TrueLabel = "CRC",
                                       PlotAUC = TRUE)





