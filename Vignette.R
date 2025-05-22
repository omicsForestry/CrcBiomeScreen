rm(list = ls())
setwd("/home/CRCscreening/CRCscreening-Workflow/")
setwd("/mnt/scratch/ngzh5554/CRCscreening-Workflow")
# Load required libraries
source("R/Environment.R")
source("R/main.R")

# Load the data_relative abundance
# toydata <- readRDS("toy_data.rds")
toydata <- curatedMetagenomicData(
            "ZellerG_2014.relative_abundance"
            , dryrun = FALSE, rownames = "short")  

CrcBiomeScreenObject <- CreateCrcBiomeScreenObject(RelativeAbundance = toydata[[1]]@assays@data@listData$relative_abundance,
                                                  TaxaData = toydata[[1]]@rowLinks$nodeLab,
                                                  SampleData = toydata[[1]]@colData)
CrcBiomeScreenObject <- SplitTaxas(CrcBiomeScreenObject)
CrcBiomeScreenObject <- KeepGenusLevel(CrcBiomeScreenObject)

CrcBiomeScreenObject <- NormalizeData(CrcBiomeScreenObject, method = "GMPR",TaskName = "Normalize_ToyData")
CrcBiomeScreenObject <- SplitDataSet(CrcBiomeScreenObject, label = c("control","CRC"), partition = 0.7)
CrcBiomeScreenObject <- TrainModels(CrcBiomeScreenObject, 
                                    model_type = "RF",
                                    TaskName = "ToyData_RF", 
                                    TrueLabel = "CRC",
                                    num_cores = 10)

CrcBiomeScreenObject <- TrainModels(CrcBiomeScreenObject, 
                                    model_type = "XGBoost",
                                    TaskName = "ToyData_XGBoost", 
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
                                    split.requirement = c(label = c("control","CRC"), partition = 0.7), 
                                    TaskName = "GMPR_NHSBCSP", num_cores = 10)



ValidationData_curated <- curatedMetagenomicData(
            paste0(available_studies[11],".","relative_abundance")
            , dryrun = FALSE, rownames = "short")  
# saveRDS(ValidationData_curated, "ValidationData.rds")
ValidationData <- CreateCrcBiomeScreenObject(RelativeAbundance = ValidationData_curated[[1]]@assays@data@listData$relative_abundance,
                                                  TaxaData = ValidationData_curated[[1]]@rowLinks$nodeLab,
                                                  SampleData = ValidationData_curated[[1]]@colData)
ValidationData <- SplitTaxas(ValidationData)
ValidationData <- KeepGenusLevel(ValidationData)
ValidationData <- NormalizeData(ValidationData, method = "GMPR",TaskName = "Normalize_ValidationData")
ValidationData_filtered <- FilterDataSet(ValidationData, 
                                 label = c("CRC","control"),
                                 condition_col = "study_condition")

CrcBiomeScreenObject <- PredictValidation(CrcBiomeScreenObject, 
                                       model_type = "RF",
                                       ValidationData = ValidationData_filtered,
                                       TaskName = "ValidationData_RF_Validation", 
                                       TrueLabel = "CRC",
                                       PlotAUC = TRUE)
                                       
CrcBiomeScreenObject <- PredictValidation(CrcBiomeScreenObject, 
                                       model_type = "XGBoost",
                                       ValidationData = ValidationData,
                                       TaskName = "ValidationData_XGBoost_Validation", 
                                       TrueLabel = "CRC",
                                       PlotAUC = TRUE)





