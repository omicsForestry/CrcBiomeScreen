setwd("/home/CRCscreening/CRCscreening-Workflow/")
# Load required libraries
source("Environment.R")
source("CRCscreening-workflow.R")
# Load the data
toy_data <- readRDS("toy_data.rds")
# Preprocess the data
## Normalize the data
toy_data_normalize <- ChooseNormalize(toy_data, method = "GMPR",taxa_col = c(3:150),TaskName = "ToyData",)
# Perform the analysis
toy_data_ModelData <- SplitDataSet(toy_data_normalize,label = "StudyCondition", partition = 0.7)

## RF model
best.params <- RF_Modeling(ModelData = toy_data_ModelData,
                        label = "StudyCondition",
                        k.rf = 10,
                        taxa_col = c(3:150),
                        TaskName = "ToyData_RF",
                        TrueLabel = "Cancer")

TestModel(ModelData = toy_data_ModelData,
                      label = "StudyCondition",
                      taxa_col = c(3:150),
                      best.params = best.params,
                      TaskName = "ToyData_RF_Test",
    ModelData <- list(train,test)
                      TrueLabel = "Cancer")

## XGBoost model
XGBoost_Modeling(ModelData = toy_data_ModelData,
                label = "StudyCondition",
                taxa_col = c(3:150),
                TaskName = "ToyData_XGBoost",
                TrueLabel = "Cancer")








