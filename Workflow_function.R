#' Set up one obeject to save the results
#' Create a list to store the results
CreateCrcBiomeScreenObject <- function(AbsoluteAbundance = NULL, TaxaData = NULL, SampleData = NULL, RelativeAbundance = NULL) {
  # If AbsoluteAbundance is NULL, check if RelativeAbundance is provided
  if (!is.null(RelativeAbundance) && is.null(AbsoluteAbundance)) {
    if (is.null(SampleData)) {
      stop("SampleData is required to convert RelativeAbundance to AbsoluteAbundance.")
    }
    if (!"number_reads" %in% colnames(SampleData)) {
      stop("SampleData must contain a 'number_reads' column to convert RelativeAbundance to AbsoluteAbundance.")
    }
      AbsoluteAbundance <- RelativeAbundance %>%
      t() %>%
      data.frame() %>%
      mutate(across(seq_len(dim(RelativeAbundance)[2]), ~ (. * SampleData$number_reads / 100))) %>%
      t() %>%
      data.frame()
  }
  
  # Set up the object
  obj <- list(
    AbsoluteAbundance = AbsoluteAbundance,
    TaxaData = TaxaData,
    SampleData = SampleData,
    RelativeAbundance = RelativeAbundance, 
    GenusLevelData = NULL,
    NormalizedData = NULL,
    ModelData = NULL,
    ScreeningResult = NULL,
    ModelResult = NULL,
    PredictResult = NULL,
    Params = list()
  )
  
  # 设置类名
  class(obj) <- "CrcBiomeScreenObject"
  return(obj)
}



SplitTaxas <- function(CrcBiomeScreenObject) {
    CrcBiomeScreenObject$TaxaData <- 
        CrcBiomeScreenObject$TaxaData %>% 
        tibble(variable = CrcBiomeScreenObject$TaxaData) %>%
        separate(variable, into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
                sep = '\\|', fill = "right") %>%
        mutate(across(everything(), ~ ifelse(. == "", NA, sub("^[a-z]__", "", .)))) %>%
        as.data.frame()
    return(CrcBiomeScreenObject)
}


KeepGenusLevel <- function(CrcBiomeScreenObject){
CrcBiomeScreenObject$GenusLevelData <- 
    CrcBiomeScreenObject$AbsoluteAbundance %>%
        as.data.frame() %>%
        mutate(genus = CrcBiomeScreenObject$TaxaData$Genus) %>%
        group_by(genus) %>%
        summarise_all(sum) %>%
        column_to_rownames("genus") %>%
        t() %>%
        as.data.frame()
        return(CrcBiomeScreenObject)
}

#' Normalization - Choose
#'
#' @Description:
#' This function is for normalization
#'
#' @param data The original Datasets(Row for samples, Column for taxa)
#' @param taxa_col Select the needed Column for normalize
#' @param method Choose the normalization method
#' @name TaskName Save the results
#' 
#' @export 
#'
NormalizeData <- function(CrcBiomeScreenObject = NULL,
                            method = NULL,
                            TaskName = NULL){
    Data <- CrcBiomeScreenObject$GenusLevelData
    if(method == "TSS"){
                # Transforming into relative abundance
                Data <-  t(normalize(t(Data), method = "TSS"))
        }else if(method == "GMPR"){
            size.factor <- GMPR(t(Data))
            size.factor[is.na(size.factor)] <- mean(size.factor,na.rm = TRUE)
            Data <- Data / size.factor 
        }

    # Add the normalization method to the Data
    attr(Data, "NormalizationMethod") <- method
    attr(Data, "TaskName") <- TaskName
    attr(Data, "Timestamp") <- Sys.time()
    CrcBiomeScreenObject$NormalizedData <- Data
    # Save the Data
    saveRDS(CrcBiomeScreenObject,paste0("CrcBiomeScreenObject_",TaskName,".rds"))
    return(CrcBiomeScreenObject)
}


#' Training, Test Splits
#'
#' @Description:
#' This function is for spliting the dataset
#' for train and test
#'
#' @param data The all Dataset
#' @param label The label are used for modeling and prediction
#' @param partition The partition for split the datasets
#'
#' @name TaskName Save the running results
#' @export List for the Datasets
#'
SplitDataSet <- function(CrcBiomeScreenObject = NULL,
                         label = NULL,
                         partition = 0.7,
                         condition_col = "study_condition") {
  # Check if the required parameters are provided
  if (is.null(CrcBiomeScreenObject)) stop("CrcBiomeScreenObject cannot be NULL.")
  if (is.null(label)) stop("Label cannot be NULL.")
  if (is.null(partition) || partition <= 0 || partition >= 1) {
    stop("Partition must be a value between 0 and 1.")
  }
  if (!condition_col %in% colnames(CrcBiomeScreenObject$SampleData)) {
    stop(paste("Condition column", condition_col, "not found in SampleData."))
  }
  
  # Select the data based on the label
  sample_condition <- CrcBiomeScreenObject$SampleData[[condition_col]]
  data <- CrcBiomeScreenObject$NormalizedData[sample_condition %in% label, ]
  sample_condition <- sample_condition[sample_condition %in% label]
  
  # Check if the data is empty
  if (nrow(data) == 0) stop("No data found for the specified label.")
  
  # Create training and test set indexes
  set.seed(123)
  trainIndex <- createDataPartition(sample_condition, p = partition, list = FALSE)
  
  # Split the data
  train <- data[trainIndex, ]
  test <- data[-trainIndex, ]
  
  # Add labels to the training and test sets
  train_label <- sample_condition[trainIndex]
  test_label <- sample_condition[-trainIndex]
  train[[condition_col]] <- as.factor(train_label)
  test[[condition_col]] <- as.factor(test_label)
  
  # Create a list to store the training and test sets
  ModelData <- list(
    Training = train,
    Test = test
  )
  
  # Add attributes to the ModelData
  attr(ModelData, "Split Partition") <- partition
  attr(ModelData, "Task Name") <- "SplitDataSet"
  attr(ModelData, "Timestamp") <- Sys.time()
  attr(ModelData, "Training Size") <- nrow(train)
  attr(ModelData, "Test Size") <- nrow(test)
  
  CrcBiomeScreenObject$ModelData <- ModelData
  saveRDS(CrcBiomeScreenObject, paste0("CrcBiomeScreenObject_", "SplitDataSet", ".rds"))
  
  return(CrcBiomeScreenObject)
}



#' Random Forest Model
#' Modeling & Cross Validation
#' @Description:
#' This function is for modeling for 
#' random forest and use the cross validation to
#' get the best parameters and model
#'
#' @param ModelData List for the Datasets
#' @param label The label are used for modeling
#' @param k.rf The parameter for cross validation
#' @param taxa_col Select the needed Column for normalize
#' @name TaskName Save the running results
#' 
#' @export 
#' 
#' @return best.params.rf The best parameters for model
#'
RF_Modeling <- function(ModelData = NULL,
                        label = NULL,
                        k.rf = NULL,
                        taxa_col = NULL,
                        TaskName = NULL,
                        TrueLabel = NULL
                        ){
    set.seed(123)
    folds.rf <- createFolds(get(label,ModelData[[1]]), k = k.rf)
    
    # Using ranger random forest for faster implementation
    grid.rf <- expand.grid(
    mtry         = seq(5,25, by = 5),
    node_size    = seq(3, 15, by = 2),
    sample_size  = c(.55, .632, .70, .80),
    num.trees    = seq(600, 1000, by = 150),
    AUC     = 0
    )                        

    for(i in 1:nrow(grid.rf)) {

    aucs <- c()  # Initialising a vector to store AUCs for each fold

    for(j in 1:k.rf) {
        val.indices <- folds.rf[[j]]  # Indices for the current fold
        val.data <- ModelData[[1]][val.indices, ]  # Validation data for the current fold
        train.fold.data <- ModelData[[1]][-val.indices, ]  # Training data for the current fold
        
        # Class weights in each fold
        class_weights <- table(get(label,train.fold.data))  # Calculate the number of samples in each class
        class_weights <- sum(class_weights) / (length(class_weights) * class_weights)  # Calculate the weight

        colnames(train.fold.data) <- make.names(colnames(train.fold.data)) # nolint
        colnames(val.data) <- make.names(colnames(val.data))
        loc <- which(label == colnames(train.fold.data))

        # Model training with the specified hyperparameters
        model <- ranger(
        formula         = as.formula(paste(label, "~ .")), 
        data            = train.fold.data[,c(taxa_col,loc)], 
        num.trees       = grid.rf$num.trees[i],
        mtry            = grid.rf$mtry[i],
        min.node.size   = grid.rf$node_size[i],
        sample.fraction = grid.rf$sample_size[i],
        class.weights   = class_weights,  
        seed            = 123,
        classification  = TRUE,
        probability     = TRUE,
        verbose         = FALSE
        )
        
        # Prediction on validation data
        predictions <- predict(model, data = val.data[,taxa_col], type = "response")$predictions

        # Computing AUC
        roc.obj <- roc(get(label,val.data), predictions[, TrueLabel])
        auc <- auc(roc.obj)
        
        # AUC on the current fold
        aucs <- c(aucs, auc)
    }

    # Averaging AUC for the current hyperparameter combination
    grid.rf$AUC[i] <- mean(aucs)
    }
    # saveRDS(grid.rf,paste0("grid.rf_",TaskName,".rds"))
    best.params.index.rf <- which.max(grid.rf$AUC)
    best.params.rf <- grid.rf[best.params.index.rf, ]
    saveRDS(best.params.rf,paste0("best.params.rf_",TaskName,".rds"))
    print("Save the result sucessfully!")
    
    return(best.params.rf)
}


#' Test the Model
#'
#' @Description:
#' This function is for testing the model
#' and get the AUC on test dataset
#'
#' @param ModelData List for the Datasets
#' @param label The label are used for modeling
#' @param taxa_col Select the needed Column
#' @param best.params The best parameters for model
#' @name TaskName Save the running results
#' @name TrueLabel The positive class
#' 
#' @export 
#'
TestModel <- function(ModelData = NULL,
                      label = NULL,
                      taxa_col = NULL,
                      best.params = NULL,
                      TaskName = NULL,
                      TrueLabel = NULL){

    set.seed(123)
    # Class weights in each fold
    class_weights <- table(get(label,ModelData[["Training"]]))  # Calculate the number of samples in each class
    class_weights <- sum(class_weights) / (length(class_weights) * class_weights)  # Calculate the weight

    colnames(ModelData[["Training"]]) <- make.names(colnames(ModelData[["Training"]])) # nolint
    colnames(ModelData[["Test"]]) <- make.names(colnames(ModelData[["Test"]]))
    loc <- which(label == colnames(ModelData[["Training"]]))

    # Retraining the model with the best hyperparameters on the entire training dataset
    rf.Model <- ranger(
    formula         = as.formula(paste(label, "~ .")), 
    data            = ModelData[["Training"]][,c(taxa_col,loc)], 
    num.trees       = best.params$num.trees,
    mtry            = best.params$mtry,
    min.node.size   = best.params$node_size,
    sample.fraction = best.params$sample_size,
    class.weights   = class_weights,
    seed            = 123,
    classification  = TRUE,
    probability     = TRUE
    )

    # Evaluate the model on the test dataset
    # Generating predictions (probabilities for Neoplasm or cancer) for positive class
    test.predictions.rf <- predict(rf.Model, data = ModelData[["Test"]][,taxa_col], type = "response")$predictions
    test.pred.prob.rf <- test.predictions.rf[, TrueLabel] 
    # Actual labels
    test.actual.classes.rf <- ModelData[["Test"]][[label]]

    # calculating the ROC Curve
    roc.curve.rf <- roc(test.actual.classes.rf, test.pred.prob.rf, levels = levels(ModelData[["Test"]][[label]]))

    # Plot
    pdf(paste0("roc.curve.rf.",TaskName,".pdf"))
    plot(roc.curve.rf, print.auc = TRUE, print.thres = TRUE)
    dev.off()

    saveRDS(roc.curve.rf,paste0("roc.curve.rf.",TaskName,".rds"))
    print("Save the result sucessfully!")
    # calculating AUC
    auc.value.roc.curve.rf <- auc(roc.curve.rf)
    
    return(print(paste("AUC: ", auc.value.roc.curve.rf)))
}



#' Modeling & Cross Validation
#' XGBoost
#'
#' @Description:
#' This function is for modeling for 
#' XGBoost and use the cross validation to
#' get the best parameters and model
#'
#' @param ModelData List for the Datasets
#' @param label The label are used for modeling
#' @param k.rf The parameter for cross validation
#' @param taxa_col Select the needed Column for normalize
#' @name TaskName Save the running results
#' 
#' @export 
#' 
#' @return best.params.rf The best parameters for model
#'
# XGBoost_Modeling <- function(ModelData = NULL,
#                         label = NULL,
#                         taxa_col = NULL,
#                         TaskName = NULL,
#                         TrueLabel = NULL
#                         ){
#     set.seed(123)
#     # Calculate the number of cores
#     no_cores <- detectCores() - 1
#     no_cores
#     # Register the cores
#     cl <- makePSOCKcluster(no_cores)
#     registerDoParallel(cl)

#     # tuneGrid for xgbTree
#     grid.xgb <- expand.grid(nrounds = c(100, 200, 300),
#                                         max_depth = c(3, 5, 7, 9),
#                                         eta = c(0.01, 0.1, 0.3),
#                                         gamma = 0,
#                                         colsample_bytree = c(0.5, 0.75, 1),
#                                         min_child_weight = 1,
#                                         subsample = c(0.5, 0.75, 1))

#     colnames(ModelData[["Training"]]) <- make.names(colnames(ModelData[["Training"]])) # nolint
#     colnames(ModelData[["Test"]]) <- make.names(colnames(ModelData[["Test"]]))
#     loc <- which(label == colnames(ModelData[["Training"]]))

#     # cross-validation control
#     trControl.xgb <- trainControl(method = 'repeatedcv', number = 5, repeats = 3,
#                                                 summaryFunction = twoClassSummary,
#                                                 classProbs = TRUE,
#                                                 allowParallel = TRUE)

#     # XGBoost Model Building
#     set.seed(123)
#     xgb <- train(as.formula(paste(label, "~ .")), data = ModelData[["Training"]][,c(taxa_col,loc)], method = 'xgbTree',
#                                 metric = 'ROC', trControl = trControl.xgb,
#                                 tuneGrid = grid.xgb)

#     ## Stopping Cluster
#     stopCluster(cl)
#     ## Resetting to sequential processing
#     registerDoSEQ()
#     saveRDS(xgb,paste0("xgb_",TaskName,".rds"))
#     # xgb <- readRDS(paste0("xgb_",TaskName,".rds"))
#     xgb$bestTune
#     test.pred.prob.xgb <- predict(xgb, newdata = ModelData[[2]], type = "prob")[[TrueLabel]]
                          
# # Actual labels
#     # Actual labels
#     test.actual.classes.xgb <- get(label,ModelData[["Test"]])
    
#     # calculating the ROC Curve
#     roc.curve.xgb <- roc(test.actual.classes.xgb, test.pred.prob.xgb, levels = levels(ModelData[["Test"]][[label]]))
#     # plotting
#     pdf(paste0("roc.curve.xgb.",TaskName,".pdf"))
#     plot(roc.curve.xgb, print.auc = TRUE, print.thres = TRUE)
#     dev.off()

#     saveRDS(roc.curve.xgb,paste0("roc.curve.xgb.",TaskName,".rds"))
#     print("Save the result sucessfully!")
#     # system("rclone copy /nobackup/ngzh5554/Thesis_project/TestOverFitting/v2/roc.curve.nometa.can.ade.xgb.pdf onedrive:HPC_Output/Thesis_project/TestOverFitting/v2/")
#     # calculating AUC
#     auc.value.xgb <- auc(roc.curve.xgb)
#     print(paste("AUC: ", auc.value.xgb))
#     }

XGBoost_Modeling <- function(ModelData = NULL,
                             label = NULL,
                             taxa_col = NULL,
                             TaskName = NULL,
                             TrueLabel = NULL) {
    set.seed(123)
    # Calculate the number of cores
    num_cores <- detectCores() - 1
    num_cores
    print(num_cores)
    # # Register the cores
    # cl <- makePSOCKcluster(no_cores)
    # registerDoParallel(cl)
    # num_cores <- detectCores() - 60
    # num_cores <- 100
    # print(num_cores)
    # num_cores <- min(15, detectCores() - 1) # 避免使用所有核心
    cl <- makePSOCKcluster(num_cores)
    registerDoParallel(cl)
    
    # tuneGrid for xgbTree
    grid.xgb <- expand.grid(nrounds = c(100, 200, 300),
                            max_depth = c(3, 5, 7, 9),
                            eta = c(0.01, 0.1, 0.3),
                            gamma = 0,
                            colsample_bytree = c(0.5, 0.75, 1),
                            min_child_weight = 1,
                            subsample = c(0.5, 0.75, 1))

    colnames(ModelData[["Training"]]) <- make.names(colnames(ModelData[["Training"]])) # nolint
    colnames(ModelData[["Test"]]) <- make.names(colnames(ModelData[["Test"]]))
    loc <- which(label == colnames(ModelData[["Training"]]))

    # Calculate class weights
    class_weights <- table(get(label, ModelData[["Training"]]))
    class_weights <- sum(class_weights) / (length(class_weights) * class_weights)
    sample_weights <- class_weights[get(label, ModelData[["Training"]])]

    # cross-validation control
    trControl.xgb <- trainControl(method = 'repeatedcv', number = 10, repeats = 5,
                                  summaryFunction = twoClassSummary,
                                  classProbs = TRUE,
                                  allowParallel = TRUE)

    # XGBoost Model Building
    set.seed(123)
    xgb <- train(as.formula(paste(label, "~ .")), 
                 data = ModelData[["Training"]][, c(taxa_col, loc)], 
                 method = 'xgbTree',
                 metric = 'ROC', 
                 trControl = trControl.xgb,
                 tuneGrid = grid.xgb,
                 weights = sample_weights)

    ## Stopping Cluster
    stopCluster(cl)
    ## Resetting to sequential processing
    registerDoSEQ()
    saveRDS(xgb, paste0("xgb_", TaskName, ".rds"))
    # xgb <- readRDS(paste0("xgb_", TaskName, ".rds"))
    # xgb <- readRDS("/nobackup/ngzh5554/Thesis_project/Preprocess/Neo_Normalc/xgb_neo.normc_GMPR_XGBoost_all.rds")
    # xgb$bestTune
    test.pred.prob.xgb <- predict(xgb, newdata = ModelData[[2]][,taxa_col], type = "prob")[[TrueLabel]]

    # Actual labels
    test.actual.classes.xgb <- get(label, ModelData[["Test"]])

    # calculating the ROC Curve
    roc.curve.xgb <- roc(test.actual.classes.xgb, test.pred.prob.xgb, levels = levels(as.factor(ModelData[["Test"]][[label]])))
    # plotting
    pdf(paste0("roc.curve.xgb.", TaskName, ".pdf"))
    plot(roc.curve.xgb, print.auc = TRUE, print.thres = TRUE)
    dev.off()

    saveRDS(roc.curve.xgb, paste0("roc.curve.xgb.", TaskName, ".rds"))
    print("Save the result successfully!")
    # calculating AUC
    auc.value.xgb <- auc(roc.curve.xgb)
    print(paste("AUC: ", auc.value.xgb))
}






RunAnalysis <- function(data, label, taxa_col, methods, task_name) {
  normalized_data <- ChooseNormalize(data, taxa_col, methods, task_name)
  model_data <- SplitDataSet(normalized_data, label, task_name, 0.7)
  rf_results <- RF_Modeling(model_data, label, 10, taxa_col, task_name, "Positive")
  return(list(NormalizedData = normalized_data, RF_Results = rf_results))
}

SelectImportantFeatures <- function(model, threshold = 0.01) {
  importance <- varImp(model)$importance
  selected <- rownames(importance)[importance$Overall > threshold]
  return(selected)
}