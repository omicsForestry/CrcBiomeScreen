
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
EvaluateModel <- function(CrcBiomeScreenObject = NULL, 
                      model_type = c("RF", "XGBoost"),
                      TaskName = NULL,
                      TrueLabel = NULL,
                      PlotAUC = NULL){
    if (is.null(CrcBiomeScreenObject$EvaluateResult)) {
        CrcBiomeScreenObject$EvaluateResult <- list()
}
    set.seed(123)
    if ("RF" %in% model_type){
       CrcBiomeScreenObject <- 
       EvaluateRF(CrcBiomeScreenObject = CrcBiomeScreenObject, 
                      TaskName = TaskName,
                      TrueLabel = TrueLabel,
                      PlotAUC = PlotAUC)

    } else if ("XGBoost" %in% model_type) {
        CrcBiomeScreenObject <- 
        EvaluateXGBoost(CrcBiomeScreenObject = CrcBiomeScreenObject, 
                      TaskName = TaskName,
                      TrueLabel = TrueLabel,
                      PlotAUC = PlotAUC)
    } else {
        stop("Invalid model type. Please choose either 'RF' or 'XGBoost'.")
    }
    # Save the result into the CrcBiomeScreenObject
    # CrcBiomeScreenObject$ModelResult <- results
    saveRDS(CrcBiomeScreenObject, paste0("CrcBiomeScreenObject_", TaskName, ".rds"))
    # print("Save the result sucessfully!")
    return(CrcBiomeScreenObject)
}



EvaluateRF <- function(CrcBiomeScreenObject = NULL, 
                      TaskName = NULL,
                      TrueLabel = NULL,
                      PlotAUC = NULL){
    # Load the best parameters
    best.params <- CrcBiomeScreenObject$ModelResult$RF$best.params
    ModelData <- CrcBiomeScreenObject$ModelData
    
    # Class weights in each fold
    class_weights <- table(ModelData$TestLabel)
    class_weights <- sum(class_weights) / (length(class_weights) * class_weights)
    
    ModelData[["Training"]]$TrainLabel <- as.factor(ModelData$TrainLabel)
    # Retraining the model with the best hyperparameters on the entire training dataset
    rf.Model <- ranger(
        formula         = as.formula(paste("TrainLabel ~ .")),
        data            = ModelData[["Training"]], 
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
    test.predictions.rf <- predict(rf.Model, data = ModelData[["Test"]], type = "response")$predictions
    test.pred.prob.rf <- test.predictions.rf[, TrueLabel] 
    # Actual labels
    test.actual.classes.rf <- as.factor(ModelData$TestLabel)

    # calculating the ROC Curve
    roc.curve.rf <- roc(test.actual.classes.rf, test.pred.prob.rf, levels = levels(as.factor(ModelData$TestLabel)))
    
    CrcBiomeScreenObject$EvaluateResult$RF <- 
    list(roc.curve = roc.curve.rf,
         AUC = auc(roc.curve.rf),
         RF.Model = rf.Model)

    # Plot
    if(PlotAUC == TRUE){
        pdf(paste0("roc.curve.rf.",TaskName,".pdf"))
        plot(roc.curve.rf, print.auc = TRUE, print.thres = TRUE)
        dev.off()
    }

    # Save the result
    saveRDS(roc.curve.rf,paste0("roc.curve.rf.",TaskName,".rds"))
    # print("Save the result sucessfully!")

    return(CrcBiomeScreenObject)
}


EvaluateXGBoost <- function(CrcBiomeScreenObject = NULL, 
                      TaskName = NULL,
                      TrueLabel = NULL,
                      PlotAUC = NULL){
                      
    best.params.xgb <- CrcBiomeScreenObject$ModelResult$XGBoost$best.params 
    # Train the final model with the best parameters
    best.params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        max_depth = best.params.xgb$max_depth,
        eta = best.params.xgb$eta,
        gamma = best.params.xgb$gamma,
        colsample_bytree = best.params.xgb$colsample_bytree,
        min_child_weight = best.params.xgb$min_child_weight,
        subsample = best.params.xgb$subsample
    )
  
    label_train <- ifelse(CrcBiomeScreenObject$ModelData$TrainLabel == TrueLabel, 1, 0)
    label_test <- ifelse(CrcBiomeScreenObject$ModelData$TestLabel == TrueLabel, 1, 0)
    dtrain <- xgb.DMatrix(data = as.matrix(CrcBiomeScreenObject$ModelData$Training), label = label_train)
    dtest <- xgb.DMatrix(data = as.matrix(CrcBiomeScreenObject$ModelData$Test), label = label_test)
    
    xgb.model <- xgb.train(
        params = best.params,
        data = dtrain,
        nrounds = best.params.xgb$nrounds,
        verbose = FALSE
    )
    
  # Test the model
  test.pred.prob.xgb <- predict(xgb.model, newdata = dtest, type = "prob")
  
  # Calculate AUC
  roc.curve.xgb <- roc(label_test, test.pred.prob.xgb)
  auc.value.xgb <- auc(roc.curve.xgb)
  
    # Plot the ROC curve
    if(PlotAUC == TRUE){
        pdf(paste0("roc.curve.xgb.",TaskName,".pdf"))
        plot(roc.curve.xgb, print.auc = TRUE, print.thres = TRUE)
        dev.off()
    }
    CrcBiomeScreenObject$EvaluateResult$XGBoost <- 
    list(roc.curve = roc.curve.xgb,
         AUC = auc(roc.curve.xgb),
         XGBoost.Model = xgb.model)
    # Save the result
    saveRDS(roc.curve.xgb,paste0("roc.curve.xgb.",TaskName,".rds"))
    print("Save the result sucessfully!")
    return(CrcBiomeScreenObject)
}
