PredictValidation <- function(
  CrcBiomeScreenObject = NULL,
  model_type = NULL,
  ValidationData = NULL,
  TaskName = NULL,
  TrueLabel = NULL,
  condition_col = "study_condition",
  PlotAUC = NULL
) {
  if (!condition_col %in% colnames(ValidationData$SampleData)) {
    stop(paste("Condition column", condition_col, "not found in SampleData."))
  }
    # Load the model
    if(model_type == "RF") {
        rf.model <- CrcBiomeScreenObject$EvaluateResult$RF$RF.Model
        # class_weights <- table(ValidationData$SampleData$study_condition)
        # class_weights <- sum(class_weights) / (length(class_weights) * class_weights)

        probs.ValidationData.rf <- predict(rf.model, data = ValidationData$NormalizedData, type = "response")$predictions
        probs.ValidationData.rf.prob <- probs.ValidationData.rf[, TrueLabel] 
        # Actual labels
        actual.classes.rf <- as.factor(ValidationData$SampleData$study_condition)
        # calculating the ROC Curve
        roc.curve.rf <- roc(actual.classes.rf, probs.ValidationData.rf.prob, levels = levels(as.factor(ValidationData$SampleData$study_condition)))
        # Plot
        if(PlotAUC == TRUE){
            pdf(paste0("roc.curve.rf.",TaskName,".pdf"))
            plot(roc.curve.rf, print.auc = TRUE, print.thres = TRUE)
            dev.off()
        }
        CrcBiomeScreenObject$PredictResult[["RF"]][[TaskName]] <- 
        list(roc.curve = roc.curve.rf,
             AUC = auc(roc.curve.rf))

    } else if(model_type == "XGBoost") {
    
        xgb.model <- CrcBiomeScreenObject$ModelResult$XGBoost$model
        
        # Test the model
        test.pred.prob.xgb <- predict(xgb.model, newdata = ValidationData$NormalizedData, type = "prob")[[TrueLabel]]
        
        # Calculate AUC
        roc.curve.xgb <- roc(ValidationData$SampleData$study_condition, test.pred.prob.xgb)
        auc.value.xgb <- auc(roc.curve.xgb)

        # xgb.model <- CrcBiomeScreenObject$EvaluateResult$XGBoost$XGBoost.Model
        # label_ValidationData <- ifelse(ValidationData$SampleData$study_condition == TrueLabel, 1, 0)
        # Valid_Data <- xgb.DMatrix(data = as.matrix(ValidationData$NormalizedData), label = as.factor(ValidationData$SampleData$study_condition))
        # # Test the model
        # pred.prob.xgb <- predict(xgb.model, newdata = Valid_Data, type = "prob")
  
        # # Calculate AUC
        # roc.curve.xgb <- roc(label_ValidationData, pred.prob.xgb)
        # auc.value.xgb <- auc(roc.curve.xgb)
        
        # Plot the ROC curve
        if(PlotAUC == TRUE){
            pdf(paste0("roc.curve.xgb.",TaskName,".pdf"))
            plot(roc.curve.xgb, print.auc = TRUE, print.thres = TRUE)
            dev.off()
        }
       CrcBiomeScreenObject$PredictResult[["XGBoost"]][[TaskName]] <- 
        list(roc.curve = roc.curve.xgb,
             AUC = auc(roc.curve.xgb))

    } else {
        stop("Invalid model type. Please choose either 'RF' or 'XGBoost'.")
    }
  # Save the result
  # saveRDS(CrcBiomeScreenObject, paste0("CrcBiomeScreenObject_", TaskName, ".rds"))
  # print("Save the result successfully!")
  return(CrcBiomeScreenObject)
}







