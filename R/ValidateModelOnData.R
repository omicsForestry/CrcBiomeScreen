#' Predict the validation data by using the trained model in CrcBiomeScreenObject
#'
#' @param CrcBiomeScreenObject A CrcBiomeScreenObject containing the model and data to be evaluated.
#' @param model_type The type of model to be evaluated, either "RF" for Random Forest or "XGBoost".
#' @param ValidationData A CrcBiomeScreenObject containing the validation data to be used for model evaluation.
#' @param TaskName A character string used to label the output files and results.
#' @param TrueLabel The true label for the classification task, which is used to evaluate the model's performance.
#' @param condition_col The column name in the SampleData that contains the study condition labels. Default is "study_condition".
#' @param PlotAUC A logical value indicating whether to plot the AUC curve. If TRUE, the AUC curve will be saved as a PDF file.
#'
#' @importFrom tidyr separate
#' @importFrom tibble tibble
#'
#' @return A CrcBiomeScreenObject with the evaluation results stored in the `PredictResult` slot for the specified model type.
#' @export
#'
#' @examples CrcBiomeScreenObject <- ValidateModelOnData(CrcBiomeScreenObject,
#'   model_type = "RF",
#'   ValidationData = ValidationData_filtered_qc,
#'   TaskName = "ValidationData_RF_Validation",
#'   TrueLabel = "CRC",
#'   PlotAUC = TRUE
#' )
#'
ValidateModelOnData <- function(
    CrcBiomeScreenObject = NULL,
    model_type = NULL,
    ValidationData = NULL,
    TaskName = NULL,
    TrueLabel = NULL,
    condition_col = "study_condition",
    PlotAUC = NULL) {
  if (!condition_col %in% colnames(ValidationData$SampleData)) {
    stop(paste("Condition column", condition_col, "not found in SampleData."))
  }
  # Load the model
  if (model_type == "RF") {
    rf.model <- CrcBiomeScreenObject$EvaluateResult$RF$RF.Model
    probs.ValidationData.rf <- predict(rf.model, data = ValidationData$NormalizedData, type = "response")$predictions
    probs.ValidationData.rf.prob <- probs.ValidationData.rf[, TrueLabel]
    rownames(probs.ValidationData.rf) <- rownames(ValidationData$NormalizedData)
    # Actual labels
    actual.classes.rf <- as.factor(ValidationData$SampleData$study_condition)

    # calculating the ROC Curve
    roc.curve.rf <- roc(actual.classes.rf, probs.ValidationData.rf.prob, levels = levels(as.factor(ValidationData$SampleData$study_condition)))

    # Plot
    if (PlotAUC == TRUE) {
      pdf(paste0("roc.curve.rf.", TaskName, ".pdf"))
      plot(roc.curve.rf, print.auc = TRUE, print.thres = TRUE)
      dev.off()
    }
    CrcBiomeScreenObject$PredictResult[["RF"]][[TaskName]] <-
      list(
        # Store the probabilities here
        predictions = probs.ValidationData.rf, 
        roc.curve = roc.curve.rf,
        AUC = auc(roc.curve.rf)
      )
  } else if (model_type == "XGBoost") {
    xgb.model <- CrcBiomeScreenObject$ModelResult$XGBoost$model

    # Test the model
    test.pred.prob.xgb <- predict(xgb.model, newdata = ValidationData$NormalizedData, type = "prob")[[TrueLabel]]

    # Calculate AUC
    roc.curve.xgb <- roc(ValidationData$SampleData$study_condition, test.pred.prob.xgb)
    auc.value.xgb <- auc(roc.curve.xgb)

    # Plot the ROC curve
    if (PlotAUC == TRUE) {
      pdf(paste0("roc.curve.xgb.", TaskName, ".pdf"))
      plot(roc.curve.xgb, print.auc = TRUE, print.thres = TRUE)
      dev.off()
    }

    predictions_df <- data.frame(
      SampleID = rownames(ValidationData$NormalizedData),
      Prediction = test.pred.prob.xgb
    )
    
    CrcBiomeScreenObject$PredictResult[["XGBoost"]][[TaskName]] <-
      list(
        predictions = predictions_df,
        roc.curve = roc.curve.xgb,
        AUC = auc(roc.curve.xgb)
      )
  } else {
    stop("Invalid model type. Please choose either 'RF' or 'XGBoost'.")
  }

  return(CrcBiomeScreenObject)
}
