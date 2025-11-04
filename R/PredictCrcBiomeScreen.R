#' Predict the class and probabilities for new data
#'
#' @param CrcBiomeScreenObject The object containing the trained model.
#' @param newdata The data frame or matrix of new features to predict on.
#' @param model_type The type of model to use for prediction ("RF" or "XGBoost").
#'
#' @return A data frame containing sample-specific predictions.
#' @export
#'
#' @examples
#' ValidationData_filtered_qc$PredictResult$XGBoost <- PredictCrcBiomeScreen(
#' CrcBiomeScreenObject,
#' newdata = ValidationData_filtered_qc$NormalizedData, # Use the appropriate data slot from your new data
#' model_type = "XGBoost")

PredictCrcBiomeScreen <- function(
    CrcBiomeScreenObject,
    newdata,
    model_type = c("RF", "XGBoost")) {

  if (model_type == "RF") {
    model <- CrcBiomeScreenObject@EvaluateResult$RF$RF.Model
    probs.ValidationData.rf <- predict(model, data = newdata, type = "response")$predictions
    rownames(probs.ValidationData.rf) <- rownames(newdata)
    predictions <- data.frame(Probability = probs.ValidationData.rf)
  } else if (model_type == "XGBoost") {
    xgb.model <- CrcBiomeScreenObject@ModelResult$XGBoost$model
    # Test the model
    test.pred.prob.xgb <- predict(xgb.model, newdata = newdata, type = "prob")
    predictions <- test.pred.prob.xgb
    rownames(predictions) <- rownames(newdata)
  }

  return(predictions)
}

