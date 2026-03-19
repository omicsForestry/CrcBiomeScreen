#' Predict the class and probabilities for new data
#'
#' @param CrcBiomeScreenObject The object containing the trained model.
#' @param newdata The data frame or matrix of new features to predict on.
#' @param model_type The type of model to use for prediction ("RF" or "XGBoost").
#'
#' @return A \linkS4class{CrcBiomeScreenObject} with a data frame containing sample-specific predictions.
#' @export
#'
#' @examples
#' # --- Minimal runnable example ---
#'
#' # Create a tiny toy dataset (2 samples × 2 features)
#' newdata <- data.frame(
#'   Feature1 = c(0.2, 0.8),
#'   Feature2 = c(0.7, 0.3)
#' )
#' rownames(newdata) <- c("S1", "S2")
#'
#' # Create a minimal CrcBiomeScreen object
#' toy_obj <- new(
#'   "CrcBiomeScreen",
#'   AbsoluteAbundance = data.frame(),
#'   RelativeAbundance = data.frame(),
#'   TaxaData = data.frame(),
#'   SampleData = data.frame(),
#'   ModelResult = list(),
#'   EvaluateResult = list()
#' )
#'
#' if (interactive()) {
#'   pred_obj <- PredictCrcBiomeScreen(
#'     CrcBiomeScreenObject = toy_obj,
#'     newdata = newdata,
#'     model_type = "RF"
#'   )
#'
#'   print(getPredictResult(pred_obj)$RF)
#' }
PredictCrcBiomeScreen <- function(
  CrcBiomeScreenObject,
  newdata,
  model_type = c("RF", "XGBoost")
) {
  colnames(newdata) <- make.names(colnames(newdata))
  if (model_type == "RF") {
    model <- CrcBiomeScreenObject@EvaluateResult$RF$RF.Model
    probs.ValidationData.rf <- predict(model, data = newdata, type = "response")$predictions
    rownames(probs.ValidationData.rf) <- rownames(newdata)
    predictions <- as.data.frame(probs.ValidationData.rf)
  } else if (model_type == "XGBoost") {
    xgb.model <- CrcBiomeScreenObject@ModelResult$XGBoost$model
    # Test the model
    test.pred.prob.xgb <- predict(xgb.model, newdata = newdata, type = "prob")
    predictions <- test.pred.prob.xgb
    rownames(predictions) <- rownames(newdata)
  }
  CrcBiomeScreenObject@PredictResult[[model_type]] <- predictions
  return(CrcBiomeScreenObject)
}
