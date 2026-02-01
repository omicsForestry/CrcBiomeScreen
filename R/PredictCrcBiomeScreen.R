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
#' # Create a minimal CrcBiomeScreen object with a fake RF model
#' # Instead of training, we attach a dummy model object whose predict()
#' # method returns fixed probabilities.
#'
#' fake_rf_model <- structure(
#'   list(),
#'   class = "fakeRF"
#' )
#'
#' # Define a simple predict method for this fake model
#' predict.fakeRF <- function(object, data, type = "response") {
#'   probs <- matrix(
#'     c(0.8, 0.2,   # S1: control=0.8, CRC=0.2
#'       0.3, 0.7),  # S2: control=0.3, CRC=0.7
#'     ncol = 2,
#'     byrow = TRUE
#'   )
#'   colnames(probs) <- c("control", "CRC")
#'   list(predictions = probs)
#' }
#'
#' toy_obj <- new(
#'   "CrcBiomeScreen",
#'   AbsoluteAbundance = data.frame(),
#'   RelativeAbundance = data.frame(),
#'   TaxaData = data.frame(),
#'   SampleData = data.frame(),
#'   ModelResult = list(),
#'   EvaluateResult = list(
#'     RF = list(RF.Model = fake_rf_model)
#'   )
#' )
#'
#' # Run prediction
#' pred_obj <- PredictCrcBiomeScreen(
#'   CrcBiomeScreenObject = toy_obj,
#'   newdata = newdata,
#'   model_type = "RF"
#' )
#'
#' pred_obj@PredictResult$RF


PredictCrcBiomeScreen <- function(
    CrcBiomeScreenObject,
    newdata,
    model_type = c("RF", "XGBoost")) {
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

