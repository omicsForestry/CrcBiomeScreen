#' Evaluate the performance of model predictions
#'
#' This function calculates performance metrics (e.g., AUC) and plots the
#' ROC curve based on prediction probabilities and true labels.
#'
#' @param predictions A data frame or matrix of model predictions, typically
#'                    containing columns for probability scores for each class.
#' @param true_labels A character vector or factor of the true class labels.
#' @param TrueLabel The positive class label (e.g., "CRC") to use for ROC/AUC calculation.
#' @param TaskName A character string used to label the output files.
#' @param PlotAUC A logical value indicating whether to plot the AUC curve.
#'
#' @importFrom pROC roc auc
#'
#' @return A list containing the ROC curve object and the AUC value.
#' @export
#'
#' @examples
#' # Assume you have obtained predictions and true labels from a previous step
#' # e.g., predictions <- PredictCrcBiomeScreen(my_object, newdata, "RF")
#' true_labels <- my_newdata_object$SampleData$study_condition
#' EvaluateCrcBiomeScreen(predictions, true_labels, TrueLabel = "CRC", TaskName = "MyEvaluation", PlotAUC = TRUE)

EvaluateCrcBiomeScreen <- function(
    predictions,
    true_labels,
    TrueLabel = NULL,
    TaskName = "ModelEvaluation",
    PlotAUC = TRUE) {

  # Check for required inputs
  if (is.null(TrueLabel)) {
    stop("Please specify the 'TrueLabel' (the positive class) for AUC calculation.")
  }

  # Ensure the predictions data frame has a column matching the TrueLabel
  if (!TrueLabel %in% colnames(predictions)) {
    stop(sprintf("The 'predictions' data frame must contain a column named:", TrueLabel))
  }

  # Extract the probabilities for the positive class
  probs <- predictions[, TrueLabel]

  # Ensure labels are a factor
  actual.classes <- as.factor(true_labels)

  # Calculate the ROC curve
  roc.curve <- pROC::roc(
    response = actual.classes,
    predictor = probs,
    levels = levels(actual.classes)
  )

  # Calculate AUC value
  auc.value <- pROC::auc(roc.curve)

  # Plot the ROC curve if requested
  if (PlotAUC) {
    pdf(paste0("roc.curve.", TaskName, ".pdf"))
    plot(roc.curve, print.auc = TRUE, print.thres = TRUE)
    dev.off()
  }

  # Return a list of evaluation results
  evaluation_results <- list(
    roc.curve = roc.curve,
    AUC = auc.value
  )

  return(evaluation_results)
}
