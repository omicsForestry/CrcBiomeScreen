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
#' # --- Minimal runnable example (no external dependencies) ---
#'
#' # Fake prediction probabilities for 4 samples and 2 classes
#' pred <- data.frame(
#'   control = c(0.8, 0.3, 0.7, 0.2),
#'   CRC     = c(0.2, 0.7, 0.3, 0.8)
#' )
#'
#' # True class labels
#' labels <- factor(c("control", "CRC", "control", "CRC"))
#'
#' # Evaluate performance using CRC as positive class
#' result <- EvaluateCrcBiomeScreen(
#'   predictions = pred,
#'   true_labels = labels,
#'   TrueLabel = "CRC",
#'   PlotAUC = FALSE  # disable plotting for speed/safety
#' )
#'
#' result$AUC

EvaluateCrcBiomeScreen <- function(
    predictions,
    true_labels,
    TrueLabel = NULL,
    TaskName = "ModelEvaluation",
    PlotAUC = FALSE) {   # default FALSE for safety

  if (is.null(TrueLabel)) {
    stop("Please specify the 'TrueLabel' (the positive class).")
  }

  if (!TrueLabel %in% colnames(predictions)) {
    stop(sprintf("The 'predictions' data frame must contain a column named '%s'.", TrueLabel))
  }

  probs <- predictions[[TrueLabel]]
  actual.classes <- as.factor(true_labels)

  roc.curve <- pROC::roc(
    response  = actual.classes,
    predictor = probs,
    levels    = levels(actual.classes)
  )

  auc.value <- pROC::auc(roc.curve)

  if (PlotAUC) {
    pdf(paste0("roc.curve.", TaskName, ".pdf"))
    plot(roc.curve, print.auc = TRUE, print.thres = TRUE)
    dev.off()
  }

  return(list(
    roc.curve = roc.curve,
    AUC = auc.value
  ))
}

