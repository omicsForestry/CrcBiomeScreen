#' Evaluate the model to select the optimal model
#'
#' @param CrcBiomeScreenObject A CrcBiomeScreenObject containing the model data and results
#' @param model_type A character vector indicating the type of model to evaluate. Options are "RF" for Random Forest and "XGBoost" for XGBoost.
#' @param outdir A character string. Path to the output directory where results (PDFs, RDS) should be saved. Defaults to tempdir().
#' @param TaskName A character string used to label the output files and results.
#' @param TrueLabel The true label for the classification task, which is used to evaluate the model's performance.
#' @param PlotAUC A logical value indicating whether to plot the AUC curve. If TRUE, the AUC curve will be saved as a PDF file.
#'
#' @importFrom dplyr mutate across
#' @importFrom pROC roc auc ci.auc coords ci.coords plot.roc
#' @importFrom withr with_seed
#' @importFrom ranger ranger
#'
#' @return A \linkS4class{CrcBiomeScreen} object with with the evaluation results stored in the `EvaluateResult` slot.
#' @export
#' @examples
#' # Minimal runnable example demonstrating EvaluateModel setup
#'
#' # Create small toy relative abundance matrix
#' rel_abund <- data.frame(S1 = 10, S2 = 20)
#' rownames(rel_abund) <- "TaxaA"
#'
#' # Sample metadata (required for CreateCrcBiomeScreenObject)
#' sample_info <- data.frame(
#'   number_reads = c(10000, 12000),
#'   condition = c("control", "CRC"),
#'   row.names = c("S1", "S2")
#' )
#'
#' # Build a minimal CrcBiomeScreen object
#' obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = rel_abund,
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = sample_info
#' )
#'
#' # Add minimal ModelData (required input shape)
#' obj@ModelData <- list(
#'   Training   = data.frame(x = c(1, 2)),
#'   Test       = data.frame(x = c(3, 4)),
#'   TrainLabel = factor(c("control", "CRC")),
#'   TestLabel  = factor(c("control", "CRC"))
#' )
#'
#' # Insert dummy model results so EvaluateModel() can run without training
#' # obj@ModelResult <- list(
#' # RF = list(best.params = list(
#' #   num.trees = 1, mtry = 1, node_size = 1, sample_size = 1
#' # )),
#' #  XGBoost = list(model = list(dummy = TRUE))
#' #)
#'
#' # NOT RUN: Real evaluation requires pROC + ranger etc.
#' # obj <- EvaluateModel(obj, model_type = "RF", TaskName = "toy", TrueLabel = "CRC")
#'
#' obj

EvaluateModel <- function(CrcBiomeScreenObject = NULL,
                          model_type = c("RF", "XGBoost"),
                          outdir = tempdir(),
                          TaskName = NULL,
                          TrueLabel = NULL,
                          PlotAUC = NULL) {
  
  # Ensure output directory exists
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  if (is.null(CrcBiomeScreenObject@EvaluateResult)) {
    CrcBiomeScreenObject@EvaluateResult <- list()
  }
  
  withr::with_seed(123, {
    if ("RF" %in% model_type) {
      CrcBiomeScreenObject <- EvaluateRF(
        CrcBiomeScreenObject = CrcBiomeScreenObject,
        TaskName = TaskName,
        TrueLabel = TrueLabel,
        PlotAUC = PlotAUC,
        outdir = outdir
      )
    } else if ("XGBoost" %in% model_type) {
      CrcBiomeScreenObject <- EvaluateXGBoost(
        CrcBiomeScreenObject = CrcBiomeScreenObject,
        TaskName = TaskName,
        TrueLabel = TrueLabel,
        PlotAUC = PlotAUC,
        outdir = outdir
      )
    } else {
      stop("Invalid model type. Please choose either 'RF' or 'XGBoost'.")
    }
  })
  # Save the result into the CrcBiomeScreenObject
  # CrcBiomeScreenObject@ModelResult <- results
  saveRDS(CrcBiomeScreenObject, file.path(outdir, paste0("CrcBiomeScreenObject_", TaskName, ".rds")))
  # print("Save the result sucessfully!")
  return(CrcBiomeScreenObject)
}



#' Evaluate the Random Forest model
#'
#' @param CrcBiomeScreenObject A CrcBiomeScreenObject containing the model data and results
#' @param TaskName A character string used to label the output files and results.
#' @param TrueLabel The true label for the classification task, which is used to evaluate the model's performance.
#' @param PlotAUC A logical value indicating whether to plot the AUC curve. If TRUE, the AUC curve will be saved as a PDF file.
#' @param outdir The output directory where plots will be saved (default: tempdir()).
#' @importFrom pROC roc auc ci.auc coords ci.coords plot.roc
#' @importFrom caret confusionMatrix
#'
#' @return A CrcBiomeScreenObject with the evaluation results stored in the `EvaluateResult$RF` slot.
#' @export
#' @examples
#' # Minimal runnable example demonstrating input structure for EvaluateRF
#'
#' # Toy training + test matrices
#' train_df <- data.frame(x = c(1, 2), TrainLabel = factor(c("control", "CRC")))
#' test_df  <- data.frame(x = c(3, 4))
#'
#' # Build minimal CrcBiomeScreen object
#' obj <- new("CrcBiomeScreen",
#'   AbsoluteAbundance = data.frame(),
#'   RelativeAbundance = data.frame(),
#'   TaxaData = data.frame(),
#'   SampleData = data.frame(),
#'   ModelData = list(
#'     Training = train_df,
#'     Test = test_df,
#'     TrainLabel = train_df$TrainLabel,
#'     TestLabel = factor(c("control", "CRC"))
#'   ),
#'   ModelResult = list(
#'     RF = list(best.params = list(
#'       num.trees = 1, mtry = 1, node_size = 1, sample_size = 1
#'     ))
#'   )
#' )
#'
#' # NOT RUN: real evaluation uses ranger + pROC (too slow for BioC builds)
#' # out <- EvaluateRF(obj, TaskName = "toy", TrueLabel = "CRC")
#'
#' obj

EvaluateRF <- function(CrcBiomeScreenObject = NULL,
                       outdir = tempdir(),
                       TaskName = NULL,
                       TrueLabel = NULL,
                       PlotAUC = NULL) {
  # Load the best parameters
  best.params <- CrcBiomeScreenObject@ModelResult$RF$best.params
  ModelData <- CrcBiomeScreenObject@ModelData

  ModelData[["Training"]] <- as.data.frame(ModelData[["Training"]])
  ModelData[["Training"]]$TrainLabel <- as.factor(ModelData$TrainLabel)

  # Retraining the model with the best hyperparameters on the entire training dataset
  rf.Model <- ranger::ranger(
    formula         = as.formula(paste("TrainLabel ~ .")),
    data            = ModelData[["Training"]],
    num.trees       = best.params$num.trees,
    mtry            = best.params$mtry,
    min.node.size   = best.params$node_size,
    sample.fraction = best.params$sample_size,
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
  auc.value.rf <- auc(roc.curve.rf)
  # Confidence Interval
  pROC::ci.auc(roc.curve.rf, conf.level = 0.95, method = "delong")
  # Finding Optimal threshold using Youden's Index
  coords.rf <- pROC::coords(roc.curve.rf, "best", ret = "all", best.method = "youden")
  optimal.threshold.rf <- coords.rf$threshold
  optimal.threshold.rf
  label <- levels(as.factor(ModelData$TestLabel))
  # Making predictions based on optimal threshold
  test.class.predictions.rf <- as.factor(ifelse(test.pred.prob.rf >= optimal.threshold.rf, TrueLabel,
    label[!label %in% TrueLabel]
  ))
  # Confusion Matrix
  conf.matrix.rf <- caret::confusionMatrix(test.class.predictions.rf, as.factor(ModelData$TestLabel), positive = TrueLabel)
  # F1-score
  f1_score.rf <- conf.matrix.rf$byClass["F1"]

  # Balanced Accuracy
  balanced_accuracy.rf <- conf.matrix.rf$byClass["Balanced Accuracy"]

  # Precision
  precision.rf <- conf.matrix.rf$byClass["Precision"]

  # Recall
  recall.rf <- conf.matrix.rf$byClass["Recall"]

  CrcBiomeScreenObject@EvaluateResult$RF <-
    list(
      predictions = test.predictions.rf,
      roc.curve = roc.curve.rf,
      AUC = auc.value.rf,
      F1 = f1_score.rf,
      BalancedAccuracy = balanced_accuracy.rf,
      Precision = precision.rf,
      Recall = recall.rf,
      RF.Model = rf.Model,
      conf.matrix = conf.matrix.rf
    )

  # Plot
  if (PlotAUC == TRUE) {
    pdf(file.path(outdir,paste0("roc.curve.rf.", TaskName, ".pdf")))
    plot(roc.curve.rf, print.auc = TRUE, print.thres = TRUE)
    dev.off()
  }

  return(CrcBiomeScreenObject)
}


#' Evaluate the XGBoost model
#'
#' @param CrcBiomeScreenObject A CrcBiomeScreenObject containing the model data and results
#' @param TaskName A character string used to label the output files and results.
#' @param TrueLabel The true label for the classification task, which is used to evaluate the model's performance.
#' @param PlotAUC A logical value indicating whether to plot the AUC curve. If TRUE, the AUC curve will be saved as a PDF file.
#' @param outdir The output directory where plots will be saved (default: tempdir()).
#'
#' @importFrom pROC roc auc ci.auc coords ci.coords plot.roc
#' @importFrom caret confusionMatrix
#'
#' @return A CrcBiomeScreenObject with the evaluation results stored in the `EvaluateResult$XGBoost` slot.
#' @export
#' @examples
#' # Minimal runnable example demonstrating input structure for EvaluateXGBoost
#'
#' # Toy data for testing
#' train_df <- data.frame(x = c(1, 2), TrainLabel = factor(c("control", "CRC")))
#' test_df  <- data.frame(x = c(3, 4))
#'
#' obj <- new("CrcBiomeScreen",
#'   AbsoluteAbundance = data.frame(),
#'   RelativeAbundance = data.frame(),
#'   TaxaData = data.frame(),
#'   SampleData = data.frame(),
#'   ModelData = list(
#'     Training = train_df,
#'     Test = test_df,
#'     TrainLabel = train_df$TrainLabel,
#'     TestLabel = factor(c("control", "CRC"))
#'   ),
#'   ModelResult = list(
#'     XGBoost = list(
#'       model = list(dummy_model = TRUE)  # placeholder model
#'     )
#'   )
#' )
#'
#' # NOT RUN: actual evaluation needs xgboost + pROC
#' # out <- EvaluateXGBoost(obj, TaskName = "toy", TrueLabel = "CRC")
#'
#' obj

EvaluateXGBoost <- function(CrcBiomeScreenObject = NULL,
                            outdir = tempdir(),
                            TaskName = NULL,
                            TrueLabel = NULL,
                            PlotAUC = NULL) {
  xgb.model <- CrcBiomeScreenObject@ModelResult$XGBoost$model

  # Test the model
  test.pred.prob.xgb <- predict(xgb.model, newdata = CrcBiomeScreenObject@ModelData$Test, type = "prob")[[TrueLabel]]

  # Calculate AUC
  roc.curve.xgb <- roc(CrcBiomeScreenObject@ModelData$TestLabel, test.pred.prob.xgb)
  auc.value.xgb <- auc(roc.curve.xgb)

  # Optimal threshold using Youden's Index
  coords.xgb <- pROC::coords(roc.curve.xgb, "best", ret = "all", best.method = "youden")
  optimal.threshold.xgb <- coords.xgb$threshold
  label <- levels(as.factor(CrcBiomeScreenObject@ModelData$TestLabel))
  test.class.predictions.xgb <- as.factor(ifelse(test.pred.prob.xgb >= optimal.threshold.xgb, TrueLabel,
    label[!label %in% TrueLabel]
  ))

  # Confusion Matrix
  conf.matrix.xgb <- caret::confusionMatrix(test.class.predictions.xgb, as.factor(CrcBiomeScreenObject@ModelData$TestLabel), positive = TrueLabel)
  # F1-score
  f1_score.xgb <- conf.matrix.xgb$byClass["F1"]

  # Balanced Accuracy
  balanced_accuracy.xgb <- conf.matrix.xgb$byClass["Balanced Accuracy"]

  # Precision
  precision.xgb <- conf.matrix.xgb$byClass["Precision"]

  # Recall
  recall.xgb <- conf.matrix.xgb$byClass["Recall"]

  # Plot the ROC curve
  if (PlotAUC == TRUE) {
    pdf(file.path(outdir,"roc.curve.xgb.", TaskName, ".pdf"))
    plot(roc.curve.xgb, print.auc = TRUE, print.thres = TRUE)
    dev.off()
  }

  # Save results
  CrcBiomeScreenObject@EvaluateResult$XGBoost <-
    list(
      predictions = test.pred.prob.xgb,
      roc.curve = roc.curve.xgb,
      AUC = auc.value.xgb,
      F1 = f1_score.xgb,
      ConfusionMatrix = conf.matrix.xgb,
      BalancedAccuracy = balanced_accuracy.xgb,
      Precision = precision.xgb,
      Recall = recall.xgb,
      XGBoost.Model = xgb.model
    )

  return(CrcBiomeScreenObject)
}
