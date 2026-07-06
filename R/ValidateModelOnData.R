#' Predict the validation data by using the trained model in CrcBiomeScreenObject
#'
#' @param CrcBiomeScreenObject A CrcBiomeScreenObject containing the model and data to be evaluated.
#' @param model_type The type of model to be evaluated, either "RF" for Random Forest or "XGBoost".
#' @param ValidationData A CrcBiomeScreenObject containing the validation data to be used for model evaluation.
#' @param TaskName A character string used to label the output files and results.
#' @param TrueLabel The true label for the classification task, which is used to evaluate the model's performance.
#' @param condition_col The column name in the SampleData that contains the study condition labels. Default is "study_condition".
#' @param PlotAUC A logical value indicating whether to plot the AUC curve. If TRUE, the AUC curve will be saved as a PDF file.
#' @param outdir The output directory where plots will be saved (default: tempdir()).
#'
#' @importFrom tidyr separate
#' @importFrom tibble tibble
#'
#' @return A CrcBiomeScreenObject with the evaluation results stored in the `PredictResult` slot for the specified model type.
#' @export
#'
#' @examples
#' set.seed(123)
#'
#' # -------------------------
#' # Toy taxonomy
#' # -------------------------
#' toy_taxa_strings <- c(
#'   "D_0__Bacteria|D_1__Firmicutes|D_2__Clostridia|D_3__OrderA|D_4__FamilyA|D_5__GenusA",
#'   "D_0__Bacteria|D_1__Firmicutes|D_2__Clostridia|D_3__OrderB|D_4__FamilyB|D_5__GenusB",
#'   "D_0__Bacteria|D_1__Firmicutes|D_2__Clostridia|D_3__OrderC|D_4__FamilyC|D_5__GenusC",
#'   "D_0__Bacteria|D_1__Bacteroidetes|D_2__Bacteroidia|D_3__OrderD|D_4__FamilyD|D_5__GenusD",
#'   "D_0__Bacteria|D_1__Bacteroidetes|D_2__Bacteroidia|D_3__OrderE|D_4__FamilyE|D_5__GenusE",
#'   "D_0__Bacteria|D_1__Proteobacteria|D_2__Gammaproteobacteria|D_3__OrderF|D_4__FamilyF|D_5__GenusF"
#' )
#'
#' toy_taxa <- data.frame(
#'   Taxa = toy_taxa_strings,
#'   stringsAsFactors = FALSE
#' )
#'
#' # -------------------------
#' # Toy training data
#' # -------------------------
#' train_samples <- paste0("S", 1:12)
#'
#' toy_abs <- matrix(
#'   c(
#'     rpois(6 * 6, lambda = 54.8887777),
#'     rpois(6 * 6, lambda = 55)
#'   ),
#'   nrow = 6,
#'   ncol = 12
#' )
#'
#' rownames(toy_abs) <- toy_taxa_strings
#' colnames(toy_abs) <- train_samples
#' toy_abs <- as.data.frame(toy_abs)
#'
#' toy_sample <- data.frame(
#'   number_reads = rep(10000, 12),
#'   study_condition = c(rep("control", 6), rep("CRC", 6)),
#'   row.names = train_samples,
#'   stringsAsFactors = FALSE
#' )
#'
#' obj <- CreateCrcBiomeScreenObject(
#'   AbsoluteAbundance = toy_abs,
#'   TaxaData = toy_taxa,
#'   SampleData = toy_sample
#' )
#'
#' obj <- SplitTaxas(obj)
#' obj <- KeepTaxonomicLevel(obj, level = "Genus")
#' obj <- NormalizeData(obj, method = "TSS", level = "Genus")
#'
#' obj <- SplitDataSet(
#'   obj,
#'   label = c("control", "CRC"),
#'   partition = 0.7
#' )
#'
#' obj <- TrainModels(
#'   obj,
#'   model_type = "RF",
#'   TaskName = "toy_rf",
#'   ClassWeights = FALSE,
#'   TrueLabel = "CRC",
#'   num_cores = 1,
#'   n_cv = 2
#' )
#'
#' obj <- EvaluateModel(
#'   obj,
#'   model_type = "RF",
#'   TaskName = "ToyData_RF_Test",
#'   TrueLabel = "CRC",
#'   PlotAUC = FALSE
#' )
#'
#' # -------------------------
#' # Toy validation data
#' # -------------------------
#' val_samples <- paste0("V", 1:8)
#'
#' val_abund <- matrix(
#'   c(
#'     rpois(6 * 4, lambda = 38),
#'     rpois(6 * 4, lambda = 48)
#'   ),
#'   nrow = 6,
#'   ncol = 8
#' )
#'
#' rownames(val_abund) <- toy_taxa_strings
#' colnames(val_abund) <- val_samples
#' val_abund <- as.data.frame(val_abund)
#'
#' val_sample <- data.frame(
#'   number_reads = rep(10000, 8),
#'   study_condition = c(rep("control", 4), rep("CRC", 4)),
#'   condition = c(rep("control", 4), rep("CRC", 4)),
#'   row.names = val_samples,
#'   stringsAsFactors = FALSE
#' )
#'
#' val_obj <- CreateCrcBiomeScreenObject(
#'   AbsoluteAbundance = val_abund,
#'   TaxaData = toy_taxa,
#'   SampleData = val_sample
#' )
#'
#' val_obj <- SplitTaxas(val_obj)
#' val_obj <- KeepTaxonomicLevel(val_obj, level = "Genus")
#' val_obj <- NormalizeData(val_obj, method = "TSS", level = "Genus")
#'
#' # -------------------------
#' # Align features
#' # -------------------------
#' train_norm <- getNormalizedData(obj)
#' val_norm <- getNormalizedData(val_obj)
#'
#' common_features <- intersect(colnames(train_norm), colnames(val_norm))
#'
#' setNormalizedData(obj) <- train_norm[, common_features, drop = FALSE]
#' setNormalizedData(val_obj) <- val_norm[, common_features, drop = FALSE]
#'
#' # -------------------------
#' # Validate model
#' # -------------------------
#' validated_obj <- ValidateModelOnData(
#'   obj,
#'   ValidationData = val_obj,
#'   model_type = "RF",
#'   TaskName = "toy_validation",
#'   TrueLabel = "CRC",
#'   PlotAUC = FALSE
#' )

ValidateModelOnData <- function(
    CrcBiomeScreenObject = NULL,
    model_type = c("RF", "XGBoost"),
    ValidationData = NULL,
    TaskName = "Validation",
    TrueLabel = NULL,
    condition_col = "study_condition",
    PlotAUC = FALSE,
    outdir = tempdir()) {

  model_type <- match.arg(model_type)

  if (is.null(CrcBiomeScreenObject)) {
    stop("Please provide a CrcBiomeScreenObject containing the trained model.")
  }

  if (is.null(ValidationData)) {
    stop("Please provide ValidationData as a CrcBiomeScreenObject.")
  }

  if (is.null(TrueLabel)) {
    stop("Please specify the 'TrueLabel' positive class.")
  }

  if (!condition_col %in% colnames(ValidationData@SampleData)) {
    stop(sprintf(
      "Condition column '%s' not found in ValidationData@SampleData.",
      condition_col
    ))
  }

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  if (is.null(CrcBiomeScreenObject@PredictResult)) {
    CrcBiomeScreenObject@PredictResult <- list()
  }

  if (is.null(CrcBiomeScreenObject@PredictResult[[model_type]])) {
    CrcBiomeScreenObject@PredictResult[[model_type]] <- list()
  }

  # Prepare validation feature matrix
  validation_data <- as.data.frame(getNormalizedData(ValidationData))
  colnames(validation_data) <- make.names(colnames(validation_data))

  ValidationData@NormalizedData <- validation_data

  # True labels for validation data
  true_labels <- as.factor(ValidationData@SampleData[[condition_col]])

  if (length(true_labels) != nrow(validation_data)) {
    stop(
      "The number of validation labels does not match the number of rows in ",
      "ValidationData@NormalizedData."
    )
  }

  if (model_type == "RF") {
    rf.model <- CrcBiomeScreenObject@EvaluateResult$RF$RF.Model

    if (is.null(rf.model)) {
      stop(
        "No evaluated RF model found in CrcBiomeScreenObject@EvaluateResult$RF$RF.Model. ",
        "Please run EvaluateModel(..., model_type = 'RF') first."
      )
    }

    probs <- predict(
      rf.model,
      data = validation_data,
      type = "response"
    )$predictions

    probs <- as.data.frame(probs)
    rownames(probs) <- rownames(validation_data)

  } else if (model_type == "XGBoost") {
    xgb.model <- CrcBiomeScreenObject@ModelResult$XGBoost$model

    if (is.null(xgb.model)) {
      stop(
        "No XGBoost model found in CrcBiomeScreenObject@ModelResult$XGBoost$model. ",
        "Please run TrainModels(..., model_type = 'XGBoost') first."
      )
    }

    probs <- predict(
      xgb.model,
      newdata = validation_data,
      type = "prob"
    )

    probs <- as.data.frame(probs)
    rownames(probs) <- rownames(validation_data)
  }

  # Evaluate predictions using the shared evaluation function
  eval_result <- EvaluateCrcBiomeScreen(
    predictions = probs,
    true_labels = true_labels,
    TrueLabel = TrueLabel,
    TaskName = TaskName,
    PlotAUC = PlotAUC,
    outdir = outdir
  )

  # Store prediction + evaluation results in PredictResult
  CrcBiomeScreenObject@PredictResult[[model_type]][[TaskName]] <- list(
    predictions = probs,
    roc.curve = eval_result$roc.curve,
    AUC = eval_result$AUC,
    AUC.CI = eval_result$AUC.CI,
    optimal.threshold = eval_result$optimal.threshold,
    F1 = eval_result$F1,
    BalancedAccuracy = eval_result$BalancedAccuracy,
    Precision = eval_result$Precision,
    Recall = eval_result$Recall,
    conf.matrix = eval_result$conf.matrix
  )

  return(CrcBiomeScreenObject)
}
