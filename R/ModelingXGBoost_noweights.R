#' The packaging function for XGBoost modeling without using class weights
#'
#' @param CrcBiomeScreenObject A \code{CrcBiomeScreenObject} containing normalized microbiome data, sample metadata, etc.
#' @param k.rf Set the number of cross validation
#' @param repeats Set the number of repeats in cross validation
#' @param TaskName A character string used to label the output
#' @param TrueLabel This label is the future prediction target
#' @param num_cores Set the number of the cores in parallel computing
#' @importFrom dplyr mutate across
#' @importFrom caret train trainControl twoClassSummary
#'
#' @return A A \code{CrcBiomeScreen} object. with the modelling results.
#' @export
#'
#' @examples
#' # Minimal runnable example for ModelingXGBoost_noweights
#'
#' rel_abund <- data.frame(S1 = 10, S2 = 20)
#' rownames(rel_abund) <- "TaxaA"
#'
#' sample_info <- data.frame(
#'   number_reads = c(10000, 12000),
#'   condition = c("control", "CRC"),
#'   row.names = c("S1", "S2")
#' )
#'
#' obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = rel_abund,
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = sample_info
#' )
#'
#'
#' # out <- ModelingXGBoost_noweights(
#' # CrcBiomeScreenObject = obj,
#' # k.rf = 2,
#' # TaskName = "toy_XGB_nw",
#' # TrueLabel = c("control", "CRC"),
#' # num_cores = 1
#' #)
#'
#' obj
ModelingXGBoost_noweights <- function(CrcBiomeScreenObject = NULL,
                                      k.rf = 2,
                                      repeats = 1,
                                      TaskName = NULL,
                                      TrueLabel = NULL,
                                      num_cores = 1) {
  train_data <- as.data.frame(CrcBiomeScreenObject@ModelData$Training)

  label_train <- CrcBiomeScreenObject@ModelData$TrainLabel

  if (is.null(TrueLabel)) {
    stop("TrueLabel must be provided for XGBoost modelling.")
  }

  other_label <- setdiff(unique(label_train), TrueLabel)

  if (length(other_label) != 1) {
    stop("XGBoost currently supports binary classification only.")
  }

  train_data$label_train <- factor(
    label_train,
    levels = c(TrueLabel, other_label)
  )

  if (any(is.na(train_data$label_train))) {
    stop("Some training labels are NA after factor conversion.")
  }

  if (length(levels(train_data$label_train)) != 2) {
    stop("XGBoost currently supports binary classification only.")
  }

  if (k.rf < 2) {
    stop("k.rf must be at least 2.")
  }

  min_class_n <- min(table(train_data$label_train))
  if (k.rf > min_class_n) {
    warning(sprintf(
      "k.rf = %d is larger than the smallest class size (%d). Setting k.rf to %d.",
      k.rf,
      min_class_n,
      min_class_n
    ))
    k.rf <- min_class_n
  }

  set.seed(123)

  folds <- caret::createFolds(
    train_data$label_train,
    k = k.rf,
    returnTrain = TRUE
  )

  ctrl <- caret::trainControl(
    method = "cv",
    number = k.rf,
    index = folds,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    allowParallel = FALSE
  )

  tune_grid <- expand.grid(
    nrounds = 10,
    max_depth = 2,
    eta = 0.1,
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample = 0.8
  )

  xgb_method <- .getCaretXgbTreeCompat()

  model_fit <- caret::train(
    label_train ~ .,
    data = train_data,
    method = xgb_method,
    metric = "ROC",
    trControl = ctrl,
    tuneGrid = tune_grid
  )

  CrcBiomeScreenObject@ModelResult$XGBoost_noweights <- list(
    model = model_fit,
    bestTune = model_fit$bestTune,
    results = model_fit$results
  )

  attr(
    CrcBiomeScreenObject@ModelResult$XGBoost_noweights,
    "TaskName"
  ) <- TaskName

  return(CrcBiomeScreenObject)
}

# ModelingXGBoost_noweights <- function(CrcBiomeScreenObject = NULL,
#                                       k.rf = 10,
#                                       repeats = 5,
#                                       TaskName = NULL,
#                                       TrueLabel = NULL,
#                                       num_cores = num_cores) {
#   # Parallel setup（memory friendly）
#   cl <- makePSOCKcluster(num_cores)
#   doParallel::registerDoParallel(cl)
#
#   # Prepare training data
#   train_data <- CrcBiomeScreenObject@ModelData$Training
#   label_train <- as.factor(CrcBiomeScreenObject@ModelData$TrainLabel)
#   label_train <- factor(label_train, levels = unique(CrcBiomeScreenObject@ModelData$TrainLabel))
#
#   # Define caret trainControl
#   ctrl <- caret::trainControl(
#     method = "repeatedcv",
#     number = k.rf,
#     repeats = repeats,
#     classProbs = TRUE,
#     summaryFunction = getFromNamespace("twoClassSummary", "caret"),
#     allowParallel = TRUE
#   )
#
#   tune_grid <- expand.grid(
#     nrounds = c(100, 200, 300),
#     max_depth = c(3, 5, 7, 9),
#     eta = c(0.01, 0.1, 0.3),
#     gamma = 0,
#     colsample_bytree = c(0.5, 0.75, 1),
#     min_child_weight = 1,
#     subsample = c(0.5, 0.75, 1)
#   )
#
#   train_data <- as.data.frame(train_data)
#   train_data$label_train <- label_train
#
#   # Train the model using caret
#   # suppressWarnings(): caret internally uses `ntree_limit`, deprecated in xgboost ≥1.6.
#   # This does not affect model behavior; warning suppressed for cleaner Bioconductor build logs.
#   withr::with_seed(123, {
#     old_warn <- getOption("warn")
#     options(warn = -1)
#     sink(tempfile())
#     on.exit({
#       sink(NULL)
#       options(warn = old_warn)
#     }, add = TRUE)
#
#     xgb_method <- .getCaretXgbTreeCompat()
#     model_fit <- caret::train(
#       label_train ~ .,
#       data = train_data,
#       method = xgb_method,
#       metric = "ROC",
#       trControl = ctrl,
#       tuneGrid = tune_grid
#     )
#   })
#   #   model_fit <- caret::train(
#   #     label_train ~ .,
#   #     data = train_data,
#   #     method = "xgbTree",
#   #     metric = "ROC",
#   #     trControl = ctrl,
#   #     tuneGrid = tune_grid,
#   #     verbose = FALSE
#   #   )
#   # })
#
#
#   parallel::stopCluster(cl)
#   foreach::registerDoSEQ()
#
#   CrcBiomeScreenObject@ModelResult$XGBoost_noweights <- list(
#     model = model_fit,
#     bestTune = model_fit$bestTune
#   )
#   attr(CrcBiomeScreenObject@ModelResult$XGBoost_noweights, "TaskName") <- TaskName
#
#   return(CrcBiomeScreenObject)
# }
