#' The packaging function for XGBoost modeling
#'
#' @param CrcBiomeScreenObject A \code{CrcBiomeScreenObject} containing normalized microbiome data, sample metadata, etc.
#' @param k.rf Set the number of cross validation
#' @param repeats Set the number of repeats in cross validation
#' @param TaskName A character string used to label the output
#' @param TrueLabel This label is the future prediction target
#' @param num_cores Set the number of the cores in parallel computing
#'
#' @importFrom dplyr mutate across
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom caret train trainControl twoClassSummary
#'
#' @return A A \code{CrcBiomeScreen} object. with the modelling results.
#' @export
#'
#' @examples
#' # Minimal runnable example for ModelingXGBoost
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
#' # out <- ModelingXGBoost(
#' # CrcBiomeScreenObject = obj,
#' #  k.rf = 2,
#' # TaskName = "toy_XGB",
#' # TrueLabel = c("control", "CRC"),
#' #  num_cores = 1
#' # )
#'
#' obj
#'


ModelingXGBoost <- function(CrcBiomeScreenObject = NULL,
                            k.rf = 10,
                            repeats = 5,
                            TaskName = NULL,
                            TrueLabel = NULL,
                            num_cores = 1) {
  # ---- Parallel setup ----
  num_cores <- as.integer(num_cores)

  if (is.na(num_cores) || num_cores < 1) {
    num_cores <- 1
  }

  if (num_cores > 1) {

    if (.Platform$OS.type == "unix" &&
        nzchar(Sys.getenv("SLURM_JOB_ID"))) {

      cl <- parallel::makeForkCluster(num_cores)

    } else {

      cl <- parallel::makePSOCKcluster(num_cores)

    }

    doParallel::registerDoParallel(cl)

    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
  }

  tune_grid <- expand.grid(
    nrounds = c(100, 200, 300),
    max_depth = c(3, 5, 7, 9),
    eta = c(0.01, 0.1, 0.3),
    gamma = 0,
    colsample_bytree = c(0.5, 0.75, 1),
    min_child_weight = 1,
    subsample = c(0.5, 0.75, 1)
  )

  # Prepare training data
  train_data <- CrcBiomeScreenObject@ModelData$Training
  label_train <- CrcBiomeScreenObject@ModelData$TrainLabel
  label_train <- factor(label_train, levels = unique(CrcBiomeScreenObject@ModelData$TrainLabel))

  class_counts <- table(label_train)
  positive_class <- names(which.max(class_counts))
  negative_class <- names(class_counts)[names(class_counts) != positive_class]

  # Suppose Positive is majority class.  Upweight Negative by ratio:
  w_pos <- 1
  w_neg <- nrow(train_data[label_train == positive_class, ]) /
    nrow(train_data[label_train == negative_class, ])
  weights <- ifelse(label_train == positive_class, w_pos, w_neg)

  # Define caret trainControl
  ctrl <- caret::trainControl(
    method = "repeatedcv",
    number = k.rf,
    repeats = repeats,
    summaryFunction = getFromNamespace("twoClassSummary", "caret"),
    classProbs = TRUE
  )

  train_data <- as.data.frame(train_data)
  train_data$label_train <- as.factor(label_train)

  # Train the model using caret
  # suppressWarnings(): caret internally uses `ntree_limit`, deprecated in xgboost ≥1.6.
  # This does not affect model behavior; warning suppressed for cleaner Bioconductor build logs.
  withr::with_seed(123, {
    old_warn <- getOption("warn")
    options(warn = -1)
    sink(tempfile())
    on.exit({
      sink(NULL)
      options(warn = old_warn)
    }, add = TRUE)

    xgb_method <- .getCaretXgbTreeCompat()
    model_fit <- caret::train(
      label_train ~ .,
      data = train_data,
      method = xgb_method,
      metric = "ROC",
      trControl = ctrl,
      weights = weights,
      tuneGrid = tune_grid
    )
  })

  parallel::stopCluster(cl)
  foreach::registerDoSEQ()

  CrcBiomeScreenObject@ModelResult$XGBoost <- list(
    model = model_fit,
    bestTune = model_fit$bestTune
  )
  attr(CrcBiomeScreenObject@ModelResult$XGBoost, "TaskName") <- TaskName

  return(CrcBiomeScreenObject)
}
