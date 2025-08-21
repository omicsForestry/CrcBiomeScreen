#' The packaging function for Random Forest modeling
#'
#' @param CrcBiomeScreenObject A \code{CrcBiomeScreenObject} containing normalized microbiome data, sample metadata, etc.
#' @param k.rf Set the number of cross validation
#' @param TaskName A character string used to label the output
#' @param TrueLabel This label is the future prediction target
#' @param num_cores Set the number of the cores in parallel computing
#'
#' @importFrom dplyr mutate across
#' @importFrom foreach foreach %dopar% %do%
#'
#' @return CrcBiomeScreenObject
#' @export
#'
#' @examples
#' \dontrun{
#' CrcBiomeScreenObject <- ModelingRF(
#'                                   CrcBiomeScreenObject = CrcBiomeScreenObject,
#'                                   k.rf = n_cv,
#'                                   TaskName = TaskName,
#'                                   TrueLabel = TrueLabel,
#'                                   num_cores = num_cores)
#' }
ModelingRF <- function(CrcBiomeScreenObject = NULL,
                       k.rf = n_cv,
                       TaskName = NULL,
                       TrueLabel = NULL,
                       num_cores = NULL) {
  # ---- Dependency checks ----
  load_Modeling_deps <- function() {
    pkgs <- c("caret", "foreach", "doParallel", "parallel", "ranger", "pROC")
    for (p in pkgs) {
      if (!requireNamespace(p, quietly = TRUE)) {
        stop(sprintf("The function ModelingRF() requires the '%s' package. Please install it with install.packages('%s').", p, p))
      } else {
        library(p, character.only = TRUE)
      }
    }
    message("All required packages for ModelingRF() are loaded.")
  }
  load_Modeling_deps()
  # ---- Main function logic ----
  set.seed(123)
  folds.rf <- caret::createFolds(CrcBiomeScreenObject$ModelData$TrainLabel, k = k.rf)

  # Calculate the number of cores
  cl <- parallel::makePSOCKcluster(num_cores)
  doParallel::registerDoParallel(cl)

  # tuneGrid for ranger
  grid.rf <- expand.grid(
    mtry = seq(5, 25, by = 5),
    node_size = seq(3, 15, by = 2),
    sample_size = c(.55, .632, .70, .80),
    num.trees = seq(300, 600, by = 100),
    AUC = 0
  )

  # Using ranger random forest for faster implementation
  grid.rf$AUC <- foreach::foreach(
    i = 1:nrow(grid.rf),
    .combine = c,
    .packages = c("ranger", "pROC", "foreach")
  ) %dopar% {
    aucs <- sapply(1:k.rf, function(j) {
      val.indices <- folds.rf[[j]]
      val.data <- CrcBiomeScreenObject$ModelData$Training[val.indices, ]
      train.fold.data <- as.data.frame(CrcBiomeScreenObject$ModelData$Training[-val.indices, ])
      train.fold.data$TrainLabel <- as.factor(CrcBiomeScreenObject$ModelData$TrainLabel[-val.indices])

      # Class weights in each fold
      class_weights <- table(train.fold.data$TrainLabel)
      class_weights <- sum(class_weights) / (length(class_weights) * class_weights)

      # Model training with the specified hyperparameters
      model <- ranger::ranger(
        formula = as.formula(paste("TrainLabel ~ .")),
        data = train.fold.data,
        num.trees = grid.rf$num.trees[i],
        mtry = grid.rf$mtry[i],
        min.node.size = grid.rf$node_size[i],
        sample.fraction = grid.rf$sample_size[i],
        class.weights = class_weights,
        seed = 123,
        classification = TRUE,
        probability = TRUE,
        verbose = FALSE
      )

      # Validation data prediction
      predictions <- predict(model, data = val.data, type = "response")$predictions
      val.Label <- CrcBiomeScreenObject$ModelData$TrainLabel[val.indices]
      roc.obj <- pROC::roc(val.Label, predictions[, TrueLabel])
      pROC::auc(roc.obj)
    })

    mean(aucs) # AUC on the current fold
  }

  # Stop the cluster
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()

  # Choose the best parameters
  best.params.index.rf <- which.max(grid.rf$AUC)
  best.params.rf <- grid.rf[best.params.index.rf, ]
  CrcBiomeScreenObject$ModelResult$RF <- list(grid.para = grid.rf, best.params = best.params.rf)
  attr(CrcBiomeScreenObject$ModelResult$RF, "TaskName") <- TaskName

  return(CrcBiomeScreenObject)
}
