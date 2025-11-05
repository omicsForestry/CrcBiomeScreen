#' The function for modeling random forest without using class weights
#'
#' @param CrcBiomeScreenObject A \code{CrcBiomeScreenObject} containing normalized microbiome data, sample metadata, etc.
#' @param k.rf Set the number of cross validation
#' @param TaskName A character string used to label the output
#' @param TrueLabel This label is the future prediction target
#' @param num_cores Set the number of the cores in parallel computing
#'
#' @importFrom dplyr mutate across
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom caret createFolds
#'
#' @return CrcBiomeScreenObject
#' @export
#'
#' @examples CrcBiomeScreenObject <- ModelingRF_noweights(
#'   CrcBiomeScreenObject = CrcBiomeScreenObject,
#'   k.rf = n_cv,
#'   TaskName = TaskName,
#'   TrueLabel = TrueLabel,
#'   num_cores = num_cores
#' )
#'
ModelingRF_noweights <- function(CrcBiomeScreenObject = NULL,
                                 k.rf = n_cv,
                                 TaskName = NULL,
                                 TrueLabel = NULL,
                                 num_cores = NULL) {

  # ---- Main function logic ----
  folds.rf <- caret::createFolds(CrcBiomeScreenObject@ModelData$TrainLabel, k = k.rf)

  # Calculate the number of cores
  # num_cores <- 10
  # num_cores <- detectCores() - 20
  cl <- makePSOCKcluster(num_cores)
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
  grid.rf$AUC <- foreach::foreach(i = seq_len(nrow(grid.rf)), .combine = c, .packages = c("ranger", "pROC", "foreach")) %dopar% {
    aucs <- vapply(seq_len(k.rf), function(j) {
      val.indices <- folds.rf[[j]]
      val.data <- CrcBiomeScreenObject@ModelData$Training[val.indices, ]
      train.fold.data <- as.data.frame(CrcBiomeScreenObject@ModelData$Training[-val.indices, ])
      train.fold.data$TrainLabel <- as.factor(CrcBiomeScreenObject@ModelData$TrainLabel[-val.indices])

      # Model training with the specified hyperparameters
      model <- ranger(
        formula = as.formula(paste("TrainLabel ~ .")),
        data = train.fold.data,
        num.trees = grid.rf$num.trees[i],
        mtry = grid.rf$mtry[i],
        min.node.size = grid.rf$node_size[i],
        sample.fraction = grid.rf$sample_size[i],
        seed = 123,
        classification = TRUE,
        probability = TRUE,
        verbose = FALSE
      )
      # Validation data prediction
      predictions <- predict(model, data = val.data, type = "response")$predictions
      val.Label <- CrcBiomeScreenObject@ModelData$TrainLabel[val.indices]
      roc.obj <- roc(val.Label, predictions[, TrueLabel])
      auc(roc.obj)
    }, FUN.VALUE = numeric(1))

    # AUC on the current fold
    mean(aucs)
  }

  # Stop the cluster
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()

  # Choose the best parameters
  best.params.index.rf <- which.max(grid.rf$AUC)
  best.params.rf <- grid.rf[best.params.index.rf, ]
  # Save the best parameters
  CrcBiomeScreenObject@ModelResult$RF_noweights <- list(grid.para = grid.rf, best.params = best.params.rf)
  attr(CrcBiomeScreenObject@ModelResult$RF_noweights, "TaskName") <- TaskName

  return(CrcBiomeScreenObject)
}
