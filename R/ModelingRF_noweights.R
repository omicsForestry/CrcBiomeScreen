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
#' @return A A \code{CrcBiomeScreen} object. with the modelling results.
#' @export
#'
#' @examples
#' # Minimal runnable example for ModelingRF_noweights
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
#' # out <- ModelingRF_noweights(
#' #   CrcBiomeScreenObject = obj,
#' #   k.rf = 2,
#' #   TaskName = "toy_RF_nw",
#' #   TrueLabel = c("control", "CRC"),
#' #   num_cores = 1
#' # )
#'
#' obj

ModelingRF_noweights <- function(CrcBiomeScreenObject = NULL,
                                 k.rf = n_cv,
                                 TaskName = NULL,
                                 TrueLabel = NULL,
                                 num_cores = NULL) {

  # ---- Main function logic ----
  folds.rf <- caret::createFolds(CrcBiomeScreenObject@ModelData$TrainLabel, k = k.rf)

  p <- ncol(CrcBiomeScreenObject@ModelData$Training)
  mtry_vals <- unique(pmin(c(1, 2, 3, 5, 10, 15, 20, 25), p))
  mtry_vals <- mtry_vals[mtry_vals >= 1]

  # Calculate the number of cores
  cl <- makePSOCKcluster(num_cores)
  doParallel::registerDoParallel(cl)

  # tuneGrid for ranger
  grid.rf <- expand.grid(
    mtry = mtry_vals,
    node_size = seq(3, 15, by = 2),
    sample_size = c(.55, .632, .70, .80),
    num.trees = seq(300, 600, by = 100),
    AUC = 0
  )

  p <- ncol(CrcBiomeScreenObject@ModelData$Training)
  if (max(grid.rf$mtry) > p) {
    warning(paste("mtry grid exceeds number of features (p =", p, "), clipping mtry to p."))
    grid.rf <- grid.rf[grid.rf$mtry <= p, ]
  }

  # Using ranger random forest for faster implementation
  grid.rf$AUC <- foreach::foreach(i = seq_len(nrow(grid.rf)), .combine = c, .packages = c("ranger", "pROC", "foreach")) %do% {
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
