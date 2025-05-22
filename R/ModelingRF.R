
ModelingRF <- function(CrcBiomeScreenObject = NULL,
                        k.rf = n_cv,
                        TaskName = NULL,
                        TrueLabel = NULL,
                        num_cores = NULL) {
  set.seed(123)
  folds.rf <- createFolds(CrcBiomeScreenObject$ModelData$TrainLabel, k = k.rf)
  
  # Calculate the number of cores
  # num_cores <- 10
  # num_cores <- detectCores() - 20
  cl <- makePSOCKcluster(num_cores)
  registerDoParallel(cl)
  
  # tuneGrid for ranger
  grid.rf <- expand.grid(
    mtry = seq(5, 25, by = 5),
    node_size = seq(3, 15, by = 2),
    sample_size = c(.55, .632, .70, .80),
    num.trees = seq(300, 600, by = 100),
    AUC = 0
  )
  
  # Using ranger random forest for faster implementation
  grid.rf$AUC <- foreach(i = 1:nrow(grid.rf), .combine = c, .packages = c("ranger", "pROC")) %dopar% {
    aucs <- sapply(1:k.rf, function(j) {
      val.indices <- folds.rf[[j]]
      val.data <- CrcBiomeScreenObject$ModelData$Training[val.indices, ]
      train.fold.data <- as.data.frame(CrcBiomeScreenObject$ModelData$Training[-val.indices, ])
      train.fold.data$TrainLabel <- as.factor(CrcBiomeScreenObject$ModelData$TrainLabel[-val.indices])

      # Class weights in each fold
      class_weights <- table(train.fold.data$TrainLabel)
      class_weights <- sum(class_weights) / (length(class_weights) * class_weights)

      # Model training with the specified hyperparameters
      model <- ranger(
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
      roc.obj <- roc(val.Label, predictions[, TrueLabel])
      auc(roc.obj)
    })
    
    # AUC on the current fold
    mean(aucs)
  }
  
  # Stop the cluster
  stopCluster(cl)
  registerDoSEQ()
  
  # Choose the best parameters
  best.params.index.rf <- which.max(grid.rf$AUC)
  best.params.rf <- grid.rf[best.params.index.rf, ]
  # Save the best parameters
  CrcBiomeScreenObject$ModelResult$RF <- list(grid.para = grid.rf, best.params = best.params.rf)

  # saveRDS(best.params.rf, paste0("best.params.rf_", TaskName, ".rds"))
  
  print("Save the result successfully!")
  return(CrcBiomeScreenObject)
}



