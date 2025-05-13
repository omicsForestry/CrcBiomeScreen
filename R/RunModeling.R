ModelingRF <- function(CrcBiomeScreenObject = NULL,
                        k.rf = NULL,
                        TaskName = NULL,
                        TrueLabel = NULL,
                        num_cores = NULL,
                        repeats = 5) {  
  set.seed(123)
  
  # Calculate the number of cores
  cl <- makePSOCKcluster(num_cores)
  registerDoParallel(cl)
  
  # tuneGrid for ranger
  grid.rf <- expand.grid(
    mtry = seq(5, 25, by = 5),
    node_size = seq(3, 15, by = 2),
    sample_size = c(.55, .632, .70, .80),
    num.trees = seq(600, 1000, by = 150),
    AUC = 0
  )
  
  # Using ranger random forest for faster implementation
  grid.rf$AUC <- foreach(i = 1:nrow(grid.rf), .combine = c, .packages = c("caret","ranger", "pROC")) %dopar% {

    repeat_aucs <- vector("numeric", repeats)
    
    for (repeat_idx in 1:repeats) {
      # Create folds for cross-validation
      folds.rf <- createFolds(CrcBiomeScreenObject$ModelData$TrainLabel, k = k.rf)
      
      # Calculate AUC for each fold
      aucs <- sapply(1:k.rf, function(j) {
        val.indices <- folds.rf[[j]]
        val.data <- CrcBiomeScreenObject$ModelData$Training[val.indices, ]
        train.fold.data <- CrcBiomeScreenObject$ModelData$Training[-val.indices, ]
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
      
      # Save the mean AUC for this parameter set
      repeat_aucs[repeat_idx] <- mean(aucs)
    }
    
    # Return the mean AUC across all repetitions
    mean(repeat_aucs)
  }
  
  # Stop the cluster
  stopCluster(cl)
  registerDoSEQ()
  
  # Choose the best parameters
  best.params.index.rf <- which.max(grid.rf$AUC)
  best.params.rf <- grid.rf[best.params.index.rf, ]
  
  # Save the best parameters
  CrcBiomeScreenObject$ModelResult$RF <- list(grid.para = grid.rf, best.params = best.params.rf)

  print("Save the result successfully!")
  return(CrcBiomeScreenObject)
}

ModelingXGBoost <- function(CrcBiomeScreenObject = NULL,
                            k.rf = NULL,
                             TaskName = NULL,
                             TrueLabel = NULL,
                             num_cores = NULL,
                             repeats = 5) {  # Add repeats parameter for 5 repetitions
  set.seed(123)

  # Check if ModelData exists in the object
  if (is.null(CrcBiomeScreenObject$ModelData)) {
    stop("ModelData is missing in CrcBiomeScreenObject. Please run SplitDataSet first.")
  }
  
  # Set the number of cores for parallel processing
  cl <- makePSOCKcluster(num_cores)
  registerDoParallel(cl)
  
  # Define the grid for hyperparameter tuning
  grid.xgb <- expand.grid(
    nrounds = c(100, 200, 300),
    max_depth = c(3, 5, 7, 9),
    eta = c(0.01, 0.1, 0.3),
    gamma = 0,
    colsample_bytree = c(0.5, 0.75, 1),
    min_child_weight = 1,
    subsample = c(0.5, 0.75, 1),
    AUC = 0
  )
  
  # Prepare training and testing labels
  label_train <- ifelse(CrcBiomeScreenObject$ModelData$TrainLabel == TrueLabel, 1, 0)
  label_test <- ifelse(CrcBiomeScreenObject$ModelData$TestLabel == TrueLabel, 1, 0)

  # Create DMatrix for XGBoost
  dtrain <- xgb.DMatrix(data = as.matrix(CrcBiomeScreenObject$ModelData$Training), label = label_train)
  dtest <- xgb.DMatrix(data = as.matrix(CrcBiomeScreenObject$ModelData$Test), label = label_test)
  
  # Perform hyperparameter tuning with 5 repetitions of 10-fold cross-validation
  for (i in seq_len(nrow(grid.xgb))) {
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = grid.xgb$max_depth[i],
      eta = grid.xgb$eta[i],
      gamma = grid.xgb$gamma[i],
      colsample_bytree = grid.xgb$colsample_bytree[i],
      min_child_weight = grid.xgb$min_child_weight[i],
      subsample = grid.xgb$subsample[i]
    )
    
    # Initialize a vector to store AUCs for each repetition
    repeat_aucs <- numeric(repeats)
    
    for (repeat_idx in 1:repeats) {
      cat("Running repetition", repeat_idx, "of", repeats, "for parameter set", i, "\n")
      
      # Perform 10-fold cross-validation
      cv <- xgb.cv(
        params = params,
        data = dtrain,
        nrounds = grid.xgb$nrounds[i],
        nfold = 10,
        verbose = FALSE,
        prediction = TRUE,
        showsd = TRUE
      )
      
      # Store the maximum AUC for this repetition
      repeat_aucs[repeat_idx] <- max(cv$evaluation_log$test_auc_mean)
    }
    
    # Calculate the mean AUC across all repetitions
    grid.xgb$AUC[i] <- mean(repeat_aucs)
  }

  # Stop the parallel cluster
  stopCluster(cl)
  registerDoSEQ()

  # Select the best parameters
  best.params.index.xgb <- which.max(grid.xgb$AUC)
  best.params.xgb <- grid.xgb[best.params.index.xgb, ]
  

  # Save the results
  CrcBiomeScreenObject$ModelResult$XGBoost <- list(
    grid.xgb = grid.xgb,
    best.params = best.params.xgb
  )
  
  return(CrcBiomeScreenObject)
}

