
#' Random Forest Model
#' Modeling & Cross Validation
#' @Description:
#' This function is for modeling for 
#' random forest and use the cross validation to
#' get the best parameters and model
#'
#' @param ModelData List for the Datasets
#' @param label The label are used for modeling
#' @param k.rf The parameter for cross validation
#' @param taxa_col Select the needed Column for normalize
#' @name TaskName Save the running results
#' 
#' @export 
#' 
#' @return best.params.rf The best parameters for model
#'
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
    num.trees = seq(600, 1000, by = 150),
    AUC = 0
  )
  
  # Using ranger random forest for faster implementation
  grid.rf$AUC <- foreach(i = 1:nrow(grid.rf), .combine = c, .packages = c("ranger", "pROC")) %dopar% {
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


ModelingXGBoost <- function(CrcBiomeScreenObject = NULL,
                             k.rf = n_cv,
                             TaskName = NULL,
                             TrueLabel = NULL,
                             num_cores = NULL) {
  set.seed(123)

  if (is.null(CrcBiomeScreenObject$ModelData)) {
    stop("ModelData is missing in CrcBiomeScreenObject. Please run SplitDataSet first.")
  }
  
  # Set the number of cores
  # num_cores <- detectCores() - 1  # 使用可用核心数减1，避免占满所有核心
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
  
  label_train <- ifelse(CrcBiomeScreenObject$ModelData$TrainLabel == TrueLabel, 1, 0)
  label_test <- ifelse(CrcBiomeScreenObject$ModelData$TestLabel == TrueLabel, 1, 0)

  dtrain <- xgb.DMatrix(data = as.matrix(CrcBiomeScreenObject$ModelData$Training), label = label_train)
  dtest <- xgb.DMatrix(data = as.matrix(CrcBiomeScreenObject$ModelData$Test), label = label_test)
  
  for(i in seq_len(nrow(grid.xgb))) {
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
    # cross-validation
    cv <- xgb.cv(
      params = params,
      data = dtrain,
      nrounds = grid.xgb$nrounds[i],
      nfold = 10,
      verbose = TRUE,
      prediction = TRUE,
      showsd = TRUE
    )
    
    # Store the AUC
    grid.xgb$AUC[i] <- max(cv$evaluation_log$test_auc_mean)
  }

  # Stop the cluster
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


