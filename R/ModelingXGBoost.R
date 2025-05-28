ModelingXGBoost <- function(CrcBiomeScreenObject = NULL,
                                   k.rf = 10,
                                   repeats = 5,
                                   TaskName = NULL,
                                   TrueLabel = NULL,
                                   num_cores = num_cores) {
  
  if (is.null(CrcBiomeScreenObject$ModelData)) {
    stop("ModelData is missing in CrcBiomeScreenObject. Please run SplitDataSet first.")
  }

  set.seed(123)
  # Parallel setup（memory friendly）
  cl <- makePSOCKcluster(num_cores)
  registerDoParallel(cl)
  getDoParWorkers()

  tune_grid <- expand.grid(nrounds = c(100, 200, 300),
                            max_depth = c(3, 5, 7, 9),
                            eta = c(0.01, 0.1, 0.3),
                            gamma = 0,
                            colsample_bytree = c(0.5, 0.75, 1),
                            min_child_weight = 1,
                            subsample = c(0.5, 0.75, 1))
  # Prepare training data
  train_data <- CrcBiomeScreenObject$ModelData$Training
  label_train <- CrcBiomeScreenObject$ModelData$TrainLabel
  # label_train <- factor(label_train, levels = unique(CrcBiomeScreenObject$ModelData$TrainLabel))
  w_pos <- 1
  w_neg <- nrow(train_data[label_train=="Blood_Negative",]) / 
          nrow(train_data[label_train=="Cancer",])
  weights <- ifelse(label_train=="Blood_Negative", w_pos, w_neg)

  # Define caret trainControl
  ctrl <- trainControl(
    method = 'repeatedcv',
    number = k.rf,
    repeats = repeats,
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    allowParallel = TRUE
  )

  # model_weights <- model_weights[CrcBiomeScreenObject$ModelData$TrainLabel]

  train_data <- as.data.frame(train_data)
  train_data$label_train <- label_train

  set.seed(123)
  # Train the model using caret
  model_fit <- train(label_train ~ .,
                data = train_data,
                method = 'xgbTree',
                metric = 'ROC', 
                trControl = ctrl,
                tuneGrid = tune_grid,
                weights = weights,
                verbose = TRUE)

  stopCluster(cl)
  registerDoSEQ()

  CrcBiomeScreenObject$ModelResult$XGBoost <- list(
    model = model_fit,
    bestTune = model_fit$bestTune
  )
  saveRDS(CrcBiomeScreenObject, paste0("CrcBiomeScreenObject_", TaskName, ".rds"))
  return(CrcBiomeScreenObject)
}




