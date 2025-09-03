test_that("TrainModels works correctly", {
  # Load the CrcBiomeScreenObject
  CrcBiomeScreenObject <- readRDS("CrcBiomeScreenObject.rds")

  # Define parameters
  model_type <- c("RF", "XGBoost")
  ClassBalance <- TRUE
  n_cv <- 10
  TaskName <- "TestTask"
  TrueLabel <- "CRC"
  num_cores <- 5

  # Call the TrainModels function
  result_RF <- TrainModels(
    CrcBiomeScreenObject = CrcBiomeScreenObject,
    model_type = model_type[1],
    ClassBalance = ClassBalance,
    n_cv = n_cv,
    TaskName = TaskName,
    TrueLabel = TrueLabel,
    num_cores = num_cores
  )

  result_XGBoost <- TrainModels(
    CrcBiomeScreenObject = CrcBiomeScreenObject,
    model_type = model_type[2],
    ClassBalance = ClassBalance,
    n_cv = n_cv,
    TaskName = TaskName,
    TrueLabel = TrueLabel,
    num_cores = num_cores
  )

  # Check if the result is not NULL and contains expected elements
  expect_true(!is.null(result_RF))
  expect_true(!is.null(result_XGBoost))
  expect_true("ModelResult" %in% names(result_RF))
  expect_true("ModelData" %in% names(result_RF))
  expect_true("ModelResult" %in% names(result_XGBoost))
  expect_true("ModelData" %in% names(result_XGBoost))

  saveRDS(result_RF, "result_RF.rds")
  saveRDS(result_XGBoost, "result_XGBoost.rds")

  # Define parameters
  model_type <- c("RF", "XGBoost")
  ClassBalance <- FALSE
  n_cv <- 10
  TaskName <- "TestTask_noweights"
  TrueLabel <- "CRC"
  num_cores <- 5

  # Call the TrainModels function
  result_RF_noweights <- TrainModels(
    CrcBiomeScreenObject = CrcBiomeScreenObject,
    model_type = model_type[1],
    ClassBalance = ClassBalance,
    n_cv = n_cv,
    TaskName = TaskName,
    TrueLabel = TrueLabel,
    num_cores = num_cores
  )

  result_XGBoost_noweights <- TrainModels(
    CrcBiomeScreenObject = CrcBiomeScreenObject,
    model_type = model_type[2],
    ClassBalance = ClassBalance,
    n_cv = n_cv,
    TaskName = TaskName,
    TrueLabel = TrueLabel,
    num_cores = num_cores
  )

  # Check if the result is not NULL and contains expected elements
  expect_true(!is.null(result_RF_noweights))
  expect_true(!is.null(result_XGBoost_noweights))
  expect_true("ModelResult" %in% names(result_RF_noweights))
  expect_true("ModelData" %in% names(result_RF_noweights))
  expect_true("ModelResult" %in% names(result_XGBoost_noweights))
  expect_true("ModelData" %in% names(result_XGBoost_noweights))

  saveRDS(result_RF_noweights, "result_RF_noweights.rds")
  saveRDS(result_XGBoost_noweights, "result_XGBoost_noweights.rds")
})
