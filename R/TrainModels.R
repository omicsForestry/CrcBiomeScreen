TrainModels <- function(CrcBiomeScreenObject = NULL,
                        model_type = c("RF", "XGBoost"),
                        ClassBalance = TRUE,
                        n_cv = 10,
                        TaskName = NULL,
                        TrueLabel = NULL,
                        num_cores = NULL) {
  # Check the input
  if (is.null(CrcBiomeScreenObject$ModelData)) {
    stop("ModelData is missing in CrcBiomeScreenObject. Please run SplitDataSet first.")
  }

  # Run the RF model
  if ("RF" %in% model_type) {
    if (ClassBalance) {
      CrcBiomeScreenObject <- ModelingRF(
        CrcBiomeScreenObject = CrcBiomeScreenObject,
        k.rf = n_cv,
        TaskName = TaskName,
        TrueLabel = TrueLabel,
        num_cores = num_cores
      )
    } else {
      CrcBiomeScreenObject <- ModelingRF_noweights(
        CrcBiomeScreenObject = CrcBiomeScreenObject,
        k.rf = n_cv,
        TaskName = TaskName,
        TrueLabel = TrueLabel,
        num_cores = num_cores
      )
    }
  }

  # Run the XGBoost model
  if ("XGBoost" %in% model_type) {
    if (ClassBalance) {
      # 使用加权的 XGBoost 模型
      CrcBiomeScreenObject <- ModelingXGBoost(
        CrcBiomeScreenObject = CrcBiomeScreenObject,
        k.rf = n_cv,
        TaskName = TaskName,
        TrueLabel = TrueLabel,
        num_cores = num_cores
      )
    } else {
      # 使用不加权的 XGBoost 模型
      CrcBiomeScreenObject <- ModelingXGBoost_noweights(
        CrcBiomeScreenObject = CrcBiomeScreenObject,
        k.rf = n_cv,
        TaskName = TaskName,
        TrueLabel = TrueLabel,
        num_cores = num_cores
      )
    }
  }

  # Save the result into the CrcBiomeScreenObject
  saveRDS(CrcBiomeScreenObject, paste0("CrcBiomeScreenObject_", TaskName, ".rds"))

  return(CrcBiomeScreenObject)
}
