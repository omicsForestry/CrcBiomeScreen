
TrainModels <- function(CrcBiomeScreenObject = NULL,
                        model_type = c("RF", "XGBoost"), 
                        n_cv = 10, 
                        TaskName = NULL, 
                        TrueLabel = NULL,
                        num_cores = NULL){
  # Check the input
  if (is.null(CrcBiomeScreenObject$ModelData)) {
    stop("ModelData is missing in CrcBiomeScreenObject. Please run SplitDataSet first.")
  }
  
  # Run the RF model
  if ("RF" %in% model_type){
    CrcBiomeScreenObject <- ModelingRF(
      CrcBiomeScreenObject = CrcBiomeScreenObject,
      k.rf = n_cv,
      TaskName = TaskName,
      TrueLabel = TrueLabel,
      num_cores = num_cores)
  }
  else if ("XGBoost" %in% model_type) {
  # Run the XGBoost model
    CrcBiomeScreenObject <- ModelingXGBoost(
      CrcBiomeScreenObject = CrcBiomeScreenObject,
      k.rf = n_cv,
      TaskName = TaskName,
      TrueLabel = TrueLabel,
      num_cores = num_cores)
  }
  # Save the result into the CrcBiomeScreenObject
  # CrcBiomeScreenObject$ModelResult <- results
  saveRDS(CrcBiomeScreenObject, paste0("CrcBiomeScreenObject_", TaskName, ".rds"))
  
  return(CrcBiomeScreenObject)
}

