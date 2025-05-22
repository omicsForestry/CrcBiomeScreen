
RunScreening <- function(obj, normalize_method = "TSS", model = "RF", split.requirement = c(label = c("control","CRC"), partition = 0.7), TaskName = "GMPR_NHSBCSP", num_cores = 10) {
  obj <- NormalizeData(obj, method = normalize_method)
  obj <- SplitDataSet(obj, split.requirement)

  obj <- TrainModels(obj, model_type = model, TaskName = TaskName, TrueLabel = "CRC", num_cores = num_cores)
  obj <- EvaluateModel(obj, model_type = model, TaskName = paste0(TaskName, "_Test"), TrueLabel = "CRC", PlotAUC = TRUE)
  return(obj)
}














