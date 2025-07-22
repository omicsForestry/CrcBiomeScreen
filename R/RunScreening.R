RunScreening <- function(obj,
                         normalize_method = "TSS",
                         model = "RF",
                         split.requirement = NULL, # c(label = c("control","CRC")
                         TaskName = TaskName,
                         ValidationData = NULL,
                         TrueLabel = TrueLabel,
                         num_cores = num_cores) {
  obj <- NormalizeData(obj, method = normalize_method)
  obj <- SplitDataSet(obj, split.requirement, partition = 0.7)

  obj <- TrainModels(obj, model_type = model, TaskName = TaskName, TrueLabel = "CRC", num_cores = num_cores)
  obj <- EvaluateModel(obj, model_type = model, TaskName = paste0(TaskName, "_Test"), TrueLabel = "CRC", PlotAUC = TRUE)

  obj <- ValidateModelOnData(obj, model_type = model, ValidationData = ValidationData, TaskName = paste0(TaskName, "_Validation"), TrueLabel = TrueLabel, PlotAUC = TRUE)

  return(obj)
}
