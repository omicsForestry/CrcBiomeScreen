#' Run the screening process for the microbiome data
#'
#' @param obj A \code{CrcBiomeScreenObject} containing normalized microbiome data, sample metadata, etc.
#' @param normalize_method Normalization method to be used, default is "TSS"
#' @param model_type Model type to be used, default is "RF"
#' @param split.requirement A list containing the label and condition column for splitting the dataset, default is NULL
#' @param TaskName A character string used to label the output
#' @param ValidationData A \code{CrcBiomeScreenObject} containing validation data for model evaluation, default is NULL
#' @param TrueLabel The true label for the classification task, which is used to evaluate the model's performance
#' @param num_cores Set the number of cores for parallel computing, default is NULL
#'
#' @return A \code{CrcBiomeScreenObject} with the results of the screening process, including model training, evaluation, and validation.
#' @export
#'
#' @examples CrcBiomeScreenObject <- RunScreening(CrcBiomeScreenObject,
#'                                                normalize_method = "GMPR",
#'                                                model = "RF",
#'                                                split.requirement =
#'                                                c(label = c("control","CRC"),
#'                                                condition_col = "study_condition"),
#'                                                TaskName = "GMPR_NHSBCSP",
#'                                                num_cores = 10,
#'                                                ValidationData = ValidationData,
#'                                                TrueLabel = "Cancer")
#'
RunScreening <- function(obj,
                         normalize_method = normalize_method, # c("TSS", "GMPR")
                         model_type = ModelType, # c("RF", "XGBoost")
                         split.requirement = NULL, # c(label = c("control","CRC"),  condition_col = "study_condition")
                         TaskName = TaskName,
                         ValidationData = NULL,
                         TrueLabel = truelabel,
                         num_cores = num_cores) {
  obj <- NormalizeData(obj, normalize_method = normalize_method)
  obj <- SplitDataSet(obj, split.requirement, partition = 0.7)

  obj <- TrainModels(obj, model_type = ModelType, TaskName = TaskName, TrueLabel = truelabel, num_cores = num_cores)
  obj <- EvaluateModel(obj, model_type = ModelType, TaskName = paste0(TaskName, "_Test"), TrueLabel = truelabel, PlotAUC = TRUE)

  obj <- ValidateModelOnData(obj, model_type = model, ValidationData = ValidationData, TaskName = paste0(TaskName, "_Validation"), TrueLabel = TrueLabel, PlotAUC = TRUE)

  return(obj)
}
