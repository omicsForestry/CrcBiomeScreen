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
#' @param partition The number of partitions for cross-validation
#' @param ClassBalance Whether to use class balancing in the model training, default is NULL
#' @param n_cv The number of cross-validation folds, default is NULL
#'
#' @importFrom dplyr mutate across
#' @importFrom tidyr separate
#' @importFrom tibble tibble
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
                         normalize_method = NULL, # c("TSS", "GMPR")
                         model_type = NULL, # c("RF", "XGBoost")
                         split.requirement = NULL, # c(label = c("control","CRC"),  condition_col = "study_condition")
                         TaskName = TaskName,
                         partition = NULL,
                         ClassBalance = NULL,
                         n_cv = NULL,
                         ValidationData = NULL,
                         TrueLabel = NULL,
                         num_cores = NULL) {

  # obj <- NormalizeData(obj, method = normalize_method)
  obj <- SplitDataSet(obj, split.requirement, partition = partition)

  obj <- TrainModels(obj, model_type = model_type, TaskName = TaskName, TrueLabel = TrueLabel, num_cores = num_cores, ClassBalance = ClassBalance, n_cv = n_cv)
  obj <- EvaluateModel(obj, model_type = model_type, TaskName = paste0(TaskName, "_Test"), TrueLabel = TrueLabel, PlotAUC = TRUE)

  obj <- ValidateModelOnData(obj, model_type = model_type, ValidationData = ValidationData, TaskName = paste0(TaskName, "_Validation"), TrueLabel = TrueLabel, PlotAUC = TRUE)

  return(obj)
}
