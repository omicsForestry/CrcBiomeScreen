#' Run the screening process for the microbiome data
#'
#' @param obj A \code{CrcBiomeScreenObject} containing normalized microbiome data, sample metadata, etc.
#' @param model_type Model type to be used, default is "RF"
#' @param split.requirement A list containing the label and condition column for splitting the dataset, default is NULL
#' @param TaskName A character string used to label the output
#' @param ValidationData A \code{CrcBiomeScreenObject} containing validation data for model evaluation, default is NULL
#' @param TrueLabel The true label for the classification task, which is used to evaluate the model's performance
#' @param num_cores Set the number of cores for parallel computing, default is NULL
#' @param partition The number of partitions for cross-validation
#' @param ClassWeights Whether to use class weights in the model training, default is NULL
#' @param n_cv The number of cross-validation folds, default is NULL
#'
#'
#' @importFrom dplyr mutate across
#' @importFrom tidyr separate
#' @importFrom tibble tibble
#'
#' @return A \linkS4class{CrcBiomeScreenObject} with the results of the screening process, including model training, evaluation, and validation.
#' @export
#'
#' @examples
#' toy_norm <- data.frame(
#'   TaxaA = runif(20),
#'   TaxaB = runif(20)
#' )
#' rownames(toy_norm) <- paste0("Sample", 1:20)
#'
#' toy_sampledata <- data.frame(
#'   study_condition = rep(c("control", "CRC"), each = 10),
#'   row.names = paste0("Sample", 1:20)
#' )
#'
#' toy_obj <- new("CrcBiomeScreen",
#'   SampleData = toy_sampledata,
#'   NormalizedData = toy_norm
#' )
#'
#' if (interactive()) {
#'   result_obj <- RunScreening(
#'     obj = toy_obj,
#'     model_type = "RF",
#'     split.requirement = list(
#'       label = c("control", "CRC"),
#'       condition_col = "study_condition"
#'     ),
#'     TaskName = "Toy_Test",
#'     num_cores = 10,
#'     partition = 0.7,
#'     ValidationData = toy_obj,
#'     TrueLabel = "CRC"
#'   )
#' }
RunScreening <- function(obj,
                         model_type = NULL, # c("RF", "XGBoost")
                         split.requirement = NULL, # list(label = c("control","CRC"),  condition_col = "study_condition")
                         TaskName = TaskName,
                         partition = NULL,
                         ClassWeights = NULL,
                         n_cv = NULL,
                         ValidationData = NULL,
                         TrueLabel = NULL,
                         num_cores = NULL) {
  obj <- SplitDataSet(
    obj,
    label = split.requirement$label,
    condition_col = split.requirement$condition_col,
    partition = partition
  )

  obj <- TrainModels(
    obj,
    model_type = model_type,
    TaskName = TaskName,
    TrueLabel = TrueLabel,
    num_cores = num_cores,
    ClassWeights = ClassWeights,
    n_cv = n_cv
  )

  obj <- EvaluateModel(
    obj,
    model_type = model_type,
    TaskName = paste0(TaskName, "_Test"),
    TrueLabel = TrueLabel,
    PlotAUC = TRUE
  )

  obj <- ValidateModelOnData(
    obj,
    model_type = model_type,
    ValidationData = ValidationData,
    TaskName = paste0(TaskName, "_Validation"),
    TrueLabel = TrueLabel,
    PlotAUC = TRUE
  )

  return(obj)
}
