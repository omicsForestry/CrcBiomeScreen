#' Train the different models
#'
#' @param CrcBiomeScreenObject  A \code{CrcBiomeScreenObject} containing normalized microbiome data, sample metadata, etc.
#' @param model_type Select the method for modeling
#' @param ClassBalance Choose using the class weights or not
#' @param n_cv Set the number of cross validation
#' @param TaskName A character string used to label the output
#' @param TrueLabel This label is the future prediction target
#' @param num_cores Set the number of the cores in parallel computing
#'
#' @importFrom dplyr mutate across
#' @importFrom parallel makePSOCKcluster
#' @importFrom tibble tibble
#' @importFrom foreach %dopar%
#'
#' @return CrcBiomeScreenObject
#' @export
#'
#' @examples CrcBiomeScreenObject <- TrainModels(CrcBiomeScreenObject,
#'                                               model_type = "RF",
#'                                               TaskName = "ToyData_RF",
#'                                               ClassBalance = TRUE,
#'                                               TrueLabel = "CRC",
#'                                               num_cores = 10)
#'
TrainModels <- function(CrcBiomeScreenObject = NULL,
                        model_type = c("RF", "XGBoost"),
                        ClassBalance = TRUE,
                        n_cv = 10,
                        TaskName = NULL,
                        TrueLabel = NULL,
                        num_cores = NULL) {
  # ---- Dependency checks ----
  required_pkgs <- c("caret", "foreach", "parallel", "ranger","xgboost")
  for(pkg in required_pkgs){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
  lapply(required_pkgs, library, character.only = TRUE)

  # For specific model types
  if ("RF" %in% model_type && !requireNamespace("ranger", quietly = TRUE)) {
    stop("The RF model in TrainModels() requires the 'ranger' package. Please install it with install.packages('ranger').")
  }
  if ("XGBoost" %in% model_type && !requireNamespace("xgboost", quietly = TRUE)) {
    stop("The XGBoost model in TrainModels() requires the 'xgboost' package. Please install it with install.packages('xgboost').")
  }

  # ---- Input check ----
  if (is.null(CrcBiomeScreenObject$ModelData)) {
    stop("ModelData is missing in CrcBiomeScreenObject. Please run SplitDataSet first.")
  }

  # ---- Run RF model ----
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

  # ---- Run XGBoost model ----
  if ("XGBoost" %in% model_type) {
    if (ClassBalance) {
      CrcBiomeScreenObject <- ModelingXGBoost(
        CrcBiomeScreenObject = CrcBiomeScreenObject,
        k.rf = n_cv,
        TaskName = TaskName,
        TrueLabel = TrueLabel,
        num_cores = num_cores
      )
    } else {
      CrcBiomeScreenObject <- ModelingXGBoost_noweights(
        CrcBiomeScreenObject = CrcBiomeScreenObject,
        k.rf = n_cv,
        TaskName = TaskName,
        TrueLabel = TrueLabel,
        num_cores = num_cores
      )
    }
  }

  return(CrcBiomeScreenObject)
}
