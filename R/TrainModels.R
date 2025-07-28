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
#' @importFrom caret createFolds
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
  # Check the input
  if (is.null(CrcBiomeScreenObject$ModelData)) {
    stop("ModelData is missing in CrcBiomeScreenObject. Please run SplitDataSet first.")
  }

  # Run the RF model
  if ("RF" %in% model_type) {
    # Using the RF with class weights
    if (ClassBalance == TRUE) {
      CrcBiomeScreenObject <- ModelingRF(
        CrcBiomeScreenObject = CrcBiomeScreenObject,
        k.rf = n_cv,
        TaskName = TaskName,
        TrueLabel = TrueLabel,
        num_cores = num_cores
      )
    } else {
      # Using the RF without class weights
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
    if (ClassBalance == TRUE) {
      # Using the XGBoost with class weights
      CrcBiomeScreenObject <- ModelingXGBoost(
        CrcBiomeScreenObject = CrcBiomeScreenObject,
        k.rf = n_cv,
        TaskName = TaskName,
        TrueLabel = TrueLabel,
        num_cores = num_cores
      )
    } else {
      # Using the XGBoost without class weights
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
