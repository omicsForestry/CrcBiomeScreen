#' Split the dataset into training and test sets
#'
#' @param CrcBiomeScreenObject From the CreateCrcBiomeScreenObject()
#' @param label Divide the data set by the binary-label
#' @param partition The ratio of dividing the data set
#' @param condition_col The colname of label in SampleData
#'
#' @importFrom dplyr mutate across
#' @importFrom tidyr separate
#'
#' @return A \linkS4class{CrcBiomeScreenObject} with CrcBiomeScreenObject@ModelData
#' @export
#'
#' @examples
#' # Minimal toy object for dataset splitting
#'
#' toy_norm <- data.frame(
#'   TaxaA = runif(10, 10, 30),
#'   TaxaB = runif(10, 5, 15)
#' )
#' rownames(toy_norm) <- paste0("S", 1:10)
#'
#' toy_sampledata <- data.frame(
#'   study_condition = rep(c("control", "CRC"), each = 5),
#'   row.names = paste0("S", 1:10)
#' )
#'
#' # Construct a minimal CrcBiomeScreen object
#' toy_obj <- new(
#'   "CrcBiomeScreen",
#'   AbsoluteAbundance = data.frame(),
#'   RelativeAbundance = data.frame(),
#'   TaxaData = data.frame(),
#'   SampleData = toy_sampledata,
#'   NormalizedData = toy_norm, # <-- IMPORTANT: SplitDataSet needs this
#'   TaxaLevelData = NULL,
#'   OrginalNormalizedData = NULL,
#'   ValidationData = NULL,
#'   ModelData = list(),
#'   ModelResult = NULL,
#'   EvaluateResult = list(),
#'   PredictResult = NULL
#' )
#'
#' # Split into training/testing sets with 70/30 ratio
#' toy_split <- SplitDataSet(
#'   toy_obj,
#'   label = c("control", "CRC"),
#'   partition = 0.7,
#'   condition_col = "study_condition"
#' )
#'
#' # Inspect training labels
#' getModelData(toy_split)$TrainLabel
#'
SplitDataSet <- function(CrcBiomeScreenObject = NULL,
                         label = NULL,
                         partition = NULL,
                         condition_col = "study_condition") {
  # Check if the required parameters are provided
  if (is.null(CrcBiomeScreenObject)) stop("CrcBiomeScreenObject cannot be NULL.")
  if (is.null(label)) stop("Label cannot be NULL.")
  if (is.null(partition) || partition <= 0 || partition >= 1) {
    stop("Partition must be a value between 0 and 1.")
  }
  if (!condition_col %in% colnames(CrcBiomeScreenObject@SampleData)) {
    stop(sprintf("Condition column '%s' not found in SampleData.", condition_col))
  }

  # Select the data based on the label
  sample_condition <- CrcBiomeScreenObject@SampleData[[condition_col]]
  data <- CrcBiomeScreenObject@NormalizedData[sample_condition %in% label, , drop = FALSE]
  sample_condition <- sample_condition[sample_condition %in% label]

  # Check if the data is empty
  if (nrow(data) == 0) stop("No data found for the specified label.")

  # Create training and test set indexes
  withr::with_seed(123, {
    trainIndex <- caret::createDataPartition(sample_condition, p = partition, list = FALSE)
  })

  # Split the data
  train <- data[trainIndex, ]
  test <- data[-trainIndex, ]

  # Add labels to the training and test sets
  train_label <- sample_condition[trainIndex]
  test_label <- sample_condition[-trainIndex]
  #   train[[condition_col]] <- as.factor(train_label)
  #   test[[condition_col]] <- as.factor(test_label)

  # Create a list to store the training and test sets
  ModelData <- list(
    Training = train,
    Test = test,
    TrainLabel = train_label,
    TestLabel = test_label
  )
  colnames(ModelData[["Training"]]) <- make.names(colnames(ModelData[["Training"]]))
  colnames(ModelData[["Test"]]) <- make.names(colnames(ModelData[["Test"]]))

  # Add attributes to the ModelData
  attr(ModelData, "Split Partition") <- partition
  attr(ModelData, "Task Name") <- "SplitDataSet"
  attr(ModelData, "Timestamp") <- Sys.time()
  attr(ModelData, "Training Size") <- nrow(train)
  attr(ModelData, "Test Size") <- nrow(test)

  CrcBiomeScreenObject@ModelData <- ModelData

  return(CrcBiomeScreenObject)
}
