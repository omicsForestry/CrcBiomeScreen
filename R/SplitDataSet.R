#' Split the dataset into training and test sets
#'
#' @param CrcBiomeScreenObject From the CreateCrcBiomeScreenObject()
#' @param label Divide the data set by the binary-label
#' @param partition The ratio of dividing the data set
#' @param condition_col The colname of label in SampleData
#'
#' @return CrcBiomeScreenObject$ModelData
#' @export
#'
#' @examples CrcBiomeScreenObject <- SplitDataSet(CrcBiomeScreenObject,
#'                                                label = c("control","CRC"),
#'                                                partition = 0.7,
#'                                                condition_col = "study_condition")
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
  if (!condition_col %in% colnames(CrcBiomeScreenObject$SampleData)) {
    stop(paste("Condition column", condition_col, "not found in SampleData."))
  }

  # Select the data based on the label
  sample_condition <- CrcBiomeScreenObject$SampleData[[condition_col]]
  data <- CrcBiomeScreenObject$NormalizedData[sample_condition %in% label, ]
  sample_condition <- sample_condition[sample_condition %in% label]

  # Check if the data is empty
  if (nrow(data) == 0) stop("No data found for the specified label.")

  # Create training and test set indexes
  set.seed(123)
  trainIndex <- createDataPartition(sample_condition, p = partition, list = FALSE)

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

  CrcBiomeScreenObject$ModelData <- ModelData

  return(CrcBiomeScreenObject)
}
