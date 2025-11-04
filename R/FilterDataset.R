#' Filter the CrcBiomeScreenObject dataset based on a specific label
#'
#' @param CrcBiomeScreenObject A \code{CrcBiomeScreenObject} containing normalized microbiome data and sample metadata.
#' @param label A character vector specifying the label(s) to filter the dataset by.
#' @param condition_col A character string indicating the column in the SampleData that contains the condition labels (default is "study_condition").
#'
#' @return A \code{CrcBiomeScreenObject} with filtered data based on the specified label.
#' @export
#'
#' @examples ValidationData_filtered <- FilterDataSet(ValidationData,
#'   label = c("CRC", "control"),
#'   condition_col = "study_condition"
#' )
#'
FilterDataSet <- function(CrcBiomeScreenObject = NULL,
                          label = NULL,
                          condition_col = "study_condition") {

  if (is.null(CrcBiomeScreenObject)) stop("CrcBiomeScreenObject cannot be NULL.")
  if (is.null(label)) stop("Label cannot be NULL.")
  if (!condition_col %in% colnames(CrcBiomeScreenObject@SampleData)) {
    stop(sprintf("Condition column", condition_col, "not found in SampleData."))
  }

  # Filter the data based on the specified label
  sample_condition <- CrcBiomeScreenObject@SampleData[[condition_col]]
  data <- CrcBiomeScreenObject@NormalizedData[sample_condition %in% label,]
  sampledata <- CrcBiomeScreenObject@SampleData[sample_condition %in% label, ]

  # Checck if any data is found
  if (nrow(data) == 0) stop("No data found for the specified label.")

  FilteredData <- list(
    NormalizedData = data,
    SampleData = sampledata
  )


  # Add attributes to the filtered data
  attr(FilteredData, "Task Name") <- "FilterDataSet"
  attr(FilteredData, "Timestamp") <- Sys.time()
  attr(FilteredData, "Filtered Size") <- nrow(data)

  CrcBiomeScreenObject@NormalizedData <- FilteredData$NormalizedData
  CrcBiomeScreenObject@SampleData <- FilteredData$SampleData

  return(CrcBiomeScreenObject)
}
