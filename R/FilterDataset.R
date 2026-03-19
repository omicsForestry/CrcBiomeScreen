#' Filter the CrcBiomeScreenObject dataset based on a specific label
#'
#' @param CrcBiomeScreenObject A \code{CrcBiomeScreenObject} containing normalized microbiome data and sample metadata.
#' @param label A character vector specifying the label(s) to filter the dataset by.
#' @param condition_col A character string indicating the column in the SampleData that contains the condition labels (default is "study_condition").
#'
#' @return A \linkS4class{CrcBiomeScreenObject} with filtered data based on the specified label.
#' @export
#' @examples
#' # Create toy normalized data (5 samples × 2 taxa)
#' norm_data <- data.frame(
#'   TaxaA = c(10, 20, 15, 30, 10),
#'   TaxaB = c(5, 7, 6, 8, 6)
#' )
#' rownames(norm_data) <- paste0("S", 1:5)
#'
#' # Create sample metadata
#' sample_info <- data.frame(
#'   study_condition = c("control", "CRC", "control", "CRC", "Adenoma"),
#'   country = c("US", "US", "UK", "UK", "US"),
#'   row.names = paste0("S", 1:5),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Construct a minimal CrcBiomeScreen object
#' toy_obj <- new(
#'   "CrcBiomeScreen",
#'   AbsoluteAbundance = data.frame(),
#'   RelativeAbundance = data.frame(),
#'   TaxaData = data.frame(),
#'   SampleData = sample_info,
#'   NormalizedData = norm_data,
#'   TaxaLevelData = NULL,
#'   OrginalNormalizedData = NULL,
#'   ValidationData = NULL,
#'   ModelData = NULL,
#'   ModelResult = NULL,
#'   EvaluateResult = list(),
#'   PredictResult = NULL
#' )
#'
#' # Filter to keep only CRC and control samples
#' filtered_obj <- FilterDataSet(
#'   toy_obj,
#'   label = c("CRC", "control"),
#'   condition_col = "study_condition"
#' )
#'
#' getSampleData(filtered_obj)
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
  data <- CrcBiomeScreenObject@NormalizedData[sample_condition %in% label, ]
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
