FilterDataSet <- function(CrcBiomeScreenObject = NULL,
                          label = NULL,
                          condition_col = "study_condition") {
  if (is.null(CrcBiomeScreenObject)) stop("CrcBiomeScreenObject cannot be NULL.")
  if (is.null(label)) stop("Label cannot be NULL.")
  if (!condition_col %in% colnames(CrcBiomeScreenObject$SampleData)) {
    stop(paste("Condition column", condition_col, "not found in SampleData."))
  }

  # Filter the data based on the specified label
  sample_condition <- CrcBiomeScreenObject$SampleData[[condition_col]]
  data <- CrcBiomeScreenObject$NormalizedData[sample_condition %in% label, ]
  sampledata <- CrcBiomeScreenObject$SampleData[sample_condition %in% label, ]

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

  CrcBiomeScreenObject <- FilteredData
  # saveRDS(CrcBiomeScreenObject, paste0("CrcBiomeScreenObject_", "FilteredDataSet", ".rds"))

  return(CrcBiomeScreenObject)
}
