#' Normalise the absolute data to relative data by using Total Sum Scaling
#'  and Geometric Mean of Pairwise Ratios (GMPR)
#'
#' @param CrcBiomeScreenObject From the CreateCrcBiomeScreenObject()
#' @param method "TSS" or "GMPR"
#'
#'
#' @return Updated CrcBiomeScreenObject with NormalizedData
#' @export
#'
#' @examples
#' # Normalize using GMPR (Geometric Mean of Pairwise Ratios)
#' CrcBiomeScreenObject <- NormalizeData(CrcBiomeScreenObject, method = "GMPR", level = "Genus")
#' # Normalize using TSS (Total Sum Scaling)
#' CrcBiomeScreenObject <- NormalizeData(CrcBiomeScreenObject, method = "TSS", level = "Genus")

NormalizeData <- function(CrcBiomeScreenObject = NULL, method = NULL, level = NULL) {
  Data <- CrcBiomeScreenObject@TaxaLevelData[[paste0(level, "LevelData")]]

  if (method == "TSS") {
    # Calculate the total number of counts in each sample
    # total_counts = rowSums(Data)

    # convert the absolute abundance to relative abundance
    Data <- Data / rowSums(Data)
  } else if (method == "GMPR") {
    if (!requireNamespace("GUniFrac", quietly = TRUE)) {
      stop("The 'GMPR' method requires the 'GUniFrac' package.\n",
        "Please install it with:\n",
        "BiocManager::install('GUniFrac')",
        call. = FALSE
      )
    }

    size.factor <- GUniFrac::GMPR(t(Data))
    size.factor[is.na(size.factor)] <- mean(size.factor, na.rm = TRUE)
    Data <- Data / size.factor
  } else {
    stop("Invalid method. Please choose either 'TSS' or 'GMPR'.")
  }

  CrcBiomeScreenObject@NormalizedData <- as.data.frame(t(Data))
  attr(CrcBiomeScreenObject@NormalizedData, "NormalizationMethod") <- method
  attr(CrcBiomeScreenObject@NormalizedData, "Timestamp") <- Sys.time()

  return(CrcBiomeScreenObject)
}
