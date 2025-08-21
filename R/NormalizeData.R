#' Normalise the absolute data to relative data by using Total Sum Scaling
#'  and Geometric Mean of Pairwise Ratios (GMPR)
#'
#' @param CrcBiomeScreenObject From the CreateCrcBiomeScreenObject()
#' @param method "TSS" or "GMPR"
#'
#' @return Updated CrcBiomeScreenObject with NormalizedData
#' @export
#'
#' @examples
#' \dontrun{
#' # Example only runs if suggested packages are installed
#' if (requireNamespace("GUniFrac", quietly = TRUE) &&
#'     requireNamespace("microbiomeMarker", quietly = TRUE)) {
#'   CrcBiomeScreenObject <- NormalizeData(CrcBiomeScreenObject,
#'                                         method = "GMPR")
#'   attributes(CrcBiomeScreenObject$NormalizedData)
#' }
#' }
NormalizeData <- function(CrcBiomeScreenObject = NULL, method = NULL) {
  Data <- CrcBiomeScreenObject$GenusLevelData

  if (method == "TSS") {
    if (!requireNamespace("microbiomeMarker", quietly = TRUE)) {
      stop("The 'TSS' method requires the 'microbiomeMarker' package.\n",
           "Please install it with:\n",
           "remotes::install_github('yiluheihei/microbiomeMarker')",
           call. = FALSE)
    }
    Data <- as.data.frame(
      t(microbiomeMarker::normalize(t(Data), method = "TSS"))
    )

  } else if (method == "GMPR") {
    if (!requireNamespace("GUniFrac", quietly = TRUE)) {
      stop("The 'GMPR' method requires the 'GUniFrac' package.\n",
           "Please install it with:\n",
           "BiocManager::install('GUniFrac')",
           call. = FALSE)
    }
    size.factor <- GUniFrac::GMPR(t(Data))
    size.factor[is.na(size.factor)] <- mean(size.factor, na.rm = TRUE)
    Data <- Data / size.factor

  } else {
    stop("Invalid method. Please choose either 'TSS' or 'GMPR'.")
  }

  CrcBiomeScreenObject$NormalizedData <- Data
  attr(CrcBiomeScreenObject$NormalizedData, "NormalizationMethod") <- method
  attr(CrcBiomeScreenObject$NormalizedData, "Timestamp") <- Sys.time()

  return(CrcBiomeScreenObject)
}
