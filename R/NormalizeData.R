#' Normalise the absolute data to relative data by using Total Sum Scaling
#'  and Geometric Mean of Pairwise Ratios (GMPR)
#'
#' @param CrcBiomeScreenObject From the CreateCrcBiomeScreenObject()
#' @param method "TSS" or "GMPR"
#' @param level Taxonomic level for normalization, e.g., "Genus"
#'
#'
#' @return A \linkS4class{CrcBiomeScreenObject} with the updated NormalizedData.
#' @export
#'
#' @examples
#' # NormalizeData expects the data to already be aggregated at a specific
#' # taxonomic level in the TaxaLevelData slot.
#' # We mock this pre-processed 'GenusLevelData' (rows = samples, cols = taxa).
#' mock_genus_data <- data.frame(
#'   G1 = c(10, 20),
#'   G2 = c(30, 20),
#'   G3 = c(60, 60),
#'   row.names = c("S1", "S2")
#' )
#'
#' # Construct toy CrcBiomeScreen object
#' toy_obj <- new(
#'   "CrcBiomeScreen",
#'   AbsoluteAbundance = data.frame(),
#'   RelativeAbundance = data.frame(),
#'   TaxaData = data.frame(),
#'   SampleData = data.frame(),
#'   TaxaLevelData = list(GenusLevelData = mock_genus_data), # <-- TUTAJ ZMIANA
#'   NormalizedData = data.frame(),
#'   OrginalNormalizedData = data.frame(),
#'   ValidationData = NULL,
#'   ModelData = NULL,
#'   ModelResult = NULL,
#'   EvaluateResult = list(),
#'   PredictResult = NULL
#' )
#'
#' # Apply TSS normalization
#' if (interactive()) {
#'   norm_obj <- NormalizeData(toy_obj, method = "TSS", level = "Genus")
#'   print(head(getNormalizedData(norm_obj)))
#' }
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
