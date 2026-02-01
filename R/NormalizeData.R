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
#' # Minimal runnable example for NormalizeData
#'
#' # Create toy absolute abundance matrix (2 samples × 3 taxa)
#' abs_abund <- data.frame(
#'   S1 = c(10, 30, 60),
#'   S2 = c(20, 20, 60)
#' )
#' rownames(abs_abund) <- c("TaxaA", "TaxaB", "TaxaC")
#'
#' # Create toy taxonomy with a valid Genus column
#' toy_taxa <- data.frame(
#'   Taxa = rownames(abs_abund),
#'   Genus = c("G1", "G2", "G3"),   # <-- REQUIRED so level="Genus" works
#'   stringsAsFactors = FALSE
#' )
#'
#' # Sample metadata
#' toy_samples <- data.frame(
#'   study_condition = c("control", "CRC"),
#'   row.names = c("S1", "S2")
#' )
#'
#' # Construct toy CrcBiomeScreen object
#' toy_obj <- new(
#'   "CrcBiomeScreen",
#'   AbsoluteAbundance   = abs_abund,
#'   RelativeAbundance   = data.frame(),
#'   TaxaData            = toy_taxa,
#'   SampleData          = toy_samples,
#'   TaxaLevelData       = NULL,
#'   NormalizedData      = NULL,
#'   OrginalNormalizedData = NULL,
#'   ValidationData      = NULL,
#'   ModelData           = NULL,
#'   ModelResult         = NULL,
#'   EvaluateResult      = list(),
#'   PredictResult       = NULL
#' )
#'
#' # Apply TSS normalization (always runnable)
#' \donttest{norm_obj <- NormalizeData(toy_obj, method = "TSS", level = "Genus")
#'
#' # Inspect normalized results
#' head(norm_obj@NormalizedData)}


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
