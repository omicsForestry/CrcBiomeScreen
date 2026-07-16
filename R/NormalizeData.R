#' Normalise the absolute data to relative data by using Total Sum Scaling
#'  and Geometric Mean of Pairwise Ratios (GMPR)
#'
#' @param CrcBiomeScreenObject From the CreateCrcBiomeScreenObject()
#' @param method "TSS" or "GMPR"
#' @param level Taxonomic level for normalization, e.g., "Genus"
#'

#' @return A A \code{CrcBiomeScreen} object. with the updated NormalizedData.
#' @export
#'
#' @examples
#' # Minimal runnable example for NormalizeData
#'
#' # Toy taxa in a simplified MetaPhlAn-like hierarchical format
#' toy_taxa <- data.frame(
#'   Taxa = c(
#'     "D_0__Bacteria|D_1__Firmicutes|D_2__Clostridia|D_3__OrderX|D_4__FamilyX|D_5__GenusA",
#'     "D_0__Bacteria|D_1__Firmicutes|D_2__Clostridia|D_3__OrderY|D_4__FamilyY|D_5__GenusB"
#'   ),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Toy abundance matrix (2 taxa, 2 samples)
#' toy_abs <- data.frame(
#'   S1 = c(10, 5),
#'   S2 = c(20, 15)
#' )
#' rownames(toy_abs) <- toy_taxa$Taxa
#'
#' # Dummy sample metadata
#' toy_sample <- data.frame(
#'   sample_id = c("S1", "S2")
#' )
#'
#' # Construct minimal CrcBiomeScreen object
#' toy_obj <- new(
#'   "CrcBiomeScreen",
#'   AbsoluteAbundance   = toy_abs,
#'   RelativeAbundance   = data.frame(),
#'   TaxaData            = toy_taxa,
#'   SampleData          = toy_sample,
#'   TaxaLevelData       = NULL,
#'   NormalizedData      = NULL,
#'   OriginalNormalizedData = NULL,
#'   ValidationData      = NULL,
#'   ModelData           = NULL,
#'   ModelResult         = NULL,
#'   EvaluateResult      = list(),
#'   PredictResult       = NULL
#' )
#'
#' # Apply taxonomy splitting + keep genus level
#' toy_obj <- SplitTaxas(toy_obj)
#' toy_obj <- KeepTaxonomicLevel(toy_obj, level = "Genus")
#' toy_obj <- NormalizeData(toy_obj, method = "TSS", level = "Genus")
#' # Inspect normalized results
#' head(getNormalizedData(toy_obj))
#'
NormalizeData <- function(CrcBiomeScreenObject = NULL, method = NULL, level = NULL) {
  Data <- t(CrcBiomeScreenObject@TaxaLevelData[[paste0(level, "LevelData")]])

  # Data format used internally:
  # rows = samples, columns = microbial features
  n_samples <- nrow(Data)
  n_features <- ncol(Data)

  if (method == "TSS") {
    # Calculate the total number of counts in each sample
    total_counts <- rowSums(Data)

    if (any(total_counts == 0)) {
      warning(
        "Some samples have zero total counts. These samples may produce NA values ",
        "during TSS normalisation.",
        call. = FALSE
      )
    }

    # Convert absolute abundance to relative abundance
    Data <- Data / total_counts

  } else if (method == "GMPR") {
    if (!requireNamespace("GUniFrac", quietly = TRUE)) {
      stop(
        "The 'GMPR' method requires the 'GUniFrac' package.\n",
        "Please install it with:\n",
        "BiocManager::install('GUniFrac')",
        call. = FALSE
      )
    }

    # GMPR expects rows = features and columns = samples.
    # Internally, Data has rows = samples and columns = features,
    # so we pass t(Data) to GUniFrac::GMPR().
    if (n_samples > n_features) {
      warning(
        "The number of samples is larger than the number of microbial features ",
        "after keeping the selected taxonomic level. GUniFrac::GMPR() may issue ",
        "a warning asking whether samples are arranged in columns. In this function, ",
        "samples are stored in rows internally and are transposed before GMPR is run. ",
        "If your CrcBiomeScreenObject was created correctly, this warning can usually ",
        "be ignored.",
        call. = FALSE
      )
    }

    gmpr_input <- t(Data)

    size.factor <- withCallingHandlers(
      GUniFrac::GMPR(gmpr_input),
      warning = function(w) {
        if (grepl("Sample size is larger than OTU number", conditionMessage(w))) {
          warning(
            "GUniFrac::GMPR() reported that the sample size is larger than the OTU ",
            "number. This can occur in datasets where the number of samples is genuinely ",
            "larger than the number of features at the selected taxonomic level. ",
            "The input has already been transposed so that rows = features and columns = samples.",
            call. = FALSE
          )
          invokeRestart("muffleWarning")
        }
      }
    )

    size.factor[is.na(size.factor)] <- mean(size.factor, na.rm = TRUE)

    Data <- Data / size.factor

  } else {
    stop("Invalid method. Please choose either 'TSS' or 'GMPR'.")
  }

  CrcBiomeScreenObject@NormalizedData <- as.data.frame(Data)
  attr(CrcBiomeScreenObject@NormalizedData, "NormalizationMethod") <- method
  attr(CrcBiomeScreenObject@NormalizedData, "Timestamp") <- Sys.time()
  attr(CrcBiomeScreenObject@NormalizedData, "TaxonomicLevel") <- level
  attr(CrcBiomeScreenObject@NormalizedData, "Samples") <- n_samples
  attr(CrcBiomeScreenObject@NormalizedData, "Features") <- n_features

  return(CrcBiomeScreenObject)
}
