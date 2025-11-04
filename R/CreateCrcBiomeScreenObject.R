#' @title Create a CrcBiomeScreen S4 object for microbiome-based CRC analysis
#' @description
#' Constructor for the \code{CrcBiomeScreen} S4 class.
#' This function creates a structured container for microbiome data,
#' including absolute and relative abundance matrices, taxonomic annotations,
#' and sample metadata. It ensures compatibility with downstream modelling
#' and evaluation functions within the CrcBiomeScreen package.
#'
#' @param AbsoluteAbundance A numeric matrix or data frame containing absolute abundance data.
#' @param TaxaData A data frame containing taxonomic information for each feature.
#' @param SampleData A data frame containing sample-level metadata.
#' @param RelativeAbundance A numeric matrix or data frame containing relative abundance data.
#'
#' @details
#' If only relative abundance data are supplied, absolute abundance is estimated
#' using the total number of reads in \code{SampleData$number_reads}.
#'
#' @return
#' A \linkS4class{CrcBiomeScreen} object containing the following slots:
#' \itemize{
#'   \item \code{AbsoluteAbundance}: Absolute abundance data.
#'   \item \code{RelativeAbundance}: Relative abundance data.
#'   \item \code{TaxaData}: Taxonomic annotations.
#'   \item \code{SampleData}: Sample metadata.
#'   \item \code{TaxaLevelData}: Optional genus-level summary data.
#'   \item \code{NormalizedData}: Normalized data.
#'   \item \code{OrginalNormalizedData}: Original normalized data.
#'   \item \code{ValidationData}: Optional validation dataset.
#'   \item \code{ModelData}, \code{ModelResult}, \code{EvaluateResult}, \code{PredictResult}: Optional model results and evaluation outputs.
#' }
#'
#' @seealso \linkS4class{CrcBiomeScreen}
#'
#' @importFrom dplyr mutate across
#' @importFrom tibble tibble
#' @export
#'
#' @examples
#' if (requireNamespace("curatedMetagenomicData", quietly = TRUE)) {
#'   toydata <- curatedMetagenomicData::curatedMetagenomicData(
#'     "ThomasAM_2018a.relative_abundance",
#'     dryrun = FALSE, rownames = "short"
#'   )
#'   CrcBiomeScreenObject <- CreateCrcBiomeScreenObject(
#'     RelativeAbundance = toydata[[1]]@assays@data@listData$relative_abundance,
#'     TaxaData = toydata[[1]]@rowLinks$nodeLab,
#'     SampleData = toydata[[1]]@colData
#'   )
#' }

CreateCrcBiomeScreenObject <- function(
    AbsoluteAbundance = NULL,
    TaxaData = NULL,
    SampleData = NULL,
    RelativeAbundance = NULL
) {
  # Check inputs -----------------------------------------------------------
  if (is.null(RelativeAbundance) && is.null(AbsoluteAbundance)) {
    stop("Either RelativeAbundance or AbsoluteAbundance must be provided.")
  }

  if (!is.null(RelativeAbundance) && is.null(AbsoluteAbundance)) {
    if (is.null(SampleData)) {
      stop("SampleData is required to convert RelativeAbundance to AbsoluteAbundance.")
    }
    if (!"number_reads" %in% colnames(SampleData)) {
      stop("SampleData must contain 'number_reads' to convert RelativeAbundance to AbsoluteAbundance.")
    }

    AbsoluteAbundance <- RelativeAbundance %>%
      t() %>%
      data.frame() %>%
      mutate(across(
        seq_len(dim(RelativeAbundance)[2]),
        ~ (. * SampleData$number_reads / 100)
      )) %>%
      t() %>%
      data.frame()
  }

  # Construct the object ---------------------------------------------------
  new("CrcBiomeScreen",
      AbsoluteAbundance = as.data.frame(AbsoluteAbundance),
      TaxaData = if (is.null(TaxaData)) data.frame() else as.data.frame(TaxaData),
      SampleData = if (is.null(SampleData)) data.frame() else as.data.frame(SampleData),
      RelativeAbundance = if (is.null(RelativeAbundance)) data.frame() else as.data.frame(RelativeAbundance),
      TaxaLevelData = NULL,
      NormalizedData = NULL,
      OrginalNormalizedData = NULL,
      ValidationData = NULL,
      ModelData = NULL,
      ModelResult = NULL,
      EvaluateResult = list(RF = NULL, XGBoost = NULL),
      PredictResult = NULL)
}


# Accessor methods ---------------------------------------------------------
#' @title CrcBiomeScreen Class
#' @description
#' An S4 container for CRC microbiome screening data, including abundance
#' matrices, taxonomy, sample metadata, and model results.
#'
#' @slot AbsoluteAbundance Absolute abundance matrix.
#' @slot TaxaData Taxonomy annotation data frame.
#' @slot SampleData Sample metadata (must include number_reads if relative abundance is used).
#' @slot RelativeAbundance Relative abundance matrix.
#' @slot TaxaLevelData Optional genus-level summary.
#' @slot NormalizedData Normalized abundance data.
#' @slot ValidationData Optional validation dataset.
#' @slot ModelData Processed training/testing data.
#' @slot ModelResult Fitted model objects.
#' @slot EvaluateResult List of evaluation metrics (RF, XGBoost, etc.).
#' @slot PredictResult Predictions for external data.
#'
#'
setClass(
  "CrcBiomeScreen",
  slots = c(
    AbsoluteAbundance = "data.frame",
    TaxaData = "data.frame",
    SampleData = "data.frame",
    RelativeAbundance = "data.frame",
    TaxaLevelData = "ANY",
    NormalizedData = "ANY",
    OrginalNormalizedData = "ANY",
    ValidationData = "ANY",
    ModelData = "ANY",
    ModelResult = "ANY",
    EvaluateResult = "list",
    PredictResult = "ANY"
  )
)

#' @title Accessor for AbsoluteAbundance slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A data.frame containing absolute abundance data.
#' @export
setGeneric("getAbsoluteAbundance", function(object) standardGeneric("getAbsoluteAbundance"))
#' @describeIn getAbsoluteAbundance Retrieve absolute abundance data from a CrcBiomeScreen object.
setMethod("getAbsoluteAbundance", "CrcBiomeScreen", function(object) object@AbsoluteAbundance)

#' @title Accessor for RelativeAbundance slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A data.frame containing relative abundance data.
#' @export
setGeneric("getRelativeAbundance", function(object) standardGeneric("getRelativeAbundance"))
#' @describeIn getRelativeAbundance Retrieve relative abundance data from a CrcBiomeScreen object.
setMethod("getRelativeAbundance", "CrcBiomeScreen", function(object) object@RelativeAbundance)

#' @title Accessor for SampleData slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A data.frame containing sample metadata.
#' @export
setGeneric("getSampleData", function(object) standardGeneric("getSampleData"))
#' @describeIn getSampleData Retrieve sample metadata from a CrcBiomeScreen object.
setMethod("getSampleData", "CrcBiomeScreen", function(object) object@SampleData)

#' @title Accessor for TaxaData slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A data.frame containing taxonomic annotations.
#' @export
setGeneric("getTaxaData", function(object) standardGeneric("getTaxaData"))
#' @describeIn getTaxaData Retrieve taxonomic annotations from a CrcBiomeScreen object.
setMethod("getTaxaData", "CrcBiomeScreen", function(object) object@TaxaData)

#' @title Accessor for ModelResult slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A list containing fitted model results.
#' @export
setGeneric("getModelResult", function(object) standardGeneric("getModelResult"))
#' @describeIn getModelResult Retrieve model results from a CrcBiomeScreen object.
setMethod("getModelResult", "CrcBiomeScreen", function(object) object@ModelResult)
