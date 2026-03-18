#' @title Create a CrcBiomeScreen S4 object for microbiome-based CRC analysis
#' @description
#' Constructor for the \code{CrcBiomeScreen} S4 class.
#' This function creates a structured container for microbiome data,
#' including absolute and relative abundance matrices, taxonomic annotations,
#' and sample metadata. It ensures compatibility with downstream modelling
#' and evaluation functions within the CrcBiomeScreen package.
#'
#' @param AbsoluteAbundance
#' A numeric matrix or data frame containing absolute abundance data.
#' @param TaxaData
#' A data frame containing taxonomic information for each feature.
#' @param SampleData
#' A data frame containing sample-level metadata.
#' @param RelativeAbundance
#' A numeric matrix or data frame containing relative abundance data.
#'
#' @details
#' If only relative abundance data are supplied, absolute abundance is estimated
#' using the total number of reads in \code{SampleData$number_reads}.
#'
#' @return A \linkS4class{CrcBiomeScreen} object.
#' \itemize{
#'   \item \code{AbsoluteAbundance}: Absolute abundance data.
#'   \item \code{RelativeAbundance}: Relative abundance data.
#'   \item \code{TaxaData}: Taxonomic annotations.
#'   \item \code{SampleData}: Sample metadata.
#'   \item \code{TaxaLevelData}: Optional genus-level summary data.
#'   \item \code{NormalizedData}: Normalized data.
#'   \item \code{OrginalNormalizedData}: Original normalized data.
#'   \item \code{ValidationData}: Optional validation dataset.
#'   \item \code{OutlierSamples}: Character vector of outlier sample names.
#'   \item \code{ModelData}, \code{ModelResult}, \code{EvaluateResult},
#'   \code{PredictResult}: Optional model results and evaluation outputs
#' }
#'
#' @seealso \linkS4class{CrcBiomeScreen}
#'
#' @importFrom rlang sym
#' @importFrom methods new validObject
#' @importFrom dplyr mutate across
#' @importFrom tibble tibble
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics abline barplot
#' @importFrom utils getFromNamespace
#' @importFrom ggplot2 aes geom_point geom_text labs theme theme_minimal ggsave
#' @export
#'
#' @examples
#' # Minimal example with tiny toy data (required for Bioconductor checks)
#'
#' # Create toy abundance matrices
#' rel_abund <- data.frame(
#'   Sample1 = c(10, 20, 70),
#'   Sample2 = c(30, 30, 40)
#' )
#' rownames(rel_abund) <- c("TaxaA", "TaxaB", "TaxaC")
#'
#' taxa_info <- data.frame(
#'   Taxa = rownames(rel_abund),
#'   stringsAsFactors = FALSE
#' )
#'
#' sample_info <- data.frame(
#'   number_reads = c(10000, 12000),
#'   condition = c("control", "CRC"),
#'   row.names = c("Sample1", "Sample2"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Create object
#' obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = rel_abund,
#'   TaxaData = taxa_info,
#'   SampleData = sample_info
#' )
#'

CreateCrcBiomeScreenObject <- function(
    AbsoluteAbundance = NULL,
    TaxaData = NULL,
    SampleData = NULL,
    RelativeAbundance = NULL
) {
  # --- Optional: handle direct TreeSummarizedExperiment input ------------
  if (inherits(AbsoluteAbundance, "TreeSummarizedExperiment")) {
    return(CreateCrcBiomeScreenObjectFromTSE(AbsoluteAbundance))
  }

  # --- Check inputs ------------------------------------------------------
  if (!is.null(RelativeAbundance) && is.null(AbsoluteAbundance)) {
    if (is.null(SampleData)) {
      stop("SampleData is required to
           convert RelativeAbundance to AbsoluteAbundance.")
    }
    if (!"number_reads" %in% colnames(SampleData)) {
      stop("SampleData must contain 'number_reads' to convert RelativeAbundance to AbsoluteAbundance.")
    }

    ra_df <- RelativeAbundance %>%
      t() %>%
      data.frame()

    ra_df <- ra_df %>%
      mutate(
        across(
          dplyr::everything(),
          ~ . * SampleData$number_reads / 100
        )
      )

    AbsoluteAbundance <- ra_df %>%
      t() %>%
      data.frame()
  }

  # --- Construct the S4 object -------------------------------------------
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
<<<<<<< Updated upstream
#' @alias CrcBiomeScreenObject-class
=======
#' @aliases CrcBiomeScreenObject-class
>>>>>>> Stashed changes
#' @description
#' An S4 container for CRC microbiome screening data, including abundance
#' matrices, taxonomy, sample metadata, and model results.
#'
#' @slot AbsoluteAbundance Absolute abundance matrix.
#' @slot TaxaData Taxonomy annotation data frame.
#' @slot SampleData Sample metadata (must include number_reads if relative
#' abundance is used).
#' @slot RelativeAbundance Relative abundance matrix.
#' @slot TaxaLevelData Optional genus-level summary.
#' @slot NormalizedData Normalized abundance data.
#' @slot ValidationData Optional validation dataset.
#' @slot ModelData Processed training/testing data.
#' @slot ModelResult Fitted model objects.
#' @slot EvaluateResult List of evaluation metrics (RF, XGBoost, etc.).
#' @slot PredictResult Predictions for external data.
#'
#' @return A \linkS4class{CrcBiomeScreen} object.
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
    PredictResult = "ANY",
    OutlierSamples = "character"
  )
)

#' @title Accessor for AbsoluteAbundance slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A \linkS4class{CrcBiomeScreen} object with A data.frame containing
#' absolute abundance data.
#'
#' @examples
#' # Construct minimal example object
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   AbsoluteAbundance = data.frame(TaxaA = c(1000, 2000)),
#'   RelativeAbundance = data.frame(TaxaA = c(10, 20)),
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = data.frame(
#'     number_reads = 10000,
#'     condition = "control"
#'   )
#' )
#'
#' # Retrieve absolute abundance
#' getAbsoluteAbundance(toy_obj)

#' @export
setGeneric("getAbsoluteAbundance",
           function(object) standardGeneric("getAbsoluteAbundance"))
#' @describeIn getAbsoluteAbundance Retrieve absolute abundance data
#' from a CrcBiomeScreen object.
setMethod("getAbsoluteAbundance", "CrcBiomeScreen",
          function(object) object@AbsoluteAbundance)

#' @title Accessor for RelativeAbundance slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A \linkS4class{CrcBiomeScreen} object with
#' a data.frame containing relative abundance data.
#' @examples
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = data.frame(TaxaA = c(10, 20)),
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = data.frame(
#'     number_reads = 10000,
#'     condition = "control"
#'   )
#' )
#'
#' getRelativeAbundance(toy_obj)
#' @export
setGeneric("getRelativeAbundance",
           function(object) standardGeneric("getRelativeAbundance"))
#' @describeIn getRelativeAbundance Retrieve relative abundance data
#' from a CrcBiomeScreen object.
setMethod("getRelativeAbundance", "CrcBiomeScreen",
          function(object) object@RelativeAbundance)

#' @title Accessor for SampleData slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A \linkS4class{CrcBiomeScreen} object with
#' a data.frame containing sample metadata.
#' @examples
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = data.frame(TaxaA = c(10, 20)),
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = data.frame(
#'     number_reads = 10000,
#'     condition = "control"
#'   )
#' )
#'
#' getSampleData(toy_obj)
#' @export
setGeneric("getSampleData", function(object) standardGeneric("getSampleData"))
#' @describeIn getSampleData Retrieve sample metadata
#' from a CrcBiomeScreen object.
setMethod("getSampleData", "CrcBiomeScreen", function(object) object@SampleData)

#' @title Accessor for ModelData slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A \linkS4class{CrcBiomeScreen} object with
#' a data.frame containing model data.
#' @examples
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = data.frame(TaxaA = c(10, 20)),
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = data.frame(
#'     number_reads = c(10000, 10000),
#'     condition = c("control", "CRC")
#'   )
#' )
#' toy_obj <- SplitDataSet(
#'   toy_obj,
#'   label = c("control", "CRC"),
#'   partition = 0.7
#' )
#' getModelData(toy_obj)
#'
#' @export
setGeneric("getModelData", function(object) standardGeneric("getModelData"))
#' @describeIn getSampleData Retrieve sample metadata
#' from a CrcBiomeScreen object.
setMethod("getModelData", "CrcBiomeScreen", function(object) object@ModelData)

#' @title Accessor for TaxaData slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A \linkS4class{CrcBiomeScreen} object with
#' a  data.frame containing taxonomic annotations.
#' @examples
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = data.frame(TaxaA = c(10, 20)),
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = data.frame(
#'     number_reads = 10000,
#'     condition = "control"
#'   )
#' )
#'
#' getTaxaData(toy_obj)
#' @export
setGeneric("getTaxaData", function(object) standardGeneric("getTaxaData"))
#' @describeIn getTaxaData Retrieve taxonomic annotations
#' from a CrcBiomeScreen object.
setMethod("getTaxaData", "CrcBiomeScreen", function(object) object@TaxaData)

#' @title Accessor for ModelResult slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A \linkS4class{CrcBiomeScreen} object with
#' a list containing fitted model results.
#' @examples
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = data.frame(TaxaA = c(10, 20)),
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = data.frame(
#'     number_reads = 10000,
#'     condition = "control"
#'   )
#' )
#'
#' getModelData(toy_obj)
#' @export
setGeneric("getModelResult", function(object) standardGeneric("getModelResult"))
#' @describeIn getModelResult Retrieve model results
#' from a CrcBiomeScreen object.
setMethod("getModelResult", "CrcBiomeScreen",
          function(object) object@ModelResult)

#' @title Accessor for OutlierSamples slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A \linkS4class{CrcBiomeScreen} object with
#' a list containing OutlierSamples results.
#' @examples
#' toy_obj <- CreateCrcBiomeScreenObject(
#' RelativeAbundance = data.frame(
#'   S1 = 10,
#'  S2 = 20,
#'   row.names = "TaxaA"
#' ),
#' TaxaData = data.frame(
#'   Taxa = "TaxaA",
#'   stringsAsFactors = FALSE
#' ),
#' SampleData = data.frame(
#'   number_reads = c(10000, 2000),
#'   condition = c("control", "control"),
#'   row.names = c("S1", "S2"),
#'   stringsAsFactors = FALSE
#' ))
#'
#' toy_obj@OutlierSamples <- c("S1")
#' getOutlierSamples(toy_obj)
#' @export
setGeneric("getOutlierSamples",
           function(object) standardGeneric("getOutlierSamples"))
#' @describeIn getOutlierSamples results
#' from a CrcBiomeScreen object.
setMethod("getOutlierSamples", "CrcBiomeScreen",
          function(object) object@OutlierSamples)


#' @title Accessor for PredictResult slot of CrcBiomeScreen object
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @return A \linkS4class{CrcBiomeScreen} object with
#' a  list containing fitted Prediction results.
#' @examples
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = data.frame(TaxaA = c(10, 20)),
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = data.frame(
#'     number_reads = 10000,
#'     condition = "control"
#'   )
#' )
#'
#' getPredictResult(toy_obj)
#' @export
setGeneric("getPredictResult",
           function(object) standardGeneric("getPredictResult"))
#' @describeIn getPredictResult Prediction results
#' from a CrcBiomeScreen object.
setMethod("getPredictResult", "CrcBiomeScreen",
          function(object) object@PredictResult)

#' @title setTaxaData<-: Setter for TaxaData slot of CrcBiomeScreen object
#' @name setTaxaData-setter
#' @aliases setTaxaData<-
#' @usage setTaxaData(object) <- value
#'
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @param value A data.frame containing updated taxonomic annotations.
#'
#' @return A \linkS4class{CrcBiomeScreen} object with
#' a modified \linkS4class{CrcBiomeScreen} object.
#' @examples
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = data.frame(TaxaA = c(10, 20)),
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = data.frame(
#'     number_reads = 10000,
#'     condition = "control"
#'   )
#' )
#'
#' setTaxaData(toy_obj) <- data.frame(Taxa = "NewTaxa")
#' getTaxaData(toy_obj)
#' @export
setGeneric("setTaxaData<-",
           function(object, value) standardGeneric("setTaxaData<-"))

#' @describeIn setTaxaData-setter Replace the TaxaData slot
#' of a CrcBiomeScreen object.
setReplaceMethod("setTaxaData", "CrcBiomeScreen", function(object, value) {
  object@TaxaData <- value
  validObject(object)
  return(object)
})

#' @title Accessor for NormalizedData slot of CrcBiomeScreen object
#' @name getNormalizedData
#'
<<<<<<< Updated upstream
#' @description Retrieve normalized abundance data
#' from a \linkS4class{CrcBiomeScreen} object.
=======
#' @description Retrieve or set normalized abundance data
#' from/in a \linkS4class{CrcBiomeScreen} object.
>>>>>>> Stashed changes
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @param value A data.frame (or matrix) containing normalized abundance data.
#'
#' @return A \linkS4class{CrcBiomeScreen} object with
<<<<<<< Updated upstream
#' a  data.frame (or matrix) containing normalized abundance data.
=======
#' a data.frame (or matrix) containing normalized abundance data.
>>>>>>> Stashed changes
#' @examples
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = data.frame(TaxaA = c(10, 20)),
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = data.frame(
#'     number_reads = 10000,
#'     condition = "control"
#'   )
#' )
#'
#' getNormalizedData(toy_obj)
#' @export
setGeneric("getNormalizedData",
           function(object) standardGeneric("getNormalizedData"))
<<<<<<< Updated upstream
=======

#' @rdname getNormalizedData
#' @export
setGeneric("getNormalizedData<-",
           function(object, value) standardGeneric("getNormalizedData<-"))
>>>>>>> Stashed changes

#' @describeIn getNormalizedData Retrieve normalized abundance data.
setMethod("getNormalizedData", "CrcBiomeScreen",
          function(object) object@NormalizedData)

<<<<<<< Updated upstream
=======
#' @describeIn getNormalizedData Set normalized abundance data.
setReplaceMethod("getNormalizedData", "CrcBiomeScreen",
                 function(object, value) {
                   object@NormalizedData <- value
                   validObject(object)
                   return(object)
                 })

>>>>>>> Stashed changes
#' @title setNormalizedData<-: Setter for NormalizedData slot
#' of CrcBiomeScreen object
#' @name setNormalizedData-setter
#' @aliases setNormalizedData<-
#' @usage setNormalizedData(object) <- value
#'
#' @param object A \linkS4class{CrcBiomeScreen} object.
#' @param value A data.frame or matrix containing normalized abundance data.
#'
#' @return A \linkS4class{CrcBiomeScreen} object with
#' a  modified \linkS4class{CrcBiomeScreen} object.
#' @examples
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = data.frame(TaxaA = c(10, 20)),
#'   TaxaData = data.frame(Taxa = "TaxaA"),
#'   SampleData = data.frame(
#'     number_reads = 10000,
#'     condition = "control"
#'   )
#' )
#'
#' setNormalizedData(toy_obj) <- data.frame(n1 = 1:2)
#' getNormalizedData(toy_obj)
#' @export
setGeneric("setNormalizedData<-",
           function(object, value) standardGeneric("setNormalizedData<-"))

#' @describeIn setNormalizedData-setter Replace the NormalizedData slot
#' of a CrcBiomeScreen object.
setReplaceMethod("setNormalizedData", "CrcBiomeScreen"
                 , function(object, value) {
  object@NormalizedData <- value
  validObject(object)
  return(object)
})
