#' @titile Create the CrcBiomeScreenObject for analysis
#' @description This function creates the CrcBiomeScreenObject to contain the data for analysis.
#' @param AbsoluteAbundance A numeric matrix or data frame containing absolute abundance data
#' @param TaxaData A data frame containing taxonomic information for each feature
#' @param SampleData A data frame containing sample metadata
#' @param RelativeAbundance A numeric matrix or data frame containing relative abundance data
#' 
#' @importFrom magrittr %>%
#' @import dplyr across mutate
#' 
#' @return An object with data
#' @export
#'
#' @examples
#' toydata <- curatedMetagenomicData("ThomasAM_2018a.relative_abundance"
#'                                  , dryrun = FALSE, rownames = "short")
#'
#'CrcBiomeScreenObject <- CreateCrcBiomeScreenObject(RelativeAbundance = toydata[[1]]@assays@data@listData$relative_abundance,
#'                                                  TaxaData = toydata[[1]]@rowLinks$nodeLab,
#'                                                  SampleData = toydata[[1]]@colData)

CreateCrcBiomeScreenObject <- function(AbsoluteAbundance = NULL, TaxaData = NULL, SampleData = NULL, RelativeAbundance = NULL) {
  # If AbsoluteAbundance is NULL, check if RelativeAbundance is provided
  if (!is.null(RelativeAbundance) && is.null(AbsoluteAbundance)) {
    if (is.null(SampleData)) {
      stop("SampleData is required to convert RelativeAbundance to AbsoluteAbundance.")
    }
    if (!"number_reads" %in% colnames(SampleData)) {
      stop("SampleData must contain a 'number_reads' column to convert RelativeAbundance to AbsoluteAbundance.")
    }
      AbsoluteAbundance <- RelativeAbundance %>%
      t() %>%
      data.frame() %>%
      mutate(across(seq_len(dim(RelativeAbundance)[2]), ~ (. * SampleData$number_reads / 100))) %>%
      t() %>%
      data.frame()
  }

  # Set up the object
  obj <- list(
    AbsoluteAbundance = AbsoluteAbundance,
    TaxaData = TaxaData,
    SampleData = SampleData,
    RelativeAbundance = RelativeAbundance,
    GenusLevelData = NULL,
    ValidationData = NULL,
    ModelData = NULL,
    ModelResult = NULL,
    EvaluateResult = list(
      RF = NULL,
      XGBoost = NULL
    ),
    PredictResult = NULL
  )

  # Set up the name of the class
  class(obj) <- "CrcBiomeScreenObject"
  return(obj)
}
