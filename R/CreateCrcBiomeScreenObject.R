#' Create the CrcBiomeScreenObject for analysis
#' This function creates the CrcBiomeScreenObject to contain the data for analysis.
#'
#' @param AbsoluteAbundance A numeric matrix or data frame containing absolute abundance data
#' @param TaxaData A data frame containing taxonomic information for each feature
#' @param SampleData A data frame containing sample metadata
#' @param RelativeAbundance A numeric matrix or data frame containing relative abundance data
#'
#' @return An object with data
#' @export
#'
#' @details
#' The `CrcBiomeScreenObject` is a list-based S3 object designed to store all data, models, and results throughout the analysis workflow. This ensures consistency and makes it easy to track changes. The object is initialized with several key slots, and other slots are populated by subsequent functions.
#'
#' \describe{
#'   \item{\strong{AbsoluteAbundance}}{A data frame containing absolute abundance data (counts).}
#'   \item{\strong{TaxaData}}{A data frame mapping features to their taxonomic information.}
#'   \item{\strong{SampleData}}{A data frame containing sample metadata, such as clinical status or the number of sequencing reads.}
#'   \item{\strong{RelativeAbundance}}{A data frame containing relative abundance data.}
#'   \item{\strong{TaxaLevelData}}{A list of data frames containing abundance data aggregated to a specific taxonomic level (e.g., genus, family). Populated by `KeepTaxonomicLevel`.}
#'   \item{\strong{NormalizedData}}{A data frame containing the normalized abundance data. Populated by `NormalizeData`.}
#'   \item{\strong{ValidationData}}{A list containing the processed external validation data.}
#'   \item{\strong{ModelData}}{A list containing the processed training and test data, ready for model training. Populated by `SplitData`.}
#'   \item{\strong{ModelResult}}{A list storing the trained model objects (e.g., `randomForest` or `xgboost` objects). Populated by `TrainModel`.}
#'   \item{\strong{EvaluateResult}}{A list storing performance metrics (e.g., AUC, ROC curve) for the models on the test set. Populated by `EvaluateModel`.}
#'   \item{\strong{PredictResult}}{A list storing prediction probabilities and evaluation metrics for external validation sets. Populated by `ValidateModelOnData`.}
#'   \item{\strong{TrainingPredictions}}{A list storing out-of-bag (OOB) or cross-validation predictions for the training data. Populated by `TrainModel`.}
#' }
#'
#' @examples
#' toydata <- curatedMetagenomicData("ThomasAM_2018a.relative_abundance"
#'                                 , dryrun = FALSE, rownames = "short")
#'
#' CrcBiomeScreenObject <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = toydata[[1]]@assays@data@listData$relative_abundance,
#'   TaxaData = toydata[[1]]@rowLinks$nodeLab,
#'   SampleData = toydata[[1]]@colData
#' )
#' 
CreateCrcBiomeScreenObject <- function(AbsoluteAbundance = NULL, TaxaData = NULL, SampleData = NULL, RelativeAbundance = NULL) {
  # If AbsoluteAbundance is NULL, check if RelativeAbundance is provided
  if (!is.null(RelativeAbundance) && is.null(AbsoluteAbundance)) {
    if (is.null(SampleData)) {
      stop("SampleData is required to convert RelativeAbundance to AbsoluteAbundance.")
    }
    if (!"number_reads" %in% colnames(SampleData)) {
      print("To convert RelativeAbundance to AbsoluteAbundance requires the total number of reads in each sample.")
    } else {
      AbsoluteAbundance <- RelativeAbundance %>%
        t() %>%
        data.frame() %>%
        mutate(across(seq_len(dim(RelativeAbundance)[2]), ~ (. * SampleData$number_reads / 100))) %>%
        t() %>%
        data.frame()
    }
  }
  
  # Set up the object
  obj <- list(
    AbsoluteAbundance = AbsoluteAbundance,
    TaxaData = TaxaData,
    SampleData = SampleData,
    RelativeAbundance = RelativeAbundance,
    TaxaLevelData = list( # To store abundance aggregated to different taxonomic levels
      PhylumLevelData = NULL,
      ClassLevelData = NULL,
      OrderLevelData = NULL,
      FamilyLevelData = NULL,
      GenusLevelData = NULL,
      SpeciesLevelData = NULL
    ),
    NormalizedData = NULL, # To store normalized abundance data
    ValidationData = NULL, # For external validation data
    ModelData = NULL, # To store processed data for modeling (train/test split)
    ModelResult = NULL, # To store trained model objects
    
    # Evaluation and prediction results
    EvaluateResult = list(
      RF = NULL,
      XGBoost = NULL
    ),
    PredictResult = NULL # To store predictions from new, external data
  )
  
  # Set up the name of the class
  class(obj) <- "CrcBiomeScreenObject"
  return(obj)
}
