

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
    NormalizedData = NULL,
    ModelData = NULL,
    ScreeningResult = NULL,
    ModelResult = NULL,
    PredictResult = NULL,
    Params = list()
  )
  
  # 设置类名
  class(obj) <- "CrcBiomeScreenObject"
  return(obj)
}
