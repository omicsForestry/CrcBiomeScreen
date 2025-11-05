#' @title Create a CrcBiomeScreen object from TreeSummarizedExperiment
#' @description
#' Wrapper function to directly convert a TreeSummarizedExperiment object
#' into a \linkS4class{CrcBiomeScreen} S4 object for downstream analysis.
#'
#' @param tse A TreeSummarizedExperiment object containing microbiome data.
#' @param assay_name Which assay to use (default: "counts" or "relative_abundance").
#' @return A \linkS4class{CrcBiomeScreen} object.
#' @importFrom SummarizedExperiment assay rowData colData assayNames
#' @export
#'
#' @examples
#' if (requireNamespace("curatedMetagenomicData", quietly = TRUE)) {
#'   tse <- curatedMetagenomicData::curatedMetagenomicData(
#'     "ThomasAM_2018a.relative_abundance",
#'     dryrun = FALSE, rownames = "short"
#'   )[[1]]
#'   obj <- CreateCrcBiomeScreenObjectFromTSE(tse)
#'   obj
#' }
CreateCrcBiomeScreenObjectFromTSE <- function(tse, assay_name = NULL) {
  if (!inherits(tse, "TreeSummarizedExperiment")) {
    stop("Input must be a TreeSummarizedExperiment object.")
  }

  if (is.null(assay_name)) {
    assay_name <- if ("relative_abundance" %in% SummarizedExperiment::assayNames(tse)) {
      "relative_abundance"
    } else {
      SummarizedExperiment::assayNames(tse)[1]
    }
  }

  rel_abund <- SummarizedExperiment::assay(tse, assay_name)
  taxa_data <- as.data.frame(SummarizedExperiment::rowData(tse))
  sample_data <- as.data.frame(SummarizedExperiment::colData(tse))

  CreateCrcBiomeScreenObject(
    RelativeAbundance = rel_abund,
    TaxaData = taxa_data,
    SampleData = sample_data
  )
}
