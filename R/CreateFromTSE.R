#' @title Create a CrcBiomeScreen object from TreeSummarizedExperiment
#' @description
#' Wrapper function to directly convert a TreeSummarizedExperiment object
#' into a A \code{CrcBiomeScreen} object. S4 object for downstream analysis.
#'
#' @param tse A TreeSummarizedExperiment object containing microbiome data.
#' @param assay_name Which assay to use (default: "counts" or "relative_abundance").
#' @return A A \code{CrcBiomeScreen} object. object.
#' @importFrom SummarizedExperiment assay rowData colData assayNames
#' @export
#'
#' @examples
#' # Runnable example using a minimal mock TreeSummarizedExperiment
#'
#' # Load required classes (SummarizedExperiment & TreeSummarizedExperiment)
#' suppressMessages({
#'   library(SummarizedExperiment)
#'   library(TreeSummarizedExperiment)
#' })
#'
#' # Create a tiny assay matrix
#' assay_mat <- matrix(
#'   c(10, 5,
#'     20, 7),
#'   nrow = 2,
#'   dimnames = list(c("Taxa1", "Taxa2"), c("S1", "S2"))
#' )
#'
#' # Create row (taxa) metadata
#' row_meta <- DataFrame(
#'   Taxa = c("A;B;C", "A;D;E")
#' )
#'
#' # Create sample metadata
#' col_meta <- DataFrame(
#'   number_reads = c(10000, 12000),
#'   condition = c("control", "CRC")
#' )
#'
#' # Build a minimal TreeSummarizedExperiment
#' tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
#'   assays = list(relative_abundance = assay_mat),
#'   rowData = row_meta,
#'   colData = col_meta
#' )
#'
#' # Convert to CrcBiomeScreen object
#' obj <- CreateCrcBiomeScreenObjectFromTSE(tse)
#'
#' # Inspect object
#' obj

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
