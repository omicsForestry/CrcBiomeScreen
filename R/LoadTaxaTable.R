#' Load a custom taxa table for ASV/OTU data
#'
#' This function allows users to load their own taxonomic assignments for ASV/OTU data.
#' The input table should map sequence IDs to their full taxonomic lineage.
#'
#' @param CrcBiomeScreenObject The CrcBiomeScreenObject to which the taxa table will be added.
#' @param taxa_table A data frame. It must contain at least two columns: one for sequence IDs (e.g., ASV or OTU names)
#'                   and another for the corresponding taxonomic lineage string.
#' @param id_column The name of the column in `taxa_table` that contains the sequence IDs (ASV/OTU names).
#' @param taxa_column The name of the column in `taxa_table` that contains the taxonomic lineage string.
#'
#' @return The \linkS4class{CrcBiomeScreenObject} with the loaded taxa table.
#' @export
#' @examples
#' ## Minimal example using CreateCrcBiomeScreenObject and LoadTaxaTable
#'
#' # Toy relative abundance matrix: 1 taxon (row) × 2 samples (columns)
#' rel_abund <- data.frame(
#'   S1 = 10,
#'   S2 = 20,
#'   row.names = "TaxaA"
#' )
#'
#' # Sample metadata with required 'number_reads' column
#' sample_info <- data.frame(
#'   number_reads = c(10000, 12000),
#'   condition    = c("control", "CRC"),
#'   row.names    = c("S1", "S2"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Simple taxa table matching the row names of rel_abund
#' taxa_info <- data.frame(
#'   Taxa = rownames(rel_abund),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Construct a minimal CrcBiomeScreen object
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = rel_abund,
#'   TaxaData          = taxa_info,
#'   SampleData        = sample_info
#' )
#'
#' # External taxonomy table to be merged in by LoadTaxaTable
#' my_taxa_table <- data.frame(
#'   ASV_ID   = "TaxaA",
#'   Taxonomy = "D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia",
#'   stringsAsFactors = FALSE
#' )
#'
#' # Load taxonomy table into the CrcBiomeScreen object
#' toy_obj <- LoadTaxaTable(
#'   CrcBiomeScreenObject = toy_obj,
#'   taxa_table  = my_taxa_table,
#'   id_column   = "ASV_ID",
#'   taxa_column = "Taxonomy"
#' )
#'
#' # Inspect updated taxonomy using the accessor
#' head(getTaxaData(toy_obj))
#'
LoadTaxaTable <- function(CrcBiomeScreenObject,
                          taxa_table,
                          id_column,
                          taxa_column) {

  # --- Input Validation ---
  if (!inherits(CrcBiomeScreenObject, "CrcBiomeScreen"))
    stop("Input must be a CrcBiomeScreen object.")

  if (!is.data.frame(taxa_table))
    stop("'taxa_table' must be a data frame.")

  if (!id_column %in% colnames(taxa_table))
    stop(sprintf("Column '%s' not found in taxa_table.", id_column))

  if (!taxa_column %in% colnames(taxa_table))
    stop(sprintf("Column '%s' not found in taxa_table.", taxa_column))


  # --- Process new taxonomy table ---
  processed_taxa <- taxa_table %>%
    dplyr::select(
      SequenceID = !!rlang::sym(id_column),
      TaxaString = !!rlang::sym(taxa_column)
    )


  # --- Existing TaxaData ---
  old_taxa <- CrcBiomeScreenObject@TaxaData

  # Case 1: No old taxonomy → direct assign
  if (nrow(old_taxa) == 0) {
    CrcBiomeScreenObject@TaxaData <- processed_taxa
    return(CrcBiomeScreenObject)
  }

  # Identify potential ID columns in old_taxa
  possible_id_cols <- c("SequenceID", "Taxa", "ASV", "FeatureID")
  merge_col <- intersect(colnames(old_taxa), possible_id_cols)

  # Case 2: No common ID column → warn and replace
  if (length(merge_col) == 0) {
    warning("No matching ID column found in TaxaData; replacing old table.")
    CrcBiomeScreenObject@TaxaData <- processed_taxa
    return(CrcBiomeScreenObject)
  }

  # Use first available matching ID column
  merge_col <- merge_col[1]

  # Align names for join
  names(old_taxa)[names(old_taxa) == merge_col] <- "SequenceID"

  # Combine tables
  merged <- dplyr::left_join(old_taxa, processed_taxa, by = "SequenceID")

  CrcBiomeScreenObject@TaxaData <- merged
  return(CrcBiomeScreenObject)
}

