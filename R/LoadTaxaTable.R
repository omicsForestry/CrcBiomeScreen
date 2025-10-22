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
#' @return The CrcBiomeScreenObject with the loaded taxa table.
#' @export
#'
#' @examples
#' CrcBiomeScreenObject <- LoadTaxaTable(CrcBiomeScreenObject,
#'                                        taxa_table = my_taxa_table,
#'                                        id_column = "ASV_ID",
#'                                        taxa_column = "Taxonomy")

LoadTaxaTable <- function(CrcBiomeScreenObject, taxa_table, id_column, taxa_column) {

  # --- Input Validation ---
  if (!inherits(CrcBiomeScreenObject, "CrcBiomeScreenObject")) {
    stop("Input 'CrcBiomeScreenObject' must be of class 'CrcBiomeScreenObject'.")
  }
  if (!is.data.frame(taxa_table)) {
    stop("Input 'taxa_table' must be a data frame.")
  }
  if (!id_column %in% colnames(taxa_table)) {
    stop(sprintf("The specified 'id_column' ('", id_column, "') was not found in 'taxa_table'."))
  }
  if (!taxa_column %in% colnames(taxa_table)) {
    stop(sprintf("The specified 'taxa_column' ('", taxa_column, "') was not found in 'taxa_table'."))
  }

  # --- Data Processing ---
  # Select only the necessary columns and rename them for consistency
  processed_taxa <- taxa_table %>%
    dplyr::select(!!rlang::sym(id_column), !!rlang::sym(taxa_column)) %>%
    dplyr::rename(SequenceID = !!rlang::sym(id_column),
                  TaxaString = !!rlang::sym(taxa_column))

  # Store the processed taxa table in the object
  CrcBiomeScreenObject$TaxaDataTable <- processed_taxa

  # Optionally, you might want to parse these TaxaStrings immediately.
  # For now, we'll store them raw and let SplitTaxas handle parsing if needed.
  # If you want to parse them here, you would call SplitTaxas(CrcBiomeScreenObject, taxa_col_name = "TaxaString", ...).
  # For this function's purpose, just loading is enough.

  # --- Return Updated Object ---
  return(CrcBiomeScreenObject)
}
