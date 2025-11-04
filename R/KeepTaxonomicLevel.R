#' Keep a specific taxonomic level
#'
#' This function aggregates abundance data to a specified taxonomic level.
#'
#' @param CrcBiomeScreenObject A list or object containing 'AbsoluteAbundance' and 'TaxaData'.
#' @param level The taxonomic level to aggregate to (e.g., "Family", "Genus", "Species").
#'
#' @importFrom magrittr %>%
#'
#' @return The CrcBiomeScreenObject with a new data frame aggregated at the specified level.
#' @export
#' @examples
#' # Keep only the genus level data
#' CrcBiomeScreenObject <- KeepTaxonomicLevel(CrcBiomeScreenObject,level = "Genus")
#' @title Summarize abundance data at a given taxonomic level
#' @description
#' Aggregate absolute abundance data in a \linkS4class{CrcBiomeScreen} object
#' to a specified taxonomic level (e.g. "Genus" or "Family").
#'
#' @param CrcBiomeScreenObject A \linkS4class{CrcBiomeScreen} object.
#' @param level Taxonomic level to summarize to. One of
#' "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species".
#'
#' @return The same \linkS4class{CrcBiomeScreen} object, updated with a slot
#' \code{@TaxaLevelData} a new data frame in
#' \code{@GenusLevelData} (or the corresponding level).
#' @export
KeepTaxonomicLevel <- function(CrcBiomeScreenObject, level = "Genus") {

  # Extract and validate data
  taxa_df <- as.data.frame(getTaxaData(CrcBiomeScreenObject), stringsAsFactors = FALSE)
  ab_df   <- as.data.frame(getAbsoluteAbundance(CrcBiomeScreenObject), stringsAsFactors = FALSE)

  valid_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (!level %in% valid_levels)
    stop("Invalid level. Choose one of: ", paste(valid_levels, collapse = ", "))

  if (!level %in% colnames(taxa_df))
    stop(paste0("Column '", level, "' not found in TaxaData."))

  # Define helper: pick most appropriate grouping name
  pick_group <- function(row) {
    idx <- match(level, valid_levels)
    for (i in idx:1) {
      val <- row[[ valid_levels[i] ]]
      if (!is.na(val) && nzchar(as.character(val)))
        return(as.character(val))
    }
    # fallback to row name if all missing
    if (!is.na(row[["OriginalTaxa"]]) && nzchar(as.character(row[["OriginalTaxa"]])))
      return(as.character(row[["OriginalTaxa"]]))
    return(NA_character_)
  }

  # Construct grouping vector
  group_vec <- apply(taxa_df[, valid_levels[valid_levels %in% colnames(taxa_df)]], 1, pick_group)

  # cleanup repeated suffixes
  group_vec <- gsub("(_uncultured)+$", "_uncultured", group_vec)
  group_vec <- gsub("(_unclassified)+$", "_unclassified", group_vec)
  group_vec <- gsub("(_unknown)+$", "_unknown", group_vec)

  if (nrow(ab_df) != length(group_vec)) {
    stop("Row count of AbsoluteAbundance does not match number of taxa entries.")
  }

  # Aggregate by selected level
  grouped <- rowsum(as.matrix(ab_df), group = group_vec, na.rm = TRUE)
  grouped <- as.data.frame(grouped, stringsAsFactors = FALSE)

  # Store summarized data into the corresponding slot
  # we store under a unified slot (e.g., GenusLevelData) for simplicity
  CrcBiomeScreenObject@TaxaLevelData <- list()
  CrcBiomeScreenObject@TaxaLevelData[[paste0(level, "LevelData")]] <- grouped

  return(CrcBiomeScreenObject)
}
