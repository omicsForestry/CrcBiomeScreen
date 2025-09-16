#' Keep a specific taxonomic level
#'
#' This function aggregates abundance data to a specified taxonomic level.
#'
#' @param CrcBiomeScreenObject A list or object containing 'AbsoluteAbundance' and 'TaxaData'.
#' @param level The taxonomic level to aggregate to (e.g., "Family", "Genus", "Species").
#'
#' @return The CrcBiomeScreenObject with a new data frame aggregated at the specified level.
#' @export
KeepTaxonomicLevel <- function(CrcBiomeScreenObject, level = "Genus") {
  # Ensure the user-provided level is valid
  valid_levels <- colnames(CrcBiomeScreenObject$TaxaData)
  if (!level %in% valid_levels) {
    stop(paste("Invalid taxonomic level provided. Please choose from:", paste(valid_levels, collapse = ", ")))
  }
  
  # Use rlang::sym() and !! to convert the string 'level' into a variable name
  # this allows dplyr to correctly use the column name specified by the user
  level_sym <- rlang::sym(level)
  
  LevelData <-
    CrcBiomeScreenObject$AbsoluteAbundance %>%
    as.data.frame() %>%
    # Use !!level_sym to dynamically select the taxonomic column
    dplyr::mutate(tax_level = CrcBiomeScreenObject$TaxaData[[level]]) %>%
    rstatix::group_by(tax_level) %>%
    dplyr::summarise_all(sum) %>%
    tibble::column_to_rownames("tax_level") %>%
    t() %>%
    as.data.frame()
  
  # Change the key name in the returned object to make it more generic
  CrcBiomeScreenObject$TaxaLevelData[[paste0(level, "LevelData")]] <- LevelData
  
  return(CrcBiomeScreenObject)
}

# KeepTaxonomicLevel <- function(CrcBiomeScreenObject) {
#   
#   CrcBiomeScreenObject$GenusLevelData <-
#     
#     CrcBiomeScreenObject$AbsoluteAbundance %>%
#     
#     as.data.frame() %>%
#     
#     mutate(genus = CrcBiomeScreenObject$TaxaData$Genus) %>%
#     
#     rstatix::group_by(genus) %>%
#     
#     summarise_all(sum) %>%
#     
#     column_to_rownames("genus") %>%
#     
#     t() %>%
#     
#     as.data.frame()
#   
#   return(CrcBiomeScreenObject)
#   
# }