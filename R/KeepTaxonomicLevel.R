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
#' 
#' 
KeepTaxonomicLevel <- function(CrcBiomeScreenObject, level = "Genus") {
  taxa_df <- as.data.frame(CrcBiomeScreenObject$TaxaData, stringsAsFactors = FALSE)
  valid_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (!level %in% valid_levels) stop("Invalid level.")
  
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
  
  # For each taxa (a row in TaxaData), select the group name: first try the target level, if NA, fall back to the previous level, and finally fall back to OriginalTaxa
  pick_group <- function(row) {
    idx <- match(level, valid_levels)
    for (i in idx:1) {
      val <- row[[ valid_levels[i] ]]
      if (!is.na(val) && nzchar(as.character(val))) return(as.character(val))
    }
    
    if (!is.na(row[["OriginalTaxa"]]) && nzchar(as.character(row[["OriginalTaxa"]]))) {
      return(as.character(row[["OriginalTaxa"]]))
    }
    return(NA_character_)
  }
  
  group_vec <- apply(taxa_df[, valid_levels], 1, pick_group)
  # Remove the _repeated suffixes
  group_vec <- gsub("(_uncultured)+$", "_uncultured", group_vec)
  group_vec <- gsub("(_unclassified)+$", "_unclassified", group_vec)
  group_vec <- gsub("(_unknown)+$", "_unknown", group_vec)
  
  ab <- as.matrix(CrcBiomeScreenObject$AbsoluteAbundance)
  if (nrow(ab) != length(group_vec)) stop("Column count of abundance does not match number of taxa entries.")
  
  grouped <- rowsum(ab, group = group_vec, na.rm = TRUE)  # result: rows = groups, cols = samples  
  LevelData <- as.data.frame(LevelData, stringsAsFactors = FALSE)

  CrcBiomeScreenObject$TaxaLevelData[[paste0(level, "LevelData")]] <- LevelData
  return(CrcBiomeScreenObject)
}



# KeepTaxonomicLevel <- function(CrcBiomeScreenObject, level = "Genus") {
#   # Ensure the user-provided level is valid
#   valid_levels <- colnames(CrcBiomeScreenObject$TaxaData)
#   if (!level %in% valid_levels) {
#     stop(paste("Invalid taxonomic level provided. Please choose from:", 
#                paste(valid_levels, collapse = ", ")))
#   }
#   
#   level_sym <- rlang::sym(level)
#   
#   LevelData <-
#     CrcBiomeScreenObject$AbsoluteAbundance %>%
#     as.data.frame() %>%
#     dplyr::mutate(tax_level = CrcBiomeScreenObject$TaxaData[[level]]) %>%
#     rstatix::group_by(tax_level) %>%
#     dplyr::summarise_all(sum) %>%
#     tibble::column_to_rownames("tax_level") %>%
#     t() %>%
#     as.data.frame()
#   
#   # Remove _uncultured/_unclassified/unknown repeatedly at the end of the names
#   rownames(LevelData) <- gsub("(_uncultured)+$", "_uncultured", rownames(LevelData))
#   rownames(LevelData) <- gsub("(_unclassified)+$", "_unclassified", rownames(LevelData))
#   rownames(LevelData) <- gsub("(_unknown)+$", "_unknown", rownames(LevelData))
#   
#   CrcBiomeScreenObject$TaxaLevelData[[paste0(level, "LevelData")]] <- LevelData
#   
#   return(CrcBiomeScreenObject)
# }


# KeepTaxonomicLevel <- function(CrcBiomeScreenObject, level = "Genus") {
#   # Ensure the user-provided level is valid
#   valid_levels <- colnames(CrcBiomeScreenObject$TaxaData)
#   if (!level %in% valid_levels) {
#     stop(paste("Invalid taxonomic level provided. Please choose from:", paste(valid_levels, collapse = ", ")))
#   }
# 
#   # Use rlang::sym() and !! to convert the string 'level' into a variable name
#   # this allows dplyr to correctly use the column name specified by the user
#   level_sym <- rlang::sym(level)
# 
#   LevelData <-
#     CrcBiomeScreenObject$AbsoluteAbundance %>%
#     as.data.frame() %>%
#     # Use !!level_sym to dynamically select the taxonomic column
#     dplyr::mutate(tax_level = CrcBiomeScreenObject$TaxaData[[level]]) %>%
#     rstatix::group_by(tax_level) %>%
#     dplyr::summarise_all(sum) %>%
#     tibble::column_to_rownames("tax_level") %>%
#     t() %>%
#     as.data.frame()
# 
#   # Change the key name in the returned object to make it more generic
#   CrcBiomeScreenObject$TaxaLevelData[[paste0(level, "LevelData")]] <- LevelData
# 
#   return(CrcBiomeScreenObject)
# }

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