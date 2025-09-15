#' Split and clean taxonomy strings
#'
#' This function automatically detects the taxonomy string format (e.g., MetaPhlAn, QIIME, SILVA, GTDB),
#' splits the string into standard taxonomic ranks (Kingdom to Species), 
#' retains the original taxonomy string in a new column (`OriginalTaxa`), 
#' and refines labels such as "uncultured" or "unclassified" by appending the parent rank.
#'
#' @param CrcBiomeScreenObject
#'
#' @importFrom dplyr mutate across
#' @importFrom tidyr separate
#' @importFrom tibble tibble
#'
#' @return CrcBiomeScreenObject$TaxaData
#' @export
#'
#' @examples CrcBiomeScreenObject <- SplitTaxas(CrcBiomeScreenObject)
# SplitTaxas <- function(CrcBiomeScreenObject) {
#   CrcBiomeScreenObject$TaxaData <-
#     CrcBiomeScreenObject$TaxaData %>%
#     tibble(variable = CrcBiomeScreenObject$TaxaData) %>%
#     separate(variable,
#       into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
#       sep = "\\|", fill = "right"
#     ) %>%
#     mutate(across(everything(), ~ ifelse(. == "", NA, sub("^[a-z]__", "", .)))) %>%
#     as.data.frame()
#   return(CrcBiomeScreenObject)
# }
# 

SplitTaxas <- function(CrcBiomeScreenObject) {
  taxa_vec <- CrcBiomeScreenObject$TaxaData
  
  # 1. Detect separator / style
  if (any(grepl("D_\\d+__", taxa_vec))) {
    sep <- "\\."
    style <- "qiime"
  } else if (any(grepl("\\|", taxa_vec))) {
    sep <- "\\|"
    style <- "metaphlan"
  } else if (any(grepl(";", taxa_vec))) {
    sep <- ";"
    style <- "qiime1"
  } else if (any(grepl("_", taxa_vec))) {
    sep <- "_"
    style <- "silva"
  } else {
    stop("❌ Cannot detect taxonomy separator. Please check input format.")
  }
  
  # 2. Split
  df <- tibble(OriginalTaxa = taxa_vec) %>%
    separate(
      OriginalTaxa,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = sep,
      fill = "right"
    )
  
  # 3. Clean names
  if (style == "qiime") {
    df <- df %>%
      mutate(across(everything(), ~ ifelse(. %in% c("", "__"), NA, sub("^D_\\d+__", "", .))))
  } else {
    df <- df %>%
      mutate(across(everything(), ~ ifelse(. %in% c("", "__"), NA, sub("^[a-z]__*", "", .))))
  }
  
  # 4. Process the uncultured / unclassified
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (i in seq_along(tax_levels)[-1]) {
    lvl <- tax_levels[i]
    parent <- tax_levels[i - 1]
    df[[lvl]] <- ifelse(
      df[[lvl]] %in% c("uncultured", "unclassified"),
      paste0(df[[parent]], "_", df[[lvl]]),
      df[[lvl]]
    )
  }
  
  # 5. Return to object
  CrcBiomeScreenObject$TaxaData <- as.data.frame(df)
  rownames(CrcBiomeScreenObject$TaxaData) <- taxa_vec
  return(CrcBiomeScreenObject)
}
