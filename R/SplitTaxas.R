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
#' @importFrom magrittr %>%
#'
#' @return CrcBiomeScreenObject$TaxaData
#' @export
#'
#' @examples CrcBiomeScreenObject <- SplitTaxas(CrcBiomeScreenObject)

SplitTaxas <- function(CrcBiomeScreenObject) {
  taxa_vec <- as.character(CrcBiomeScreenObject$TaxaData)

  # detect QIIME D_ style
  if (any(grepl("D_\\d+__", taxa_vec))) {
    taxa_vec2 <- gsub("\\.(?=__|D_\\d+__)", "|", taxa_vec, perl = TRUE)
    sep <- "\\|"
    style <- "qiime"
  } else if (any(grepl("\\|", taxa_vec))) {
    taxa_vec2 <- taxa_vec; sep <- "\\|"; style <- "metaphlan"
  } else if (any(grepl(";", taxa_vec))) {
    taxa_vec2 <- taxa_vec; sep <- ";"; style <- "semi"
  } else if (any(grepl("__", taxa_vec))) {
    taxa_vec2 <- taxa_vec; sep <- "\\."; style <- "dot"
  } else if (any(grepl("_", taxa_vec))) {
    taxa_vec2 <- taxa_vec; sep <- "_"; style <- "underscore"
  } else {
    taxa_vec2 <- taxa_vec; sep <- "\\|"; style <- "fallback"
  }
  
  # split into ranks
  suppressWarnings(taxa_df <- tibble::tibble(OriginalTaxa = taxa_vec, .rows = length(taxa_vec)) %>%
    dplyr::mutate(tmp = taxa_vec2) %>%
    tidyr::separate(tmp,
                    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                    sep = sep, fill = "right", remove = TRUE))
  
  # cleanup prefixes
  remove_prefix <- function(x) {
    if (is.na(x) || x == "" || x == "__") return(NA_character_)
    x <- trimws(x)
    x <- sub("^D_\\d+__", "", x)      # D_0__ style
    x <- sub("^[a-zA-Z]{1,2}__", "", x) # k__, p__, g__ ...
    x <- sub("^[a-zA-Z]_", "", x)     # g_ style
    x <- trimws(x)
    if (x == "" || x == "__") return(NA_character_)
    return(x)
  }
  taxa_df <- taxa_df %>% mutate(across(Kingdom:Species, ~ sapply(., remove_prefix, USE.NAMES = FALSE)))
  
  # handle bad labels: attach to parent but avoid duplicate suffixes
  bad_labels <- c("uncultured","unclassified","unknown")
  clean_parent <- function(x) {
    if (is.na(x)) return(NA_character_)
    # remove existing suffixes for parent candidate
    x <- gsub("(_uncultured)+$", "_uncultured", x)
    x <- gsub("(_unclassified)+$", "_unclassified", x)
    x <- gsub("(_unknown)+$", "_unknown", x)
    # if parent itself is NA, return NA
    if (is.na(x) || x == "") return(NA_character_)
    # strip trailing suffix for building new suffix (so we don't get a_b_uncultured_uncultured)
    x <- sub("(_uncultured|_unclassified|_unknown)$", "", x)
    return(x)
  }
  
  for (i in 2:length(c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))) {
    lvl <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")[i]
    parent <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")[i-1]
    cur <- taxa_df[[lvl]]
    par <- taxa_df[[parent]]
    replace_idx <- which(!is.na(cur) & tolower(cur) %in% bad_labels & !is.na(par))
    if (length(replace_idx) > 0) {
      for (ri in replace_idx) {
        pclean <- clean_parent(par[ri])
        if (!is.na(pclean)) {
          taxa_df[[lvl]][ri] <- paste0(pclean, "_", tolower(cur[ri]))
        } else {
          taxa_df[[lvl]][ri] <- tolower(cur[ri])
        }
      }
    }
  }
  
  # final cleanup: collapse repeated suffixes like _uncultured_uncultured -> _uncultured
  taxa_df <- taxa_df %>%
    mutate(across(Kingdom:Species, ~ ifelse(is.na(.), NA_character_,
                                            gsub("(_uncultured)+$", "_uncultured",
                                                 gsub("(_unclassified)+$", "_unclassified",
                                                      gsub("(_unknown)+$", "_unknown", .))))))

  CrcBiomeScreenObject$TaxaData <- as.data.frame(taxa_df, stringsAsFactors = FALSE)
  rownames(CrcBiomeScreenObject$TaxaData) <- taxa_df$OriginalTaxa
  CrcBiomeScreenObject$TaxaData <- CrcBiomeScreenObject$TaxaData[,-1]
  return(CrcBiomeScreenObject)
}


