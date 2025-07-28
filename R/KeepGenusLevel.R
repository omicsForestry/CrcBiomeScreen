#' Keep the taxa data at genus level
#'
#' @param CrcBiomeScreenObject
#'
#' @importFrom dplyr mutate across summarise_all
#' @importFrom tidyr separate %>%
#' @importFrom tibble column_to_rownames
#' @importFrom rstatix group_by
#' @importFrom dplyr
#'
#' @return CrcBiomeScreenObject$GenusLevelData
#' @export
#'
#' @examples CrcBiomeScreenObject <- KeepGenusLevel(CrcBiomeScreenObject)
KeepGenusLevel <- function(CrcBiomeScreenObject) {
  CrcBiomeScreenObject$GenusLevelData <-
    CrcBiomeScreenObject$AbsoluteAbundance %>%
    as.data.frame() %>%
    mutate(genus = CrcBiomeScreenObject$TaxaData$Genus) %>%
    group_by(genus) %>%
    summarise_all(sum) %>%
    column_to_rownames("genus") %>%
    t() %>%
    as.data.frame()
  return(CrcBiomeScreenObject)
}
