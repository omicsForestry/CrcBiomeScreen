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
