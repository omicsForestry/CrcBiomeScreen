



SplitTaxas <- function(CrcBiomeScreenObject) {
    CrcBiomeScreenObject$TaxaData <- 
        CrcBiomeScreenObject$TaxaData %>% 
        tibble(variable = CrcBiomeScreenObject$TaxaData) %>%
        separate(variable, into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
                sep = '\\|', fill = "right") %>%
        mutate(across(everything(), ~ ifelse(. == "", NA, sub("^[a-z]__", "", .)))) %>%
        as.data.frame()
    return(CrcBiomeScreenObject)
}