test_that("SplitTaxas and KeepTaxonomicLevel handle uncultured/unclassified/unknown correctly", {
  taxa_vec <- c(
    "D_0__Bacteria.D_1__Firmicutes.D_2__Clostridia.D_3__uncultured.D_4__uncultured",
    "D_0__Bacteria.D_1__Proteobacteria.D_2__Gammaproteobacteria.D_3__Enterobacterales.D_4__Enterobacteriaceae.D_5__unclassified",
    "D_0__Bacteria.D_1__Actinobacteria.D_2__Actinobacteria.D_3__unknown"
  )
  
  ab <- data.frame(
    t1 = c(10, 5),
    t2 = c(20, 15),
    t3 = c(30, 25)
  )
  rownames(ab) <- c("Sample1","Sample2")
  
  CrcBiomeScreenObject <- list(
    TaxaData = taxa_vec,            # order must match ab columns (t1,t2,t3)
    AbsoluteAbundance = ab,
    TaxaLevelData = list()
  )

  CrcBiomeScreenObject <- SplitTaxas(CrcBiomeScreenObject)
  CrcBiomeScreenObject <- KeepTaxonomicLevel(CrcBiomeScreenObject, level = "Genus")

  expect_true("OriginalTaxa" %in% colnames(CrcBiomeScreenObject$TaxaData))
  
  genus_data <- CrcBiomeScreenObject$TaxaLevelData$GenusLevelData
  
  expect_true("Clostridia_uncultured" %in% rownames(genus_data))
  expect_true("Enterobacteriaceae_unclassified" %in% rownames(genus_data))
  expect_true("Actinobacteria_unknown" %in% rownames(genus_data))
  
  expect_equal(genus_data["Clostridia_uncultured", "Sample1"], 10)
  expect_equal(genus_data["Enterobacteriaceae_unclassified", "Sample2"], 15)
  expect_equal(genus_data["Actinobacteria_unknown", "Sample1"], 30)
})
