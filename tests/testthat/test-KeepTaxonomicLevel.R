test_that("SplitTaxas and KeepTaxonomicLevel handle uncultured/unclassified/unknown correctly", {
  taxa_vec <- c(
    "D_0__Bacteria.D_1__Firmicutes.D_2__Clostridia.D_3__uncultured.D_4__uncultured",
    "D_0__Bacteria.D_1__Proteobacteria.D_2__Gammaproteobacteria.D_3__Enterobacterales.D_4__Enterobacteriaceae.D_5__unclassified"
  )

  ab <- data.frame(
    t1 = c(10, 5, 30, 25),
    t2 = c(20, 15, 40, 35)
  )
  colnames(ab) <- taxa_vec

  obj <- CreateCrcBiomeScreenObject(
    TaxaData = taxa_vec,
    AbsoluteAbundance = t(ab)
  )

  obj2 <- SplitTaxas(obj)
  obj3 <- KeepTaxonomicLevel(obj2, level = "Genus")

  # print to inspect
  print(colnames(getTaxaData(obj2)))
  print(rownames(obj3@TaxaLevelData$GenusLevelData))
  expect_equal(rownames(obj3@TaxaLevelData$GenusLevelData), c(
    "Clostridia_uncultured",
    "Enterobacteriaceae_unclassified"
  ))
})
