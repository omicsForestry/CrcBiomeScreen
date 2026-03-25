test_that("CreateCrcBiomeScreenObject works correctly", {
  rel_abund <- data.frame(
    S1 = c(10, 20, 70),
    S2 = c(30, 30, 40)
  )
  rownames(rel_abund) <- c("TaxaA", "TaxaB", "TaxaC")

  taxa_info <- data.frame(
    Taxa = rownames(rel_abund),
    stringsAsFactors = FALSE
  )

  sample_info <- data.frame(
    number_reads = c(10000, 12000),
    study_condition = c("control", "CRC"),
    row.names = c("S1", "S2"),
    stringsAsFactors = FALSE
  )

  obj <- CreateCrcBiomeScreenObject(
    RelativeAbundance = rel_abund,
    TaxaData = taxa_info,
    SampleData = sample_info
  )

  expect_s4_class(obj, "CrcBiomeScreen")
  expect_true(nrow(getTaxaData(obj)) == 3)
  expect_true(nrow(getSampleData(obj)) == 2)
})
