test_that("NormalizeData works correctly", {
  # Construct test input data
  rel_abund <- matrix(runif(9), nrow = 3, dimnames = list(
    c(
      "k__Bacteria|p__Firmicutes|c__Clostridia|o__order1|f__family1|g__test",
      "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__order2|f__family2|g__test1",
      "k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__order3|f__family3|g__test2"
    ),
    c("sample1", "sample2", "sample3")
  ))

  taxa_data <- rownames(rel_abund)

  sample_data <- data.frame(
    sample_id = c("sample1", "sample2", "sample3"),
    number_reads = c(100000, 120000, 110000),
    row.names = c("sample1", "sample2", "sample3")
  )

  # Create the object
  crc_obj <- CreateCrcBiomeScreenObject(
    RelativeAbundance = rel_abund,
    TaxaData = taxa_data,
    SampleData = sample_data
  )

  # Now test SplitTaxas
  split_result <- SplitTaxas(crc_obj)
  # Keep the Genus level
  genus_result <- KeepGenusLevel(split_result)
  # Normalize the data
  normalized_GMPR <- NormalizeData(genus_result, method = "GMPR")
  normalized_TSS <- NormalizeData(genus_result, method = "TSS")
  # Check result format
  expect_s3_class(normalized_GMPR$NormalizedData, "data.frame")
  expect_s3_class(normalized_TSS$NormalizedData, "data.frame")
  expect_equal(attributes(normalized_GMPR$NormalizedData)$NormalizationMethod, "GMPR")
  expect_equal(attributes(normalized_TSS$NormalizedData)$NormalizationMethod, "TSS")
})
