test_that("SplitTaxas works correctly", {
  # Construct test input data
  rel_abund <- matrix(runif(9), nrow = 3, dimnames = list(
    c(
      "k__Bacteria|p__Firmicutes|c__Clostridia",
      "k__Bacteria|p__Bacteroidetes|c__Bacteroidia",
      "k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria"
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

  # Check result format
  expect_s3_class(split_result$TaxaData, "data.frame")
  expect_true(all(c("Kingdom", "Phylum", "Class") %in% colnames(split_result$TaxaData)))
})
