test_that("SplitTaxas works correctly", {
  set.seed(123)
  library(CrcBiomeScreen)
  toy_taxa_strings <- c(
    "D_0__Bacteria|D_1__Firmicutes|D_2__Clostridia|D_3__OrderA|D_4__FamilyA|D_5__GenusA",
    "D_0__Bacteria|D_1__Firmicutes|D_2__Clostridia|D_3__OrderB|D_4__FamilyB|D_5__GenusB",
    "D_0__Bacteria|D_1__Firmicutes|D_2__Clostridia|D_3__OrderC|D_4__FamilyC|D_5__GenusC",
    "D_0__Bacteria|D_1__Bacteroidetes|D_2__Bacteroidia|D_3__OrderD|D_4__FamilyD|D_5__GenusD",
    "D_0__Bacteria|D_1__Bacteroidetes|D_2__Bacteroidia|D_3__OrderE|D_4__FamilyE|D_5__GenusE",
    "D_0__Bacteria|D_1__Proteobacteria|D_2__Gammaproteobacteria|D_3__OrderF|D_4__FamilyF|D_5__GenusF"
  )

  toy_taxa <- data.frame(
    Taxa = toy_taxa_strings,
    stringsAsFactors = FALSE
  )
  train_samples <- paste0("S", 1:12)

  toy_abs <- matrix(
    c(
      rpois(6 * 6, lambda = 54.8887777),   # controls
      rpois(6 * 6, lambda = 55)    # CRC
    ),
    nrow = 6,
    ncol = 12
  )

  rownames(toy_abs) <- toy_taxa_strings
  colnames(toy_abs) <- train_samples
  toy_abs <- as.data.frame(toy_abs)

  toy_sample <- data.frame(
    number_reads = rep(10000, 12),
    study_condition = c(rep("control", 6), rep("CRC", 6)),
    row.names = train_samples,
    stringsAsFactors = FALSE
  )
  toy_obj <- CreateCrcBiomeScreenObject(
    AbsoluteAbundance = toy_abs,
    TaxaData = toy_taxa,
    SampleData = toy_sample
  )

  # Now test SplitTaxas
  split_result <- SplitTaxas(toy_obj)

  # Check result format
  expect_s3_class(getTaxaData(split_result), "data.frame")
  expect_true(all(c("Kingdom", "Phylum", "Class") %in% colnames(getTaxaData(split_result))))
})
