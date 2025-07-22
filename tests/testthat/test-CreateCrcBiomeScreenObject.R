test_that("CreateCrcBiomeScreenObject works correctly", {
  toydata <- curatedMetagenomicData(
    "ThomasAM_2018a.relative_abundance",
    dryrun = FALSE, rownames = "short"
  )

  rel_abund <- toydata[[1]]@assays@data@listData$relative_abundance
  taxa_data <- toydata[[1]]@rowLinks$nodeLab
  sample_data <- toydata[[1]]@colData

  result <- CreateCrcBiomeScreenObject(
    RelativeAbundance = rel_abund,
    TaxaData = taxa_data,
    SampleData = sample_data
  )

  expect_s3_class(result, "CrcBiomeScreenObject")

  expect_error(CreateCrcBiomeScreenObject(
    RelativeAbundance = NULL,
    TaxaData = NULL,
    SampleData = sample_data
  ))
})
