test_that("checkClassBalance works correctly", {
  skip_on_bioc()
  library(curatedMetagenomicData)
  toydata <- curatedMetagenomicData(
    "ThomasAM_2018a.relative_abundance",
    dryrun = FALSE, rownames = "short"
  )

  toydata_object <- CreateCrcBiomeScreenObject(
    RelativeAbundance = toydata[[1]]@assays@data@listData$relative_abundance,
    TaxaData = toydata[[1]]@rowLinks$nodeLab,
    SampleData = toydata[[1]]@colData
  )
  toydata_object <- SplitTaxas(toydata_object)
  toydata_object <- KeepTaxonomicLevel(toydata_object, level = "Genus")

  toydata_object <- NormalizeData(toydata_object, method = "GMPR", level = "Genus")
  k <- 0.6
  toydata_object <- SplitDataSet(toydata_object, label = c("control", "CRC"), partition = k)
  checkClassBalanceResult <- checkClassBalance(getModelData(toydata_object)$TrainLabel)

  pdf_name <- "class_balance_plot.pdf"
  # Check result format
  expect_equal(class(checkClassBalanceResult), "list")
  expect_true(file.exists(pdf_name))
})
