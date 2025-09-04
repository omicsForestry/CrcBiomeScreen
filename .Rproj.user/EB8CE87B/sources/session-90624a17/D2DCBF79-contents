test_that("NormalizeData works correctly", {

  toydata <- curatedMetagenomicData(
    "ThomasAM_2018a.relative_abundance"
    , dryrun = FALSE, rownames = "short")

  toydata_object <- CreateCrcBiomeScreenObject(RelativeAbundance = toydata[[1]]@assays@data@listData$relative_abundance,
                                                     TaxaData = toydata[[1]]@rowLinks$nodeLab,
                                                     SampleData = toydata[[1]]@colData)
  toydata_object <- SplitTaxas(toydata_object)
  toydata_object <- KeepGenusLevel(toydata_object)

  toydata_object <- NormalizeData(toydata_object, method = "GMPR")
  k = 0.6
  toydata_object <- SplitDataSet(toydata_object, label = c("control","CRC"), partition = k)

  # Check result format
  expect_equal(class(toydata_object$ModelData), "list")
  expect_s3_class(toydata_object$ModelData[["Training"]], "data.frame")
  expect_s3_class(toydata_object$ModelData[["Test"]], "data.frame")
  expect_equal(class(toydata_object$ModelData[["TrainLabel"]]), "character")
  expect_equal(class(toydata_object$ModelData[["TestLabel"]]), "character")
  expect_equal(attributes(toydata_object$ModelData)$names,c("Training","Test","TrainLabel","TestLabel"))
  expect_equal(attributes(toydata_object$ModelData)$`Split Partition`,k)

})
