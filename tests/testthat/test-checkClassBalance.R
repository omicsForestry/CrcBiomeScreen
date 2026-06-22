test_that("checkClassBalance returns expected structure", {

  # minimal labels (2 classes, imbalanced)
  train_labels <- factor(c(rep("control", 8), rep("CRC", 2)))

  res <- checkClassBalance(train_labels, plot = TRUE)

  expect_type(res, "list")
  expect_true(all(c("class_counts","class_proportions","is_imbalanced","suggestion") %in% names(res)))
  expect_true(setequal(names(res$class_counts), levels(train_labels)))
  expect_equal(sum(res$class_counts), length(train_labels))
  expect_equal(sum(res$class_proportions), 1)
})
