checkClassBalance <- function(labels, threshold = 0.3, plot = TRUE) {
  # Convert to factor
  labels <- as.factor(labels)
  class_counts <- table(labels)
  class_props <- prop.table(class_counts)

  # Identify minority classes
  is_imbalanced <- any(class_props < threshold | class_props > (1 - threshold))
  suggestion <- if (is_imbalanced) {
    "Unbalanced: Using class weights is recommended."
  } else {
    "Balanced: Class weights are likely unnecessary."
  }

  # Optional plot
  if (plot) {
    pdf("class_balance_plot.pdf", width = 8, height = 6)
    barplot(class_counts,
      col = "steelblue",
      main = "Class Distribution",
      ylab = "Sample Count",
      xlab = "Class Label"
    )
    abline(h = mean(class_counts), col = "red", lty = 2)
    dev.off()
  }

  # Print summary
  cat("Class counts:\n")
  print(class_counts)
  cat("Class proportions:\n")
  print(round(class_props, 3))
  cat("Balance status:\n", suggestion, "\n")

  # Return result
  return(list(
    class_counts = class_counts,
    class_proportions = class_props,
    is_imbalanced = is_imbalanced,
    suggestion = suggestion
  ))
}
