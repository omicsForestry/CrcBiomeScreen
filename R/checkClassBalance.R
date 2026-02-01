#' Check the sample distribution of the dataset and
#' give the suggestion if need the class weight or not
#'
#' @param labels The label for distribution
#' @param threshold The threshold for the ratio (0.5) if it is the
#' imbalanced dataset
#' @param plot Choose to have the figures or not
#' @param outdir The output directory where plots will be saved
#' (default: tempdir()).
#'
#' @importFrom tibble tibble
#'
#' @return A \linkS4class{CrcBiomeScreen} object with updated slots.
#' @export
#'
#' @examples
#' # Small toy example for runnable demonstration
#'
#' # Fake labels for the training set
#' train_labels <- factor(c("control", "CRC", "control", "CRC"))
#'
#' # Create a tiny relative abundance matrix (1 taxa × 4 samples)
#' rel_abund <- data.frame(
#'   S1 = 10, S2 = 20, S3 = 30, S4 = 40
#' )
#' rownames(rel_abund) <- "TaxaA"
#'
#' # Minimal taxonomy and sample metadata
#' taxa_info <- data.frame(Taxa = "TaxaA", stringsAsFactors = FALSE)
#' sample_info <- data.frame(
#'   number_reads = c(10000, 10000, 12000, 12000),
#'   condition     = train_labels,
#'   row.names     = paste0("S", 1:4),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Build CrcBiomeScreen object
#' toy_obj <- CreateCrcBiomeScreenObject(
#'   RelativeAbundance = rel_abund,
#'   TaxaData          = taxa_info,
#'   SampleData        = sample_info
#' )
#'
#' # Add minimal ModelData needed for class balance check
#' toy_obj@ModelData <- list(
#'   TrainLabel = train_labels
#' )
#'
#' # Check class balance
#' checkClassBalance(toy_obj@ModelData$TrainLabel)
checkClassBalance <- function(labels,
                              outdir = tempdir(),
                              threshold = 0.5,
                              plot = TRUE) {
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
    pdf(file.path(outdir, "class_balance_plot.pdf"), width = 8, height = 6)
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
  message("Class counts:\n")
  message(class_counts)
  message("Class proportions:\n")
  message(round(class_props, 3))
  message("Balance status:\n", suggestion, "\n")

    # Return result
  return(list(
    class_counts = class_counts,
    class_proportions = class_props,
    is_imbalanced = is_imbalanced,
    suggestion = suggestion
  ))
}
